from collections import Counter
from collections.abc import Mapping
from typing import Any, Literal

from tacular import (
    ELEMENT_LOOKUP,
    FRAGMENT_ION_LOOKUP,
    ElementInfo,
    IonType,
    IonTypeProperty,
)

from ..constants import ELECTRON_MASS, ModType
from ..diagnostics import PeptacularError
from ..proforma_components import (
    ChargedFormula,
    GlobalChargeCarrier,
)
from .mod import Mods
from .positions import validate_position

__all__ = [
    "Fragment",
]

# Maps internal ion type value tuples to their neutral loss diff relative to "by" (the default internal fragment).
# None means no difference from "by". Derived from tacular's internal(F,B) = deltaF + deltaB
# offsets (tacular>=1.1.0, itself derived from mzPAF's own primary-ion formulas: a=b-CO,
# c=b+NH3, x=y+CO2-H2O, z=y-NH3) and written using mzPAF's neutral-loss conventions: strung
# together signed tokens (mzPAF section 4.5), an ordinal prefix for a repeated named atom
# (e.g. "-2H", not "H2"), and the canonical group names the spec requires when they apply
# (e.g. "NH3" per "do not write an ammonia loss (NH3) as H3N"; "HCONH2"/Formamide for the
# combined CO+NH3 magnitude).
_INTERNAL_MASS_DIFFS: dict[tuple[str, str], str | None] = {
    ("a", "x"): "-2H",
    ("b", "x"): "+CO-2H",
    ("c", "x"): "+CHNO",
    ("a", "y"): "-CO",
    ("b", "y"): None,
    ("c", "y"): "+NH3",
    ("a", "z"): "-HCONH2",
    ("b", "z"): "-NH3",
    ("c", "z"): None,
}

# Maps peptacular variant ion types to their canonical mzPAF series string.
_ION_TYPE_TO_MZPAF_SERIES: dict[IonType, str] = {
    IonType.W_VALINE: "w",
    IonType.D_VALINE: "d",
    IonType.WB_ISOLEUCINE: "wb",
    IonType.WB_THREONINE: "wb",
    IonType.WA_ISOLEUCINE: "wa",
    IonType.WA_THREONINE: "wa",
    IonType.DB_ISOLEUCINE: "db",
    IonType.DB_THREONINE: "db",
    IonType.DA_ISOLEUCINE: "da",
    IonType.DA_THREONINE: "da",
    # mzPAF 1.0.1 "z" is the z-dot radical (sum + H2O - NH2). The other z and c variants
    # have no letter of their own, so they are written as that series plus a hydrogen
    # delta (_MZPAF_SERIES_DELTA): Biemann z = "z-H", z+H = "z+H", c-H = "c-H".
    IonType.Z: "z",
    IonType.Z_RADICAL: "z",
    IonType.Z_PLUS_H: "z",
    IonType.C_MINUS_H: "c",
}

_MZPAF_SERIES_DELTA: dict[IonType, str] = {
    IonType.Z: "-H",
    IonType.Z_PLUS_H: "+H",
    IonType.C_MINUS_H: "-H",
}


class Fragment:
    """One theoretical ion: a fragment or precursor with its ion type, position, charge and mass.

    Returned by :meth:`ProFormaAnnotation.frag`, :meth:`ProFormaAnnotation.fragment` and
    :func:`peptacular.fragment`. ``mass`` is the mass of the charged ion (adducts included), so
    ``mz`` is ``mass / abs(charge_state)`` and ``neutral_mass`` removes the charge carriers.
    Use :meth:`to_mzpaf` for an mzPAF annotation string.

    >>> import peptacular as pt
    >>> frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2)
    >>> frag.ion_type, frag.position, frag.charge_state
    (<IonType.B: 'b'>, 2, 1)
    >>> round(frag.mz, 4)
    227.1026

    :param ion_type: Ion type (``b``, ``y``, precursor, immonium, ...).
    :type ion_type: IonType
    :param position: Ion ordinal (residue count) for terminal ions, a one-based inclusive
        ``(start, end)`` residue range for internal ions, or None for the precursor.
    :type position: int | tuple[int, int] | None
    :param mass: Mass of the charged ion in Da, charge carriers included.
    :type mass: float
    :param monoisotopic: True for a monoisotopic mass, False for an average mass.
    :type monoisotopic: bool
    :param charge_state: Signed total charge.
    :type charge_state: int
    :param charge_adducts: Charge carrier strings (e.g. ``("Na:z+1",)``); None means protonated.
    :type charge_adducts: tuple[str, ...] | None
    :param external_charge: Part of ``charge_state`` carried by external adducts rather than a
        charged modification. Defaults to ``charge_state``.
    :type external_charge: int | None
    :param isotopes: Isotope shifts as ``{isotope: count}``, or an int for a number of 13C.
    :type isotopes: Mapping[str, int] | int | None
    :param deltas: Neutral losses or gains as ``{formula_or_mass: count}``.
    :type deltas: Mapping[str | float, int] | None
    :param composition: Precomputed elemental composition, if known.
    :type composition: Counter[ElementInfo] | None
    :param parent_sequence: ProForma string of the parent peptidoform, used to compute
        ``composition`` and ``sequence`` lazily.
    :type parent_sequence: str | None
    :param parent_sequence_length: Residue count of the parent peptidoform.
    :type parent_sequence_length: int | None
    """

    def __init__(
        self,
        ion_type: IonType,
        position: int | tuple[int, int] | None,
        mass: float,
        monoisotopic: bool,
        charge_state: int,
        charge_adducts: tuple[str, ...] | None = None,
        external_charge: int | None = None,
        isotopes: Mapping[str, int] | int | None = None,
        deltas: Mapping[str | float, int] | None = None,
        composition: Counter[ElementInfo] | None = None,
        parent_sequence: str | None = None,
        parent_sequence_length: int | None = None,
    ) -> None:
        self.ion_type: IonType = ion_type
        self.position: int | tuple[int, int] | None = position
        self.mass: int | float = mass
        self.monoisotopic: bool = monoisotopic
        self.charge_state: int = charge_state
        # If None and charge_state != 0: means protonated
        self._charge_adducts: tuple[str, ...] | None = charge_adducts
        # The portion of charge_state that comes from real external adducts/charge carriers,
        # as opposed to charge intrinsic to an internal formula modification (e.g. [Formula:...:z+N]).
        # Used to reconstruct the default proton adduct when charge_adducts is None, so that
        # internal charge is never mistaken for extra external protons. Defaults to charge_state
        # (i.e. "assume it's all external protonation") when not given explicitly, matching direct
        # construction of a Fragment outside the internal internal+external charge-splitting pipeline.
        self.external_charge: int = external_charge if external_charge is not None else charge_state
        # int means 13C count
        self._isotopes: Mapping[str, int] | int | None = isotopes
        self._losses: Mapping[str | float, int] | None = deltas
        # Optional composition cache
        self._composition: Counter[ElementInfo] | None = composition
        self.parent_sequence: str | None = parent_sequence
        self.parent_sequence_length: int | None = parent_sequence_length

    @property
    def composition(self) -> Counter[ElementInfo] | None:
        if self._composition is not None:
            return self._composition

        if self.parent_sequence is None:
            raise PeptacularError("Cannot calculate composition without parent sequence or explicit composition")

        if self.parent_sequence_length is None:
            raise PeptacularError("Cannot calculate composition without parent sequence length")

        from .annotation import ProFormaAnnotation

        annot = ProFormaAnnotation.parse(self.parent_sequence)

        pos = validate_position(self.ion_type, self.position, self.parent_sequence_length)
        if pos is not None:
            start, end = pos
            annot = annot[slice(start, end)]

        # Apply this fragment's ion-type offset to the (sliced) sub-sequence. Without it the
        # composition defaulted to the sub-sequence's *precursor* composition, which is heavier
        # than the actual fragment by the ion-type offset (e.g. +H2O for a b-ion), disagreeing
        # with the fragment's own `.mass`. y-ions coincidentally matched (y neutral == precursor).
        charge = self.external_charge if self._charge_adducts is None else self.charge_adducts
        return annot.comp(ion_type=self.ion_type, isotopes=self.isotopes, deltas=self.losses, charge=charge)  # type: ignore

    @property
    def mz(self) -> float:
        return self.mass / abs(self.charge_state) if self.charge_state != 0 else self.mass

    @property
    def neutral_mass(self) -> float:
        # subtract adduct masses and add back the electrons removed by the charge:
        # self.mass == neutral + adduct_atoms - charge*electron, so the electron term
        # must be undone to recover the true neutral mass.
        total_adduct_mass = 0.0
        for adduct in self.charge_adducts:
            total_adduct_mass += adduct.get_mass(self.monoisotopic)
        return self.mass - total_adduct_mass + self.charge_state * ELECTRON_MASS

    @property
    def charge_adducts(self) -> Mods[GlobalChargeCarrier]:
        if self._charge_adducts is None:
            if self.external_charge != 0:
                return Mods[GlobalChargeCarrier](
                    mod_type=ModType.CHARGE,
                    _mods={GlobalChargeCarrier.charged_proton(self.external_charge).serialize(): 1},
                )
            # no real external adducts (charge is entirely internal, or there is no charge at all)
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={})

        # we have adducts, convert to Mods object
        else:
            # Tally identical adduct strings into counts so repeated carriers (e.g. two
            # 'Na:z+1') are not collapsed to one; Mod scales mass/composition by count.
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods=dict(Counter(self._charge_adducts)))

    @property
    def is_protonated(self) -> bool:
        if self._charge_adducts is None and self.external_charge != 0:
            return True
        return False

    @property
    def isotopes(self) -> Mapping[ElementInfo, int]:
        if self._isotopes is not None:
            if isinstance(self._isotopes, int):
                return {ELEMENT_LOOKUP["13C"]: self._isotopes}

            if isinstance(self._isotopes, dict):
                return {ELEMENT_LOOKUP[elem]: count for elem, count in self._isotopes.items()}

        return {}

    @property
    def is_c13(self) -> bool:
        if self._isotopes is not None and isinstance(self._isotopes, int):
            return True
        return False

    @property
    def losses(self) -> Mapping[ChargedFormula | float, int]:
        if self._losses is not None:
            losses = {}
            for loss_name, count in self._losses.items():
                if isinstance(loss_name, float | int):
                    losses[loss_name] = count
                    continue
                loss_formula = ChargedFormula.from_string(loss_name, require_formula_prefix=False)
                losses[loss_formula] = count
            return losses
        return {}

    def asdict(self) -> dict[str, Any]:
        return {
            "ion_type": self.ion_type,
            "position": self.position,
            "mass": self.mass,
            "charge_state": self.charge_state,
            "monoisotopic": self.monoisotopic,
            "charge_adducts": self.charge_adducts,
            "isotopes": self.isotopes,
            "losses": self.losses,
        }

    def serialize(self, format: Literal["default", "mzpaf"] = "default", include_sequence: bool = True) -> str:
        """Serialize the fragment to a string representation.

        :param format: Output format. ``"default"`` returns the human-readable representation,
            ``"mzpaf"`` returns the mzPAF (Peak Annotation Format) label string.
        :type format: Literal["default", "mzpaf"]
        :param include_sequence: If True, include the peptide sequence in the mzPAF label. Only used for ``"mzpaf"`` format.
        :type include_sequence: bool
        :return: The serialized fragment string.
        :rtype: str
        """
        if format == "default":
            return str(self)
        elif format == "mzpaf":
            return self._serialize_mzpaf(include_sequence=include_sequence)
        else:
            raise PeptacularError(f"Unknown format: {format!r}. Use 'default' or 'mzpaf'.")

    def to_mzpaf(self, include_sequence: bool = True) -> str:
        """Serialize the fragment to an mzPAF (Peak Annotation Format) label string.

        mzPAF 1.0.1 ``z`` is the z-dot radical (``IonType.Z_RADICAL``). The Biemann ``z``
        (``IonType.Z``), ``z+H`` and ``c-H`` ions have no letter of their own and are
        written as ``z``/``c`` with a hydrogen delta: ``z3{IDE}-H``, ``z3{IDE}+H``,
        ``c3{PEP}-H``, so the label parses back to the same m/z.

        :param include_sequence: If True, include the peptide sequence in the label.
        :type include_sequence: bool
        :return: The mzPAF label string (e.g. ``"y3{IDE}^2"``).
        :rtype: str
        """
        return self._serialize_mzpaf(include_sequence=include_sequence)

    def _serialize_mzpaf(self, include_sequence: bool = True) -> str:
        """Build the mzPAF label string for this fragment."""
        from .annotation import ProFormaAnnotation

        parts: list[str] = []
        internal_loss: str | None = None
        series_delta: str | None = None

        if self.ion_type is None:
            parts.append("?")
        else:
            ion_info = FRAGMENT_ION_LOOKUP[self.ion_type]

            if ion_info.properties & (IonTypeProperty.FORWARD | IonTypeProperty.BACKWARD):
                # Peptide series ions (a, b, c, x, y, z, d, w, da, db, wa, wb)
                series_str = _ION_TYPE_TO_MZPAF_SERIES.get(ion_info.ion_type, ion_info.ion_type.value)
                series_delta = _MZPAF_SERIES_DELTA.get(ion_info.ion_type)
                position = self.position if isinstance(self.position, int) else -1
                parts.append(f"{series_str}{position}")

                if include_sequence and self.parent_sequence is not None:
                    seq = self.sequence
                    if seq is not None:
                        seq_no_charge = ProFormaAnnotation.parse(seq).serialize(exclude_charge=True)
                        parts.append(f"{{{seq_no_charge}}}")

            elif ion_info.properties & IonTypeProperty.INTERNAL:
                if ion_info.id == IonType.IMMONIUM:
                    # Immonium ion: I{amino_acid}[{modification}]
                    if self.parent_sequence is not None:
                        seq = self.sequence
                        if seq is not None:
                            annot = ProFormaAnnotation.parse(seq)
                            parts.append(f"I{annot.sequence}")

                            if annot.has_internal_mods_at_index(0):
                                internal_mods = annot.get_internal_mods_at_index(0)
                                if len(internal_mods) > 1:
                                    raise PeptacularError(f"Multiple internal mods on immonium ion not supported in mzPAF, got {internal_mods}")
                                if len(internal_mods) == 1 and internal_mods.mods[0].count > 1:
                                    raise PeptacularError(f"Multiple occurrences of internal mod on immonium ion not supported in mzPAF, got {internal_mods}")
                                mods_str = internal_mods.serialize()[1:-1]  # remove surrounding brackets
                                if mods_str == "":
                                    raise PeptacularError(f"Empty modification string for immonium ion is not valid in mzPAF. Internal mods: {internal_mods}")
                                parts.append(f"[{mods_str}]")
                        else:
                            raise PeptacularError("Immonium ion must have a sequence annotation.")
                    else:
                        raise PeptacularError("Immonium ion must have a parent sequence.")
                else:
                    # Internal fragment: m{start}:{end}[{sequence}]
                    if isinstance(self.position, tuple) and len(self.position) == 2:
                        start, end = self.position
                    else:
                        start, end = -1, -1

                    parts.append(f"m{start}:{end}")

                    if include_sequence and self.parent_sequence is not None:
                        seq = self.sequence
                        if seq is not None:
                            seq_no_charge = ProFormaAnnotation.parse(seq).serialize(exclude_charge=True)
                            parts.append(f"{{{seq_no_charge}}}")

                    # Add internal mass diff neutral loss
                    ion_value = ion_info.ion_type.value
                    internal_ion_key = (ion_value[0], ion_value[1]) if len(ion_value) == 2 else None
                    if internal_ion_key is not None and internal_ion_key in _INTERNAL_MASS_DIFFS:
                        internal_loss = _INTERNAL_MASS_DIFFS[internal_ion_key]
                    else:
                        raise PeptacularError(f"Internal ion type {ion_info.ion_type} not supported in mzPAF.")

            elif ion_info.properties & IonTypeProperty.INTACT:
                if ion_info.ion_type == IonType.PRECURSOR:
                    parts.append("p")
                else:
                    raise PeptacularError(f"Cannot convert intact ion type {ion_info.id} to mzPAF.")
            else:
                raise PeptacularError(f"Cannot convert fragment with ion type {self.ion_type} to mzPAF.")

        # Hydrogen delta of a z/c variant that mzPAF writes as its parent series
        if series_delta is not None:
            parts.append(series_delta)

        # Neutral losses from self.losses
        if self._losses is not None:
            for loss_key, count in self._losses.items():
                if isinstance(loss_key, float | int):
                    # mzPAF's neutral_loss grammar only accepts a chemical formula or a
                    # bracketed reference-group name after the sign (spec section 4.5);
                    # there is no representation for an arbitrary unnamed mass delta.
                    raise PeptacularError(
                        f"Cannot convert numeric neutral loss/gain delta ({loss_key!r}) to mzPAF: "
                        "mzPAF neutral losses must be a chemical formula or a named reference group, "
                        "not a bare mass delta."
                    )
                else:
                    loss_formula = ChargedFormula.from_string(loss_key, require_formula_prefix=False)
                    paf_formula = loss_formula.to_mz_paf()
                    sign = paf_formula[0]
                    if sign not in ("+", "-"):
                        raise PeptacularError(f"Invalid formula sign in loss: {paf_formula}")
                    mult = 1 if sign == "+" else -1
                    abs_count = abs(count * mult)
                    count_str = str(abs_count) if abs_count > 1 else ""
                    parts.append(f"{sign}{count_str}{paf_formula[1:]}")

        # Internal mass diff loss (from internal fragment type)
        if internal_loss is not None:
            parts.append(internal_loss)

        # Isotopes
        if self._isotopes is not None:
            if isinstance(self._isotopes, int):
                count_str = str(self._isotopes) if self._isotopes > 1 else ""
                parts.append(f"+{count_str}i")
            elif isinstance(self._isotopes, dict):
                for elem, count in self._isotopes.items():
                    count_str = str(count) if count > 1 else ""
                    parts.append(f"+{count_str}i{elem}")

        # Adducts
        if self._charge_adducts is not None:
            adduct_parts: list[str] = []
            for mod in self.charge_adducts.mods:
                carrier: GlobalChargeCarrier = mod.value
                # A repeated carrier is tallied into mod.count (e.g. two 'Na:z+1' list
                # entries -> count=2), so fold that into the carrier's own occurance;
                # otherwise the mzPAF repeat-count prefix would show only one copy while
                # the mass/charge (which scale by mod.count) show all of them.
                if mod.count != 1:
                    carrier = GlobalChargeCarrier(charged_formula=carrier.charged_formula, occurance=carrier.occurance * mod.count)
                # to_mz_paf() returns "M+Na", we strip the "M" prefix
                paf_str = carrier.to_mz_paf()
                adduct_parts.append(paf_str[1:])  # strip "M", keep "+Na"

            def _adduct_sort_key(part: str) -> str:
                # mzPAF section 4.7: "If there are multiple types of atoms/molecules,
                # alphabetical order SHOULD be followed, e.g. [M+2H+Na] rather than
                # [M+Na+2H]." Sort on the element/molecule name, ignoring the leading
                # sign and any repeat-count digits (e.g. "+2H" sorts as "H").
                name = part[1:].lstrip("0123456789")
                return name

            adduct_parts.sort(key=_adduct_sort_key)
            parts.append(f"[M{''.join(adduct_parts)}]")

        # Charge: mzPAF omits the component only for +1 (implicit); everything else,
        # including negative charges, is written as a bare magnitude with no sign
        # (mzPAF spec section 4.8: "The charge state component ... MUST NOT include
        # the minus sign").
        if self.charge_state != 0 and self.charge_state != 1:
            parts.append(f"^{abs(self.charge_state)}")

        return "".join(parts)

    def __str__(self) -> str:
        parts = []
        parts.append(f"ion_type={self.ion_type}")

        if self.position is not None:
            parts.append(f"position={self.position}")

        parts.append(f"mass={self.mass:.4f}")
        parts.append(f"charge={self.charge_state}")

        # Only show charge_adducts if not simple protonation
        if not self.is_protonated and self.charge_state != 0:
            parts.append(f"charge_adducts={self.charge_adducts}")

        if self.isotopes:
            parts.append(f"isotopes={dict(self.isotopes)}")

        if self.losses:
            parts.append(f"losses={dict(self.losses)}")

        return f"Fragment({', '.join(parts)})"

    def __repr__(self) -> str:
        return (
            f"Fragment(ion_type={self.ion_type!r}, position={self.position!r}, "
            f"mass={self.mass}, monoisotopic={self.monoisotopic}, "
            f"charge_state={self.charge_state}, charge_adducts={self._charge_adducts!r}, "
            f"isotopes={self._isotopes!r}, deltas={self._losses!r}, "
            f"composition={self._composition!r}, parent_sequence={self.parent_sequence!r}, "
            f"parent_sequence_length={self.parent_sequence_length})"
        )

    @property
    def sequence(self) -> str | None:
        if self.parent_sequence is None:
            raise PeptacularError("Cannot determine fragment sequence without parent sequence")

        if self.parent_sequence_length is None:
            raise PeptacularError("Cannot determine fragment sequence without parent sequence length")

        pos = validate_position(self.ion_type, self.position, self.parent_sequence_length)
        if pos is None:
            return self.parent_sequence
        elif isinstance(pos, tuple):
            from .annotation import ProFormaAnnotation

            start, end = pos
            return (
                ProFormaAnnotation.parse(self.parent_sequence)[slice(start, end)]
                .set_charge(self.external_charge if self._charge_adducts is None else self.charge_adducts)
                .serialize()
            )

        raise PeptacularError("Invalid position format for fragment sequence extraction")

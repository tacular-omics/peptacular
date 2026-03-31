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

from ..constants import ModType
from ..proforma_components import (
    ChargedFormula,
    GlobalChargeCarrier,
)
from .mod import Mods
from .positions import validate_position

# Maps internal ion type value tuples to their neutral loss diff relative to "by" (the default internal fragment).
# None means no difference from "by".
_INTERNAL_MASS_DIFFS: dict[tuple[str, str], str | None] = {
    ("a", "x"): None,
    ("b", "x"): "+CO",
    ("c", "x"): "+CHNO",
    ("a", "y"): "-CO",
    ("b", "y"): None,
    ("c", "y"): "+NH",
    ("a", "z"): "-CHNO",
    ("b", "z"): "-NH",
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
}


class Fragment:
    def __init__(
        self,
        ion_type: IonType,
        position: int | tuple[int, int] | None,
        mass: float,
        monoisotopic: bool,
        charge_state: int,
        charge_adducts: tuple[str, ...] | None = None,
        isotopes: Mapping[str, int] | int | None = None,
        deltas: Mapping[str | float, int] | None = None,
        composition: Mapping[ElementInfo, int] | None = None,
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
            raise ValueError("Cannot calculate composition without parent sequence or explicit composition")

        if self.parent_sequence_length is None:
            raise ValueError("Cannot calculate composition without parent sequence length")

        from .annotation import ProFormaAnnotation

        annot = ProFormaAnnotation.parse(self.parent_sequence)

        pos = validate_position(self.ion_type, self.position, self.parent_sequence_length)
        if pos is not None:
            start, end = pos
            annot = annot[slice(start, end)]

        return annot.comp(isotopes=self.isotopes, deltas=self.losses, charge=self.charge_state if self._charge_adducts is None else self.charge_adducts)  # type: ignore

    @property
    def mz(self) -> float:
        return self.mass / abs(self.charge_state) if self.charge_state != 0 else self.mass

    @property
    def neutral_mass(self) -> float:
        # subract adduct masses
        total_adduct_mass = 0.0
        for adduct in self.charge_adducts:
            total_adduct_mass += adduct.get_mass(self.monoisotopic)
        return self.mass - total_adduct_mass

    @property
    def charge_adducts(self) -> Mods[GlobalChargeCarrier]:
        if self._charge_adducts is None:
            if self.charge_state != 0:
                return Mods[GlobalChargeCarrier](
                    mod_type=ModType.CHARGE,
                    _mods={GlobalChargeCarrier.charged_proton(self.charge_state).serialize(): 1},
                )
            # no adducts no charge
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={})

        # we have adducts, convert to Mods object
        else:
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={k: 1 for k in self._charge_adducts})

    @property
    def is_protonated(self) -> bool:
        if self._charge_adducts is None and self.charge_state != 0:
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
            raise ValueError(f"Unknown format: {format!r}. Use 'default' or 'mzpaf'.")

    def to_mzpaf(self, include_sequence: bool = True) -> str:
        """Serialize the fragment to an mzPAF (Peak Annotation Format) label string.

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

        if self.ion_type is None:
            parts.append("?")
        else:
            ion_info = FRAGMENT_ION_LOOKUP[self.ion_type]

            if ion_info.properties & (IonTypeProperty.FORWARD | IonTypeProperty.BACKWARD):
                # Peptide series ions (a, b, c, x, y, z, d, w, da, db, wa, wb)
                series_str = _ION_TYPE_TO_MZPAF_SERIES.get(ion_info.ion_type, ion_info.ion_type.value)
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
                                    raise ValueError(f"Multiple internal mods on immonium ion not supported in mzPAF, got {internal_mods}")
                                if len(internal_mods) == 1 and internal_mods.mods[0].count > 1:
                                    raise ValueError(f"Multiple occurrences of internal mod on immonium ion not supported in mzPAF, got {internal_mods}")
                                mods_str = internal_mods.serialize()[1:-1]  # remove surrounding brackets
                                if mods_str == "":
                                    raise ValueError(f"Empty modification string for immonium ion is not valid in mzPAF. Internal mods: {internal_mods}")
                                parts.append(f"[{mods_str}]")
                        else:
                            raise ValueError("Immonium ion must have a sequence annotation.")
                    else:
                        raise ValueError("Immonium ion must have a parent sequence.")
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
                    internal_ion_key = tuple(list(ion_info.ion_type.value))
                    if internal_ion_key in _INTERNAL_MASS_DIFFS:
                        internal_loss = _INTERNAL_MASS_DIFFS[internal_ion_key]
                    else:
                        raise ValueError(f"Internal ion type {ion_info.ion_type} not supported in mzPAF.")

            elif ion_info.properties & IonTypeProperty.INTACT:
                if ion_info.ion_type == IonType.PRECURSOR:
                    parts.append("p")
                else:
                    raise ValueError(f"Cannot convert intact ion type {ion_info.id} to mzPAF.")
            else:
                raise ValueError(f"Cannot convert fragment with ion type {self.ion_type} to mzPAF.")

        # Neutral losses from self.losses
        if self._losses is not None:
            for loss_key, count in self._losses.items():
                if isinstance(loss_key, float | int):
                    mass_val = float(loss_key) * count
                    parts.append(f"{mass_val:+.5f}")
                else:
                    loss_formula = ChargedFormula.from_string(loss_key, require_formula_prefix=False)
                    paf_formula = loss_formula.to_mz_paf()
                    sign = paf_formula[0]
                    if sign not in ("+", "-"):
                        raise ValueError(f"Invalid formula sign in loss: {paf_formula}")
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
                # to_mz_paf() returns "M+Na", we strip the "M" prefix
                paf_str = carrier.to_mz_paf()
                adduct_parts.append(paf_str[1:])  # strip "M", keep "+Na"
            parts.append(f"[M{''.join(adduct_parts)}]")

        # Charge: omit for +1 (implicit); include for everything else
        if self.charge_state != 0 and self.charge_state != 1:
            parts.append(f"^{self.charge_state}")

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
            raise ValueError("Cannot determine fragment sequence without parent sequence")

        if self.parent_sequence_length is None:
            raise ValueError("Cannot determine fragment sequence without parent sequence length")

        pos = validate_position(self.ion_type, self.position, self.parent_sequence_length)
        if pos is None:
            return self.parent_sequence
        elif isinstance(pos, tuple):
            from .annotation import ProFormaAnnotation

            start, end = pos
            return (
                ProFormaAnnotation.parse(self.parent_sequence)[slice(start, end)]
                .set_charge(self.charge_state if self._charge_adducts is None else self.charge_adducts)
                .serialize()
            )

        raise ValueError("Invalid position format for fragment sequence extraction")

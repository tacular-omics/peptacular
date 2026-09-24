from collections import Counter
from collections.abc import Mapping
from dataclasses import FrozenInstanceError
from functools import cache
from typing import Any, Literal

from tacular import (
    ELEMENT_LOOKUP,
    FRAGMENT_ION_LOOKUP,
    NEUTRAL_DELTA_LOOKUP,
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


_FormulaKey = frozenset[tuple[str, int | None, int]]


def _formula_key(formula: ChargedFormula) -> _FormulaKey:
    """The formula's element magnitudes, ignoring sign and written order (``HCONH2`` == ``CH3NO``)."""
    counts: Counter[tuple[str, int | None]] = Counter()
    for fe in formula.formula:
        counts[(fe.element.value, fe.isotope)] += fe.occurance
    return frozenset((element, isotope, abs(count)) for (element, isotope), count in counts.items())


@cache
def _named_neutral_deltas() -> dict[_FormulaKey, str]:
    """Composition -> the canonical mzPAF name of each known neutral delta (``NH3``, not ``H3N``)."""
    named: dict[_FormulaKey, str] = {}
    for info in NEUTRAL_DELTA_LOOKUP.values():
        formula = ChargedFormula.from_string(info.formula, require_formula_prefix=False)
        named.setdefault(_formula_key(formula), info.formula)
    return named


def _hill_rank(formula: ChargedFormula) -> Any:
    has_carbon = any(fe.element.value == "C" for fe in formula.formula)

    def rank(fe: Any) -> tuple[int, str, int]:
        symbol = fe.element.value
        if has_carbon and symbol in ("C", "H"):
            return (0 if symbol == "C" else 1, "", fe.isotope or 0)
        return (2, symbol, fe.isotope or 0)

    return rank


def _mzpaf_formula(formula: ChargedFormula) -> str:
    """A loss/gain formula as a signed mzPAF token: the canonical name of a known neutral delta
    (``-NH3``, ``-H2O``, ``-H3PO4``), else the formula in Hill order. mzPAF section 4.5 forbids
    ``H3N`` for ammonia, and tacular stores compositions H-first, so the written order cannot be used.
    """
    signs = {fe.occurance > 0 for fe in formula.formula}
    if len(signs) != 1:
        return formula.to_mz_paf()  # raises the mixed-sign / empty error
    sign = "+" if signs.pop() else "-"
    name = _named_neutral_deltas().get(_formula_key(formula))
    if name is not None:
        return sign + name
    ordered = sorted(formula.formula, key=_hill_rank(formula))
    return sign + "".join(str(fe.abs()) for fe in ordered)


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


_FRAGMENT_FIELDS: dict[str, str] = {
    "ion_type": "ion_type",
    "position": "position",
    "mass": "mass",
    "monoisotopic": "monoisotopic",
    "charge_state": "charge_state",
    "charge_adducts": "_charge_adducts",
    "external_charge": "external_charge",
    "isotopes": "_isotopes",
    "deltas": "_deltas",
    "composition": "_composition",
    "parent_sequence": "parent_sequence",
    "parent_sequence_length": "parent_sequence_length",
}
_COMPOSITION_INPUTS = frozenset(_FRAGMENT_FIELDS) - {"mass", "monoisotopic", "composition"}


def _freeze(value: Mapping[Any, int] | int | None) -> Any:
    """Hashable form of an isotope or delta mapping."""
    if isinstance(value, Mapping):
        return frozenset(value.items())
    return value


def _mzpaf_mass(value: float) -> str:
    """A signed mass for mzPAF: fixed-point, 6 decimals, trailing zeros stripped (``-34.0``, ``+1e-05`` -> ``+0.00001``)."""
    text = f"{value:+.6f}".rstrip("0")
    if text.endswith("."):
        text += "0"
    return "+0.0" if text == "-0.0" else text  # the caller drops a zero delta


def _delta_keys(deltas: Mapping[Any, int] | None) -> Mapping[str | float, int] | None:
    """Store delta keys as strings or masses, so ``Fragment(deltas=frag.deltas)`` round-trips."""
    if not deltas:
        return None  # ``frag.deltas`` reports "no deltas" as {}
    if all(isinstance(key, str | float | int) for key in deltas):
        return deltas
    out: dict[str | float, int] = {}
    for key, count in deltas.items():
        if isinstance(key, ChargedFormula):
            key = key.serialize().removeprefix("Formula:")
        out[key] = out.get(key, 0) + count
    return out


def _isotope_keys(isotopes: Mapping[Any, int] | int | None) -> Mapping[str, int] | int | None:
    """Store isotope keys as strings (``"15N"``), so ``Fragment(isotopes=frag.isotopes)`` round-trips."""
    if isinstance(isotopes, int) or isotopes is None:
        return isotopes
    if not isotopes:
        return None  # ``frag.isotopes`` reports "no isotopes" as {}
    return {str(key): count for key, count in isotopes.items()}


def _adduct_strings(adducts: Any) -> tuple[str, ...] | None:
    """Store charge adducts as a tuple of strings; a :class:`Mods` (``frag.charge_adducts``) is expanded."""
    if isinstance(adducts, Mods):
        return tuple(key for key, count in (adducts._mods or {}).items() for _ in range(count))
    return adducts


class Fragment:
    """One theoretical ion: a fragment or precursor with its ion type, position, charge and mass.

    Returned by :meth:`ProFormaAnnotation.frag`, :meth:`ProFormaAnnotation.fragment` and
    :func:`peptacular.fragment`. ``mass`` is the mass of the charged ion (adducts included), so
    ``mz`` is ``mass / abs(charge_state)`` and ``neutral_mass`` removes the charge carriers.
    Use :meth:`to_mzpaf` for an mzPAF annotation string. Fragments are immutable: assigning
    to an attribute raises :class:`dataclasses.FrozenInstanceError`; :meth:`replace` returns a
    changed copy. Fragments compare and hash by value.

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

    __slots__ = (
        "ion_type",
        "position",
        "mass",
        "monoisotopic",
        "charge_state",
        "_charge_adducts",
        "external_charge",
        "_isotopes",
        "_deltas",
        "_composition",
        "parent_sequence",
        "parent_sequence_length",
    )

    ion_type: IonType
    position: int | tuple[int, int] | None
    mass: int | float
    monoisotopic: bool
    charge_state: int
    _charge_adducts: tuple[str, ...] | None
    external_charge: int
    _isotopes: Mapping[str, int] | int | None
    _deltas: Mapping[str | float, int] | None
    _composition: Counter[ElementInfo] | None
    parent_sequence: str | None
    parent_sequence_length: int | None

    def __init__(
        self,
        ion_type: IonType,
        position: int | tuple[int, int] | None,
        mass: float,
        monoisotopic: bool,
        charge_state: int,
        *,
        charge_adducts: tuple[str, ...] | None = None,
        external_charge: int | None = None,
        isotopes: Mapping[str, int] | int | None = None,
        deltas: Mapping[str | float, int] | None = None,
        composition: Counter[ElementInfo] | None = None,
        parent_sequence: str | None = None,
        parent_sequence_length: int | None = None,
    ) -> None:
        _set = object.__setattr__
        _set(self, "ion_type", ion_type)
        _set(self, "position", position)
        _set(self, "mass", mass)
        _set(self, "monoisotopic", monoisotopic)
        _set(self, "charge_state", charge_state)
        # If None and charge_state != 0: means protonated
        _set(self, "_charge_adducts", _adduct_strings(charge_adducts))
        # The portion of charge_state that comes from real external adducts/charge carriers,
        # as opposed to charge intrinsic to an internal formula modification (e.g. [Formula:...:z+N]).
        # Used to reconstruct the default proton adduct when charge_adducts is None, so that
        # internal charge is never mistaken for extra external protons. Defaults to charge_state
        # (i.e. "assume it's all external protonation") when not given explicitly.
        _set(self, "external_charge", external_charge if external_charge is not None else charge_state)
        # int means 13C count
        _set(self, "_isotopes", _isotope_keys(isotopes))
        _set(self, "_deltas", _delta_keys(deltas))
        # Optional composition cache
        _set(self, "_composition", composition)
        _set(self, "parent_sequence", parent_sequence)
        _set(self, "parent_sequence_length", parent_sequence_length)

    def __setattr__(self, name: str, value: object) -> None:
        raise FrozenInstanceError(f"cannot assign to field {name!r}: Fragment is immutable")

    def __delattr__(self, name: str) -> None:
        raise FrozenInstanceError(f"cannot delete field {name!r}: Fragment is immutable")

    def __getstate__(self) -> dict[str, Any]:
        return {name: getattr(self, name) for name in self.__slots__}

    def __setstate__(self, state: dict[str, Any]) -> None:
        for name, value in state.items():
            object.__setattr__(self, name, value)

    def _replace(self, **changes: Any) -> "Fragment":
        """Return a copy with the given slot values replaced (internal helper)."""
        state = self.__getstate__()
        state.update(changes)
        new = Fragment.__new__(Fragment)
        new.__setstate__(state)
        return new

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
        return annot.comp(ion_type=self.ion_type, isotopes=self.isotopes, deltas=self.deltas, charge=charge)  # type: ignore

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
            total_adduct_mass += adduct.get_mass(monoisotopic=self.monoisotopic)
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
    def deltas(self) -> Mapping[ChargedFormula | float, int]:
        """Neutral losses and gains as ``{ChargedFormula or mass: count}``.

        A named loss such as ``H2O`` is stored as a negative formula (``H-2O-1``); a plain
        formula or a positive mass is a gain.
        """
        if self._deltas is not None:
            deltas: dict[ChargedFormula | float, int] = {}
            for key, count in self._deltas.items():
                if isinstance(key, float | int):
                    deltas[key] = count
                    continue
                deltas[ChargedFormula.from_string(key, require_formula_prefix=False)] = count
            return deltas
        return {}

    def _value_key(self) -> tuple[Any, ...]:
        """The fields that define a fragment's value (the composition cache is excluded)."""
        return (
            self.ion_type,
            self.position,
            self.mass,
            self.monoisotopic,
            self.charge_state,
            self._charge_adducts,
            self.external_charge,
            _freeze(self._isotopes),
            _freeze(self._deltas),
            self.parent_sequence,
            self.parent_sequence_length,
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Fragment):
            return NotImplemented
        return self._value_key() == other._value_key()

    def __hash__(self) -> int:
        return hash(self._value_key())

    def replace(self, **changes: Any) -> "Fragment":
        """Return a copy with the given constructor arguments replaced.

        Takes the :class:`Fragment` constructor names (``mass``, ``charge_state``, ``deltas``,
        ...). The cached composition is dropped when a field it depends on changes, unless
        ``composition`` is passed too. Changing ``charge_state`` of a fragment whose charge is
        all external also moves ``external_charge``. The values the properties return
        (``frag.deltas``, ``frag.isotopes``, ``frag.charge_adducts``) are accepted as is.

        >>> import peptacular as pt
        >>> frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=2)
        >>> frag.replace(mass=100.0).mass
        100.0
        >>> frag.replace(mass=frag.mass) == frag
        True

        :raises TypeError: If a key is not a constructor argument.
        :return: A new fragment.
        :rtype: Fragment
        """
        unknown = set(changes) - _FRAGMENT_FIELDS.keys()
        if unknown:
            raise TypeError(f"Fragment.replace() got unexpected field(s) {sorted(unknown)}; expected any of {sorted(_FRAGMENT_FIELDS)}")
        if "isotopes" in changes:
            isotopes = changes["isotopes"]
            # `isotopes` reports a 13C count as {13C: n}; passing that back keeps the count
            changes["isotopes"] = self._isotopes if isotopes == self.isotopes else _isotope_keys(isotopes)
        if "deltas" in changes:
            changes["deltas"] = _delta_keys(changes["deltas"])
        if "charge_adducts" in changes:
            adducts = changes["charge_adducts"]
            # the default protons that `charge_adducts` reports for a protonated ion stay implicit
            keep_default = self._charge_adducts is None and isinstance(adducts, Mods) and adducts == self.charge_adducts
            changes["charge_adducts"] = None if keep_default else _adduct_strings(adducts)
        slot_changes = {_FRAGMENT_FIELDS[name]: value for name, value in changes.items()}
        if "external_charge" in changes and changes["external_charge"] is None:
            slot_changes["external_charge"] = changes.get("charge_state", self.charge_state)
        elif "charge_state" in changes and "external_charge" not in changes and self.external_charge == self.charge_state:
            slot_changes["external_charge"] = changes["charge_state"]
        if "composition" not in changes and _COMPOSITION_INPUTS & changes.keys():
            slot_changes["_composition"] = None
        return self._replace(**slot_changes)

    def asdict(self) -> dict[str, Any]:
        return {
            "ion_type": self.ion_type,
            "position": self.position,
            "mass": self.mass,
            "charge_state": self.charge_state,
            "monoisotopic": self.monoisotopic,
            "charge_adducts": self.charge_adducts,
            "isotopes": self.isotopes,
            "deltas": self.deltas,
        }

    def serialize(self, *, format: Literal["default", "mzpaf"] = "default", include_sequence: bool = True, signed_charge: bool = True) -> str:
        """Serialize the fragment to a string representation.

        :param format: Output format. ``"default"`` returns the human-readable representation,
            ``"mzpaf"`` returns the mzPAF (Peak Annotation Format) label string.
        :type format: Literal["default", "mzpaf"]
        :param include_sequence: If True, include the peptide sequence in the mzPAF label. Only used for ``"mzpaf"`` format.
        :type include_sequence: bool
        :param signed_charge: If True, write a negative charge as ``^-n`` in the mzPAF label. Only used for ``"mzpaf"`` format.
        :type signed_charge: bool
        :return: The serialized fragment string.
        :rtype: str
        """
        if format == "default":
            return str(self)
        elif format == "mzpaf":
            return self._serialize_mzpaf(include_sequence=include_sequence, signed_charge=signed_charge)
        else:
            raise PeptacularError(f"Unknown format: {format!r}. Use 'default' or 'mzpaf'.")

    def to_mzpaf(self, *, include_sequence: bool = True, signed_charge: bool = True) -> str:
        """Serialize the fragment to an mzPAF (Peak Annotation Format) label string.

        mzPAF 1.0.1 ``z`` is the z-dot radical (``IonType.Z_RADICAL``). The Biemann ``z``
        (``IonType.Z``), ``z+H`` and ``c-H`` ions have no letter of their own and are
        written as ``z``/``c`` with a hydrogen delta: ``z3{IDE}-H``, ``z3{IDE}+H``,
        ``c3{PEP}-H``, so the label parses back to the same m/z.

        Each delta is a signed formula or mass, added ``count`` times. Named losses such as
        H2O are stored as negative formulas (``H-2O-1``), so ``{"H2O": 1}`` is written
        ``-H2O`` and ``{"H2O": -1}`` ``+H2O``; a plain formula such as ``C2H2O`` is a gain
        (``+C2H2O``). A numeric delta is written as a signed fixed-point mass rounded to 6
        decimals (``b2-34.0`` for ``{-17.0: 2}``); one that rounds to zero is left out. A formula with both positive and negative
        element counts (``CH-2``) cannot be written and raises :class:`PeptacularError`.

        A negative charge is written signed (``y3{IDE}^-1``) so the label parses back to the
        same m/z, as paftacular does. mzPAF 1.0.1 section 4.8 says the charge MUST NOT include
        the minus sign (negative mode is a property of the spectrum); pass
        ``signed_charge=False`` to write only the magnitude.

        An immonium ion takes at most one modification (``IP[Oxidation]``). A terminal
        modification on the residue is written there as well, so ``[Acetyl]-PEP`` at position
        1 gives ``IP[Acetyl]``; more than one modification raises :class:`PeptacularError`.

        :param include_sequence: If True, include the peptide sequence in the label.
        :type include_sequence: bool
        :param signed_charge: If True (default), write a negative charge as ``^-n``;
            if False, write ``^n``.
        :type signed_charge: bool
        :return: The mzPAF label string (e.g. ``"y3{IDE}^2"``).
        :rtype: str
        """
        return self._serialize_mzpaf(include_sequence=include_sequence, signed_charge=signed_charge)

    def _serialize_mzpaf(self, *, include_sequence: bool = True, signed_charge: bool = True) -> str:
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

                            # mzPAF allows one modification on an immonium ion. A terminal
                            # modification of the residue (e.g. an N-terminal acetyl) adds the
                            # same mass, so it is written there too (matches paftacular 2.0).
                            tags: list[str] = []
                            for has_mods, get_mods in (
                                (annot.has_internal_mods_at_index(0), lambda: annot.get_internal_mods_at_index(0)),
                                (annot.has_nterm_mods, lambda: annot.nterm_mods),
                                (annot.has_cterm_mods, lambda: annot.cterm_mods),
                            ):
                                if has_mods:
                                    for mod in get_mods().mods:
                                        tags.extend([str(mod.value)] * mod.count)
                            if len(tags) > 1:
                                raise PeptacularError(f"mzPAF allows one modification on an immonium ion, got {', '.join(tags)}")
                            if tags:
                                if tags[0] == "":
                                    raise PeptacularError("Empty modification string for immonium ion is not valid in mzPAF.")
                                parts.append(f"[{tags[0]}]")
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

        # Neutral losses and gains from self._deltas
        if self._deltas is not None:
            for loss_key, count in self._deltas.items():
                if isinstance(loss_key, float | int):
                    # mzPAF (spec section 4.5) writes an unnamed mass delta as a signed number,
                    # ``y8-17.0265``. A number takes no repeat count, so the count is folded in.
                    mass = _mzpaf_mass(loss_key * count)
                    if mass != "+0.0":  # a delta that rounds to zero changes nothing; leave it out
                        parts.append(mass)
                    continue
                loss_formula = ChargedFormula.from_string(loss_key, require_formula_prefix=False)
                paf_formula = _mzpaf_formula(loss_formula)
                sign = paf_formula[0]
                if sign not in ("+", "-"):
                    raise PeptacularError(f"Invalid formula sign in loss: {paf_formula}")
                if count < 0:  # a negative count turns a loss into a gain and back
                    sign = "+" if sign == "-" else "-"
                count_str = str(abs(count)) if abs(count) > 1 else ""
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

        # Charge: mzPAF omits the component only for +1 (implicit). A negative charge is
        # written signed (``^-1``) so the label parses back to the same m/z. mzPAF 1.0.1
        # section 4.8 says the charge MUST NOT include the minus sign (negative mode is a
        # property of the spectrum); ``signed_charge=False`` writes only the magnitude.
        charge = self.charge_state if signed_charge else abs(self.charge_state)
        if charge != 0 and charge != 1:
            parts.append(f"^{charge}")

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

        if self.deltas:
            parts.append(f"deltas={dict(self.deltas)}")

        return f"Fragment({', '.join(parts)})"

    def __repr__(self) -> str:
        return (
            f"Fragment(ion_type={self.ion_type!r}, position={self.position!r}, "
            f"mass={self.mass}, monoisotopic={self.monoisotopic}, "
            f"charge_state={self.charge_state}, charge_adducts={self._charge_adducts!r}, "
            f"isotopes={self._isotopes!r}, deltas={self._deltas!r}, "
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

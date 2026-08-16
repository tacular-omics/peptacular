"""
Core component data structures.

This module contains the basic NamedTuple definitions without any circular dependencies.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Protocol, Self, runtime_checkable

from tacular import (
    AA_LOOKUP,
    ELEMENT_LOOKUP,
    GNO_LOOKUP,
    MONOSACCHARIDE_LOOKUP,
    PSIMOD_LOOKUP,
    RESID_LOOKUP,
    UNIMOD_LOOKUP,
    XLMOD_LOOKUP,
    AminoAcid,
    Element,
    ElementInfo,
    GnoInfo,
    Monosaccharide,
    PsimodInfo,
    ResidInfo,
    UnimodInfo,
    XlModInfo,
)

from ..constants import CV, Terminal

# Reusable hint appended to "unknown modification" errors so callers (including AI
# agents) can immediately see how to specify a resolvable modification.
_MOD_SPEC_HINT = (
    "Specify one of: a known modification name (e.g. 'Oxidation'), a CV accession "
    "(e.g. 'UNIMOD:35' or 'MOD:00046'), a chemical formula (e.g. '[Formula:HO3P]'), "
    "a glycan (e.g. '[Glycan:HexNAc]'), or a delta mass (e.g. '[+15.9949]')."
)


@runtime_checkable
class HasMassComp(Protocol):
    """Protocol for objects that have mass and composition."""

    def get_mass(self, monoisotopic: bool = True) -> float: ...

    def get_composition(self) -> Counter[ElementInfo]: ...


class MassPropertyMixin(ABC):
    """Mixin to add mass properties to classes that implement get_mass()"""

    @abstractmethod
    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get the mass of this component.

        Args:
            monoisotopic: If True, return monoisotopic mass; otherwise average mass.

        Returns:
            Mass in Daltons.
        """
        ...

    @property
    def monoisotopic_mass(self) -> float:
        return self.get_mass(monoisotopic=True)

    @property
    def average_mass(self) -> float:
        return self.get_mass(monoisotopic=False)

    def to_dict(self) -> dict[str, Any]:
        """Return the versioned, JSON-compatible representation of this component."""
        from ..proforma_json import to_proforma_dict

        return to_proforma_dict(self)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Restore this component from its versioned JSON-compatible representation."""
        from ..proforma_json import from_proforma_dict

        return from_proforma_dict(data, expected_type=cls)

    def to_json(self, *, indent: int | None = None) -> str:
        """Return deterministic JSON text for this component."""
        from ..proforma_json import to_proforma_json

        return to_proforma_json(self, indent=indent)

    @classmethod
    def from_json(cls, data: str | bytes | bytearray) -> Self:
        """Restore this component from versioned JSON text."""
        from ..proforma_json import from_proforma_json

        return from_proforma_json(data, expected_type=cls)


def sum_masses(components: Iterable[HasMassComp], monoisotopic: bool = True) -> float:
    """Sum masses from multiple components."""
    return sum(comp.get_mass(monoisotopic=monoisotopic) for comp in components)


def add_composition(total: Counter[ElementInfo], other: Mapping[ElementInfo, int]) -> None:
    """Merge ``other`` into ``total`` in place, preserving negative counts.

    Counter's ``+=`` / ``+`` discard non-positive results, which silently drops
    atom-removing modifications (e.g. ``Formula:H-2`` or Amidated's ``O:-1``). Merging
    element-by-element keeps them. This is the single canonical composition-merge helper;
    every additive merge of element compositions should route through it (or the
    ``merge_compositions`` wrapper below) rather than re-implementing the workaround.
    """
    for element, count in other.items():
        total[element] += count


def merge_compositions(components: Iterable[HasMassComp]) -> Counter[ElementInfo]:
    """Merge compositions from multiple components (negative-count safe)."""
    total: Counter[ElementInfo] = Counter()
    for comp in components:
        add_composition(total, comp.get_composition())
    return total


@runtime_checkable
class HasPositionScore(Protocol):
    """Protocol for objects that have position_id and score attributes."""

    position_id: str | None
    score: float | None

    def serialize_position_score(self) -> str:
        """Serialize the position_id and score components."""
        ...


class PositionScoreMixin(ABC):
    """Mixin to add position/score serialization to classes."""

    position_id: str | None
    score: float | None

    def serialize_position_score(self) -> str:
        """Serialize the position_id and score components.

        Returns empty string if neither position_id nor score are present.
        Returns #<position_id> if only position_id is present.
        Returns #<position_id>(<score>) if both are present.
        Returns (<score>) if only score is present.
        """

        if self.position_id is not None and self.score is not None:
            return f"#{self.position_id}({self.score})"
        elif self.position_id is not None:
            return f"#{self.position_id}"
        elif self.score is not None:
            return f"({self.score})"
        else:
            return ""


@dataclass(frozen=True, slots=True)
class FormulaElement(MassPropertyMixin):
    """A single element in a molecular formula like [13C]2 or H2"""

    element: Element
    occurance: int
    isotope: int | None = None

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        try:
            _ = self.get_mass()
            return None
        except Exception as e:
            return str(e)

    def get_mass(self, monoisotopic: bool = True) -> float:
        if self.isotope is not None:
            monoisotopic = True
        return ELEMENT_LOOKUP[(self.element, self.isotope)].get_mass(monoisotopic=monoisotopic) * self.occurance

    def get_element_count(self) -> tuple[ElementInfo, int]:
        return (ELEMENT_LOOKUP[(self.element, self.isotope)], self.occurance)

    @staticmethod
    def from_element_info(elem_info: ElementInfo, occurance: int) -> FormulaElement:
        s = f"{elem_info}{occurance if occurance != 1 else ''}"
        if elem_info.mass_number is not None:
            s = f"[{s}]"
        return FormulaElement.from_string(s)

    def get_composition(self) -> Counter[ElementInfo]:
        elem_info, occurance = self.get_element_count()
        return Counter({elem_info: occurance})

    @staticmethod
    def from_string(s: str, allow_zero: bool = False) -> FormulaElement:
        from ..proforma_components.parsers import parse_formula_element

        return parse_formula_element(s, allow_zero=allow_zero)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_formula_element

        return serialize_formula_element(self)

    def __str__(self) -> str:
        return self.serialize()

    def abs(self) -> FormulaElement:
        """Get a FormulaElement with absolute occurance."""
        return FormulaElement(element=self.element, occurance=abs(self.occurance), isotope=self.isotope)


@dataclass(frozen=True, slots=True)
class ChargedFormula(MassPropertyMixin, PositionScoreMixin):
    """A formula that can be charged (``<formula>:z<charge>``).

    As a localised residue modification it carries the ``Formula:`` prefix, e.g.
    ``[Formula:C2H6:z+2]``; as a charge carrier it is written bare, e.g. ``C2H6:z+2``
    (see ProForma 2.1 sections 11.1 and 11.5).
    """

    formula: tuple[FormulaElement, ...]
    charge: int | None = None
    position_id: str | None = None
    score: float | None = None

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        try:
            _ = self.get_mass()
            return None
        except Exception as e:
            return str(e)

    def formula_dict(self) -> dict[str, int]:
        from .serializers import get_element_key

        formula_dict: dict[str, int] = {}
        for fe in self.formula:
            key = get_element_key(fe)
            formula_dict[key] = formula_dict.get(key, 0) + fe.occurance
        return formula_dict

    def get_mass(self, monoisotopic: bool = True) -> float:
        mass: float = 0.0
        for elem in self.formula:
            mass += elem.get_mass(monoisotopic=monoisotopic)
        return mass

    def get_composition(self) -> Counter[ElementInfo]:
        composition: Counter[ElementInfo] = Counter()
        for elem in self.formula:
            elem_info, occurance = elem.get_element_count()
            composition[elem_info] += occurance
        return composition

    def get_dict_composition(self) -> dict[str, int]:
        """Get the composition as a dict of element symbols to counts"""
        return {str(elem_info): count for elem_info, count in self.get_composition().items()}

    @staticmethod
    def from_composition(
        composition: Mapping[ElementInfo | str, int] | Counter[ElementInfo],
        charge: int | None = None,
    ) -> ChargedFormula:
        formula_elements: list[FormulaElement] = []
        for elem_info, occurance in composition.items():
            if isinstance(elem_info, str):
                elem_info = ELEMENT_LOOKUP[elem_info]
            formula_elements.append(FormulaElement.from_element_info(elem_info, occurance))
        return ChargedFormula(formula=tuple(formula_elements), charge=charge)

    @staticmethod
    def from_string(
        s: str,
        allow_zero: bool = False,
        require_formula_prefix: bool = True,
        sep: str = "",
    ) -> ChargedFormula:
        from ..proforma_components.parsers import parse_charged_formula

        return parse_charged_formula(
            s,
            allow_zero=allow_zero,
            require_formula_prefix=require_formula_prefix,
            sep=sep,
        )

    def serialize(
        self,
        sep: str = "",
        hill_order: bool = False,
        include_formula_prefix: bool = True,
    ) -> str:
        from ..proforma_components.serializers import serialize_charged_formula

        return serialize_charged_formula(
            self,
            space=sep,
            hill_order=hill_order,
            include_formula_prefix=include_formula_prefix,
        )

    def __str__(self) -> str:
        return self.serialize()

    def to_mz_paf(self) -> str:
        """Convert to mzPAF format string."""
        # mzPAF's chemical-formula notation (secs. 4.4.9/4.5/4.6) explicitly reuses
        # ProForma's own molecular-formula notation -- atom then count (e.g. "H2O",
        # "[13C1]") -- for both plain and isotope-tagged elements. There is no
        # count-before-atom convention; the only leading integer mzPAF defines is a
        # separate repeat-count multiplier for an entire repeated loss group (e.g.
        # "-2H2O" for a double water loss), not a per-atom prefix.
        pos_parts = [str(fe.abs()) for fe in self.formula if fe.occurance > 0]
        neg_parts = [str(fe.abs()) for fe in self.formula if fe.occurance < 0]

        if pos_parts and neg_parts:
            raise ValueError("Cannot convert to mzPAF: contains both positive and negative elements")

        if pos_parts:
            return "+" + "".join(pos_parts)
        if neg_parts:
            return "-" + "".join(neg_parts)

        raise ValueError("Cannot convert to mzPAF: no elements present")

    @staticmethod
    def from_mz_paf(s: str) -> ChargedFormula:
        """Parse from mzPAF format string."""
        # mzPAF chemical formulas are never "Formula:"-prefixed (unlike a ProForma
        # residue modification), so require_formula_prefix must be False here or
        # to_mz_paf()'s own output (e.g. "+H2O") fails to round-trip.
        # split on + and -
        if s.startswith("+"):
            formula = ChargedFormula.from_string(s[1:], require_formula_prefix=False)
            # assert all are positive
            for fe in formula.formula:
                if fe.occurance < 0:
                    raise ValueError("Invalid mzPAF format: negative occurance in positive part")
            return formula
        if s.startswith("-"):
            # Parse the bare (unsigned) formula, then negate every element's count.
            # (The previous "0" + s prepend trick never actually worked: e.g.
            # "0-H2O" isn't a parseable formula either way.)
            formula = ChargedFormula.from_string(s[1:], require_formula_prefix=False)
            # assert none are already negative (would double-negate)
            for fe in formula.formula:
                if fe.occurance < 0:
                    raise ValueError("Invalid mzPAF format: negative occurance in negative part")
            negated = tuple(FormulaElement(element=fe.element, occurance=-fe.occurance, isotope=fe.isotope) for fe in formula.formula)
            return ChargedFormula(
                formula=negated,
                charge=formula.charge,
                position_id=formula.position_id,
                score=formula.score,
            )
        raise ValueError("Invalid mzPAF format: must start with + or -")

    def __add__(self, other: ChargedFormula) -> ChargedFormula:
        """Add two ChargedFormulas together."""
        combined_comp: Counter[ElementInfo] = self.get_composition() + other.get_composition()
        combined_charge = None
        if self.charge is not None and other.charge is not None:
            combined_charge = self.charge + other.charge
        return ChargedFormula.from_composition(combined_comp, charge=combined_charge)

    def __sub__(self, other: ChargedFormula) -> ChargedFormula:
        """Subtract one ChargedFormula from another."""
        combined_comp: Counter[ElementInfo] = self.get_composition() - other.get_composition()
        combined_charge = None
        if self.charge is not None and other.charge is not None:
            combined_charge = self.charge - other.charge
        return ChargedFormula.from_composition(combined_comp, charge=combined_charge)

    @property
    def is_neutral(self) -> bool:
        """Check if the formula is neutral (charge == 0 or None)."""
        return self.charge is None or self.charge == 0

    @property
    def is_charged(self) -> bool:
        """Check if the formula is charged (charge != 0 and not None)."""
        return self.charge is not None and self.charge != 0

    @property
    def is_protonated(self) -> bool:
        """Check if the formula is protonated (charge == 1)."""
        if len(self.formula) == 1:
            if self.formula[0].element == Element.H:
                return True
        return False


@dataclass(frozen=True, slots=True)
class PositionRule:
    terminal: Terminal
    amino_acid: AminoAcid | None = None

    @staticmethod
    def from_string(s: str) -> PositionRule:
        from ..proforma_components.parsers import parse_position_rule

        return parse_position_rule(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_position_rule

        return serialize_position_rule(self)

    def __str__(self) -> str:
        return self.serialize()

    def to_dict(self) -> dict[str, Any]:
        """Return the versioned, JSON-compatible representation of this rule."""
        from ..proforma_json import to_proforma_dict

        return to_proforma_dict(self)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Restore this rule from its versioned JSON-compatible representation."""
        from ..proforma_json import from_proforma_dict

        return from_proforma_dict(data, expected_type=cls)

    def to_json(self, *, indent: int | None = None) -> str:
        """Return deterministic JSON text for this rule."""
        from ..proforma_json import to_proforma_json

        return to_proforma_json(self, indent=indent)

    @classmethod
    def from_json(cls, data: str | bytes | bytearray) -> Self:
        """Restore this rule from versioned JSON text."""
        from ..proforma_json import from_proforma_json

        return from_proforma_json(data, expected_type=cls)


@dataclass(frozen=True, slots=True)
class TagAccession(MassPropertyMixin, PositionScoreMixin):
    """The accession for a modification"""

    accession: str
    cv: CV
    position_id: str | None = None
    score: float | None = None

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        try:
            mod_info = self._get_mod_info_by_accession()
            if mod_info is None:
                return f"Unknown accession: {self.accession} for CV: {self.cv}"
            return None
        except Exception as e:
            return str(e)

    def _get_mod_info_by_accession(
        self,
    ) -> UnimodInfo | PsimodInfo | ResidInfo | GnoInfo | XlModInfo | None:
        match self.cv:
            case CV.UNIMOD:
                return UNIMOD_LOOKUP.query_id(self.accession)
            case CV.PSI_MOD:
                return PSIMOD_LOOKUP.query_id(self.accession)
            case CV.RESID:
                return RESID_LOOKUP.query_id(self.accession)
            case CV.GNOME:
                return GNO_LOOKUP.query_id(self.accession)
            case CV.XL_MOD:
                return XLMOD_LOOKUP.query_id(self.accession)
            case _:
                raise ValueError(f"Modification lookup by accession not implemented for CV: {self.cv}")

        return None

    def get_mass(self, monoisotopic: bool = True) -> float:
        mod_info = self._get_mod_info_by_accession()
        if mod_info is not None:
            mass = mod_info.monoisotopic_mass if monoisotopic else mod_info.average_mass
            if mass is None:
                kind = "monoisotopic" if monoisotopic else "average"
                raise ValueError(f"Modification '{self}' was found but has no {kind} mass in its controlled vocabulary.")
            return mass
        raise ValueError(f"Unknown modification accession '{self}': not found in the '{self.cv}' controlled vocabulary. {_MOD_SPEC_HINT}")

    def get_charge(self) -> int | None:
        return None

    def get_composition(self) -> Counter[ElementInfo]:
        mod_info = self._get_mod_info_by_accession()
        if mod_info is not None:
            comp = mod_info.composition
            if comp is None:
                raise ValueError(f"Modification '{self}' was found but has no elemental composition in its controlled vocabulary.")
            return Counter(comp)
        raise ValueError(f"Unknown modification accession '{self}': not found in the '{self.cv}' controlled vocabulary. {_MOD_SPEC_HINT}")

    @staticmethod
    def from_string(s: str) -> TagAccession:
        from ..proforma_components.parsers import parse_tag_accession

        return parse_tag_accession(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_tag_accession

        return serialize_tag_accession(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class TagMass(MassPropertyMixin, PositionScoreMixin):
    """A mass modification"""

    mass_str: str | float
    cv: CV | None = None
    position_id: str | None = None
    score: float | None = None

    def __post_init__(self):
        # Ensure mass_str includes +/- sign
        mass_str = str(self.mass_str)

        # Add sign if not present
        if not mass_str.startswith(("+", "-")):
            mass_str = f"+{mass_str}"

        object.__setattr__(self, "mass_str", mass_str)

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        try:
            _ = self.get_mass()
            return None
        except Exception as e:
            return str(e)

    @property
    def mass(self) -> float:
        return float(self.mass_str)

    def get_mass(self, monoisotopic: bool = True) -> float:
        return self.mass

    def get_composition(self) -> Counter[ElementInfo]:
        match self.cv:
            case CV.UNIMOD:
                mod_info = UNIMOD_LOOKUP.query_mass(self.mass, monoisotopic=True, tolerance=0.005)
            case CV.PSI_MOD:
                mod_info = PSIMOD_LOOKUP.query_mass(self.mass, monoisotopic=True, tolerance=0.005)
            case CV.RESID:
                mod_info = RESID_LOOKUP.query_mass(self.mass, monoisotopic=True, tolerance=0.005)
            case CV.GNOME:
                mod_info = GNO_LOOKUP.query_mass(self.mass, monoisotopic=True, tolerance=0.005)
            case CV.XL_MOD:
                mod_info = XLMOD_LOOKUP.query_mass(self.mass, monoisotopic=True, tolerance=0.005)
            case _:
                raise ValueError(f"Modification lookup by mass not implemented for CV: {self.cv}")

        if len(mod_info) > 1:
            raise ValueError(f"Multiple modifications found for mass: {self.mass} in CV: {self.cv}")

        return Counter(mod_info[0].composition) if mod_info and mod_info[0].composition else Counter()

    def get_charge(self) -> int | None:
        return None

    @staticmethod
    def from_string(s: str) -> TagMass:
        from ..proforma_components.parsers import parse_tag_mass

        return parse_tag_mass(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_tag_mass

        return serialize_tag_mass(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class PositionScore(MassPropertyMixin, PositionScoreMixin):
    position_id: str
    score: float | None = None

    def get_mass(self, monoisotopic: bool = True) -> float:
        return 0.0

    def get_composition(self) -> Counter[ElementInfo]:
        return Counter()

    def validate(self) -> str | None:
        return None


@dataclass(frozen=True, slots=True)
class TagName(MassPropertyMixin, PositionScoreMixin):
    """A named modification"""

    name: str
    cv: CV | None = None
    position_id: str | None = None
    score: float | None = None

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        try:
            mod_info = self._get_mod_info_by_name()
            if mod_info is None:
                return f"Unknown modification name: {self.name}"
            return None
        except Exception as e:
            return str(e)

    def _get_mod_info_by_name(
        self,
    ) -> UnimodInfo | PsimodInfo | ResidInfo | GnoInfo | XlModInfo | None:
        match self.cv:
            case CV.UNIMOD:
                return UNIMOD_LOOKUP.query_name(self.name)
            case CV.PSI_MOD:
                return PSIMOD_LOOKUP.query_name(self.name)
            case CV.RESID:
                return RESID_LOOKUP.query_name(self.name)
            case CV.GNOME:
                return GNO_LOOKUP.query_name(self.name)
            case CV.XL_MOD:
                return XLMOD_LOOKUP.query_name(self.name)
            case None:
                unimod = UNIMOD_LOOKUP.query_name(self.name)
                if unimod is not None:
                    return unimod
                psimod = PSIMOD_LOOKUP.query_name(self.name)
                if psimod is not None:
                    return psimod
            case _:
                raise ValueError(f"Modification lookup by name not implemented for CV: {self.cv}")

        return None

    def get_mass(self, monoisotopic: bool = True) -> float:
        mod_info = self._get_mod_info_by_name()
        if mod_info is not None:
            mass = mod_info.monoisotopic_mass if monoisotopic else mod_info.average_mass
            if mass is None:
                kind = "monoisotopic" if monoisotopic else "average"
                raise ValueError(f"Modification '{self}' was found but has no {kind} mass in its controlled vocabulary.")
            return mass
        raise ValueError(f"Unknown modification name '{self}': not found in any controlled vocabulary. {_MOD_SPEC_HINT}")

    def get_composition(self) -> Counter[ElementInfo]:
        mod_info = self._get_mod_info_by_name()
        if mod_info is not None:
            comp = mod_info.composition
            if comp is None:
                raise ValueError(f"Modification '{self}' was found but has no elemental composition in its controlled vocabulary.")
            return Counter(comp)
        raise ValueError(f"Unknown modification name '{self}': not found in any controlled vocabulary. {_MOD_SPEC_HINT}")

    def get_charge(self) -> int | None:
        return None

    @staticmethod
    def from_string(s: str) -> TagName:
        from ..proforma_components.parsers import parse_tag_name

        return parse_tag_name(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_tag_name

        return serialize_tag_name(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True)
class TagInfo(MassPropertyMixin, PositionScoreMixin):
    """An INFO tag modification"""

    info: str

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        return None

    def get_mass(self, monoisotopic: bool = True) -> float:
        return 0.0

    def get_composition(self) -> Counter[ElementInfo]:
        return Counter()

    def get_charge(self) -> int | None:
        return None

    @staticmethod
    def from_string(s: str) -> TagInfo:
        from ..proforma_components.parsers import parse_tag_info

        return parse_tag_info(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_tag_info

        return serialize_tag_info(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True)
class TagCustom(MassPropertyMixin, PositionScoreMixin):
    """A custom/user-defined modification"""

    name: str
    position_id: str | None = None
    score: float | None = None

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        return None

    def get_mass(self, monoisotopic: bool = True) -> float:
        return 0.0

    def get_composition(self) -> Counter[ElementInfo]:
        return Counter()

    def get_charge(self) -> int | None:
        return None

    @staticmethod
    def from_string(s: str) -> TagCustom:
        from ..proforma_components.parsers import parse_tag_custom

        return parse_tag_custom(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_tag_custom

        return serialize_tag_custom(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class GlycanComponent(MassPropertyMixin):
    """A single component of a glycan composition.

    A component is a named monosaccharide, or -- for components not in the supported list --
    a molecular formula or a monoisotopic mass wrapped in curly braces (ProForma 2.1 §10.2),
    e.g. ``{C8H13N1O5}1``, ``{C8H13N1O5Na1:z+1}1`` (a charged formula) or ``{+203.079}1`` (a
    bare mass). A mass component contributes mass but has no elemental composition.
    """

    monosaccharide: Monosaccharide | ChargedFormula | float
    occurance: int

    @property
    def is_mass(self) -> bool:
        """True when this component is a bare monoisotopic mass (no elemental composition)."""
        return isinstance(self.monosaccharide, (int, float))

    def validate(self) -> str | None:
        try:
            _ = self.get_mass()
            return None
        except Exception as e:
            return str(e)

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def get_mass(self, monoisotopic: bool = True) -> float:
        value = self.monosaccharide
        if isinstance(value, (int, float)):
            return value * self.occurance
        elif isinstance(value, ChargedFormula):
            return value.get_mass(monoisotopic=monoisotopic) * self.occurance
        else:
            monosaccharide = MONOSACCHARIDE_LOOKUP.proforma(value)
            mass = monosaccharide.mass(monoisotopic=monoisotopic)
            if mass is None:
                raise ValueError(f"Unknown mass for monosaccharide: {value}")
            return mass * self.occurance

    def get_composition(self) -> Counter[ElementInfo]:
        # Must multiply by occurance to match get_mass (e.g. Glycan:Hex3 is three Hex units).
        value = self.monosaccharide
        if isinstance(value, (int, float)):
            # A bare mass has no elemental composition; callers that need the whole glycan's
            # composition route this through GlycanTag.get_composition_and_delta_mass instead.
            raise ValueError(f"Glycan mass component {{{value:+f}}} has no elemental composition")
        elif isinstance(value, ChargedFormula):
            composition = value.get_composition()
        else:
            monosaccharide = MONOSACCHARIDE_LOOKUP.proforma(value)
            comp = monosaccharide.composition
            if comp is None:
                raise ValueError(f"Unknown composition for monosaccharide: {value}")
            composition = Counter(comp)
        if self.occurance != 1:
            composition = Counter({element: count * self.occurance for element, count in composition.items()})
        return composition

    def get_charge(self) -> int | None:
        if isinstance(self.monosaccharide, ChargedFormula) and self.monosaccharide.charge is not None:
            return self.monosaccharide.charge * self.occurance
        return None

    @staticmethod
    def from_string(s: str) -> GlycanComponent:
        from ..proforma_components.parsers import parse_glycan_component

        return parse_glycan_component(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_glycan_component

        return serialize_glycan_component(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class GlycanTag(MassPropertyMixin, PositionScoreMixin):
    """A glycan composition tag"""

    components: tuple[GlycanComponent, ...]
    position_id: str | None = None
    score: float | None = None

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        try:
            _ = self.get_mass()
            return None
        except Exception as e:
            return str(e)

    def get_mass(self, monoisotopic: bool = True) -> float:
        return sum_masses(self.components, monoisotopic=monoisotopic)

    def get_composition(self) -> Counter[ElementInfo]:
        return merge_compositions(self.components)

    def get_composition_and_delta_mass(self, monoisotopic: bool = True) -> tuple[Counter[ElementInfo], float]:
        """Split this glycan into an elemental composition plus a residual delta mass.

        Formula and named-monosaccharide components contribute to the composition; bare-mass
        components (e.g. ``{+203.079}``) have no composition and contribute to the delta mass.
        Used by the composition path so a glycan that mixes the two still resolves.

        :param monoisotopic: Use monoisotopic masses when ``True``, average masses otherwise.
        :type monoisotopic: bool
        :return: A ``(composition, delta_mass)`` tuple.
        :rtype: tuple[Counter[ElementInfo], float]
        """
        composition: Counter[ElementInfo] = Counter()
        delta_mass = 0.0
        for component in self.components:
            if component.is_mass:
                delta_mass += component.get_mass(monoisotopic=monoisotopic)
            else:
                add_composition(composition, component.get_composition())
        return composition, delta_mass

    def get_charge(self) -> int | None:
        total = sum(component.get_charge() or 0 for component in self.components)
        return total or None

    def __len__(self) -> int:
        """Get the number of components in this glycan tag"""
        return len(self.components)

    def __getitem__(self, index: int) -> GlycanComponent:
        """Get a component by index"""
        return self.components[index]

    @staticmethod
    def from_string(s: str) -> GlycanTag:
        from ..proforma_components.parsers import parse_glycan

        return GlycanTag(parse_glycan(s))

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_glycan_tag

        return serialize_glycan_tag(self)

    def __str__(self) -> str:
        return self.serialize()


class PlacementTagMixin(MassPropertyMixin):
    """Mixin to add mass properties to classes that implement get_mass()"""

    def get_mass(self, monoisotopic: bool = True) -> float:
        raise ValueError("PlacementTag has no mass")

    def get_composition(self) -> Counter[ElementInfo]:
        return Counter()

    def get_charge(self) -> int | None:
        return None

    def validate(self) -> str | None:
        return None


@dataclass(frozen=True, slots=True)
class PositionTag(PlacementTagMixin):
    residues: tuple[PositionRule, ...]

    @staticmethod
    def from_string(s: str) -> PositionTag:
        s = s.strip().lower()
        if s.startswith("position:"):
            s = s[len("position:") :]
        else:
            raise ValueError("PositionTag string must start with 'Position:'")
        return PositionTag(tuple(PositionRule.from_string(part.strip()) for part in s.split(",")))

    def serialize(self) -> str:
        return "Position:" + ",".join(str(r) for r in self.residues)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class LimitTag(PlacementTagMixin):
    limit: int

    @staticmethod
    def from_string(s: str) -> LimitTag:
        s = s.strip().lower()
        if s.startswith("limit:"):
            s = s[len("limit:") :]
        else:
            raise ValueError("LimitTag string must start with 'Limit:'")
        return LimitTag(limit=int(s.strip()))

    def serialize(self) -> str:
        return f"Limit:{self.limit}"

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class ComkpTag(PlacementTagMixin):
    @staticmethod
    def from_string(s: str) -> ComkpTag:
        if not s.strip().lower() == "comkp":
            raise ValueError("ComkpTag string must be 'Comkp'")
        return ComkpTag()

    def serialize(self) -> str:
        return "CoMKP"

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class ComupTag(PlacementTagMixin):
    @staticmethod
    def from_string(s: str) -> ComupTag:
        if not s.strip().lower() == "comup":
            raise ValueError("ComupTag string must be 'Comup'")
        return ComupTag()

    def serialize(self) -> str:
        return "CoMUP"

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True)
class IsotopeReplacement(MassPropertyMixin):
    """A global isotope replacement"""

    element: Element
    isotope: int

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        try:
            _ = self.get_isotope_replacements()
            return None
        except Exception as e:
            return str(e)

    @staticmethod
    def from_string(s: str) -> IsotopeReplacement:
        from ..proforma_components.parsers import parse_isotope_replacement

        return parse_isotope_replacement(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_isotope_replacement

        return serialize_isotope_replacement(self)

    def get_isotope_replacements(self) -> tuple[ElementInfo, ElementInfo]:
        """Get a tuple of (original ElementInfo, replaced ElementInfo)"""
        original = ELEMENT_LOOKUP[(self.element, None)]
        replaced = ELEMENT_LOOKUP[(self.element, self.isotope)]
        return (original, replaced)

    def __str__(self) -> str:
        return self.serialize()

    def get_mass(self, monoisotopic: bool = True) -> float:
        return 0.0

    def get_composition(self) -> Counter[ElementInfo]:
        return Counter()

    def get_charge(self) -> int | None:
        return None


@dataclass(frozen=True)
class GlobalChargeCarrier(MassPropertyMixin):
    """A charge carrier specification, a bare charged formula like 'Na:z+1' or 'H:z+1^2'.

    Per ProForma 2.1 section 11.5 charge carriers are written without a ``Formula:`` prefix
    (that prefix is only used for localised residue modifications).
    """

    charged_formula: ChargedFormula
    occurance: int

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        return self.charged_formula.validate()

    def get_mass(self, monoisotopic: bool = True) -> float:
        return self.charged_formula.get_mass(monoisotopic=monoisotopic) * self.occurance

    def get_charge(self) -> int:
        if self.charged_formula.charge is None:
            raise ValueError("Charge carrier has no defined charge")
        return int(self.charged_formula.charge * self.occurance)

    def get_composition(self) -> Counter[ElementInfo]:
        composition = self.charged_formula.get_composition()
        # Multiply all counts by occurance
        for elem_info, count in composition.items():
            composition[elem_info] = int(count * self.occurance)
        return composition

    @staticmethod
    def from_string(s: str) -> GlobalChargeCarrier:
        from ..proforma_components.parsers import parse_global_charge_carrier

        return parse_global_charge_carrier(s)

    @staticmethod
    def charged_proton(charge: int) -> GlobalChargeCarrier:
        """Create a GlobalChargeCarrier representing protons."""
        return GlobalChargeCarrier(charged_formula=PROTON_FORMULA, occurance=charge)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_global_charge_carrier

        return serialize_global_charge_carrier(self)

    def __str__(self) -> str:
        return self.serialize()

    def to_mz_paf(self) -> str:
        """Convert to mzPAF format string."""
        # mzPAF adduct notation (spec section 4.7) prefixes a repeated adduct with its
        # count, e.g. "[M+2Na]" for two sodium atoms; a count of 1 is omitted. This is
        # self.occurance -- how many instances of this charge carrier are present --
        # not to be confused with a count baked into charged_formula's own elements.
        #
        # The +/- direction of the adduct is the combination of TWO signs: the charged
        # formula's own sign (e.g. a removed proton is stored as "H-1", serialized "-H")
        # and the sign of occurance (a negative occurance, e.g. charged_proton(-2) for a
        # doubly deprotonated ion, flips the direction). Deriving the sign purely from the
        # formula and pasting str(occurance) after it produced malformed output like
        # "M+-2H"; XOR-ing the two signs and using the magnitude gives "M-2H".
        paf_formula = self.charged_formula.to_mz_paf()
        formula_sign, rest = paf_formula[0], paf_formula[1:]
        negative = (formula_sign == "-") ^ (self.occurance < 0)
        count = abs(self.occurance)
        count_str = str(count) if count != 1 else ""
        return f"M{'-' if negative else '+'}{count_str}{rest}"

    @property
    def is_protonated(self) -> bool:
        if self.charged_formula.is_protonated:
            return True
        return False


# Type alias for modification tags
MODIFICATION_TAG_TYPE = (
    TagAccession | ChargedFormula | GlycanTag | TagInfo | TagMass | TagName | TagCustom | PositionScore | PositionTag | LimitTag | ComkpTag | ComupTag
)


@dataclass(frozen=True, slots=True)
class ModificationTags(MassPropertyMixin):
    """A modification"""

    tags: tuple[MODIFICATION_TAG_TYPE, ...]

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self, all_tags: bool = False) -> str | None:
        if not self.tags:
            return "ModificationTags cannot be empty"

        tags_to_check = self.tags if all_tags else (self.tags[0],)

        for tag in tags_to_check:
            if hasattr(tag, "validate"):
                error = tag.validate()
                if error is not None:
                    return error

        return None

    @property
    def first_tag(self) -> MODIFICATION_TAG_TYPE:
        """Get the first tag in this modification"""
        return self.tags[0]

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get the mass from this modification"""
        return self.first_tag.get_mass(monoisotopic=monoisotopic)

    def get_composition(self) -> Counter[ElementInfo]:
        """Get the composition from this modification"""
        return self.first_tag.get_composition()

    def get_charge(self) -> int | None:
        """Get the charge of this modification, if any."""
        if isinstance(self.first_tag, ChargedFormula):
            return self.first_tag.charge
        if isinstance(self.first_tag, GlycanTag):
            # A glycan may carry a charged formula component (e.g. {C8H13N1O5Na1:z+1}).
            return self.first_tag.get_charge()
        return None

    def __len__(self) -> int:
        """Get the number of tags in this modification"""
        return len(self.tags)

    def __getitem__(self, index: int) -> MODIFICATION_TAG_TYPE:
        """Get a tag by index"""
        return self.tags[index]

    @staticmethod
    def from_string(s: str) -> ModificationTags:
        from ..proforma_components.parsers import parse_modification_tags

        return parse_modification_tags(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_modification_tags

        return serialize_modification_tags(self)

    def __str__(self) -> str:
        return self.serialize()

    @property
    def has_placement_tags(self) -> bool:
        """Check if any tags specify placement rules (PositionTag, LimitTag, ComkpTag, ComupTag)"""
        for tag in self.tags:
            if isinstance(tag, (PositionTag, LimitTag, ComkpTag, ComupTag)):
                return True
        return False

    @property
    def placement_tags(
        self,
    ) -> tuple[PositionTag | LimitTag | ComkpTag | ComupTag, ...]:
        """Get all placement tags (PositionTag, LimitTag, ComkpTag, ComupTag)"""
        placement_tags: list[PositionTag | LimitTag | ComkpTag | ComupTag] = []
        for tag in self.tags:
            if isinstance(tag, (PositionTag, LimitTag, ComkpTag, ComupTag)):
                placement_tags.append(tag)
        return tuple(placement_tags)

    @property
    def position_tag(self) -> PositionTag | None:
        """Get the PositionTag if present, else None"""
        for tag in self.tags:
            if isinstance(tag, PositionTag):
                return tag
        return None

    @property
    def limit_tag(self) -> LimitTag | None:
        """Get the LimitTag if present, else None"""
        for tag in self.tags:
            if isinstance(tag, LimitTag):
                return tag
        return None

    @property
    def comkp_tag(self) -> ComkpTag | None:
        """Get the ComkpTag if present, else None"""
        for tag in self.tags:
            if isinstance(tag, ComkpTag):
                return tag
        return None

    @property
    def comup_tag(self) -> ComupTag | None:
        """Get the ComupTag if present, else None"""
        for tag in self.tags:
            if isinstance(tag, ComupTag):
                return tag
        return None


@dataclass(frozen=True, slots=True)
class ModificationAmbiguousPrimary(MassPropertyMixin):
    """The primary definition of an ambiguous modification"""

    label: str
    tags: ModificationTags
    score: float | None = None
    position: tuple[PositionRule, ...] | None = None
    limit: int | None = None
    comkp: bool | None = None
    comup: bool | None = None

    def __post_init__(self):
        """Validate constraints."""
        if self.score is not None and not (0 <= self.score <= 1):
            raise ValueError(f"Score must be between 0 and 1, got {self.score}")

        if self.limit is not None and self.limit < 1:
            raise ValueError(f"Limit must be positive, got {self.limit}")

        if not self.tags:
            raise ValueError("tags cannot be empty")

    def get_mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError()

    def get_composition(self) -> Counter[ElementInfo]:
        raise NotImplementedError()

    def __len__(self) -> int:
        """Get the number of tags in this modification"""
        return len(self.tags)

    @staticmethod
    def from_string(s: str) -> ModificationAmbiguousPrimary:
        from ..proforma_components.parsers import parse_modification_ambiguous_primary

        return parse_modification_ambiguous_primary(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import (
            serialize_modification_ambiguous_primary,
        )

        return serialize_modification_ambiguous_primary(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class ModificationAmbiguousSecondary(MassPropertyMixin):
    """A reference to an ambiguous modification"""

    label: str
    score: float | None = None

    def __post_init__(self):
        """Validate constraints."""
        if self.score is not None and not (0 <= self.score <= 1):
            raise ValueError(f"Score must be between 0 and 1, got {self.score}")

    def get_mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError()

    def get_composition(self) -> Counter[ElementInfo]:
        raise NotImplementedError()

    @staticmethod
    def from_string(s: str) -> ModificationAmbiguousSecondary:
        from ..proforma_components.parsers import parse_modification_ambiguous_secondary

        return parse_modification_ambiguous_secondary(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import (
            serialize_modification_ambiguous_secondary,
        )

        return serialize_modification_ambiguous_secondary(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class ModificationCrossLinker(MassPropertyMixin):
    """A cross-linked modification"""

    label: str | None = None
    tags: ModificationTags | None = None

    def get_mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError()

    def get_composition(self) -> Counter[ElementInfo]:
        raise NotImplementedError()

    def __len__(self) -> int:
        """Get the number of tags in this modification"""
        if self.tags is None:
            return 0
        return len(self.tags)

    @staticmethod
    def from_string(s: str) -> ModificationCrossLinker:
        from ..proforma_components.parsers import parse_modification_cross_linker

        return parse_modification_cross_linker(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import (
            serialize_modification_cross_linker,
        )

        return serialize_modification_cross_linker(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class FixedModification(MassPropertyMixin):
    """A fixed modification"""

    modifications: ModificationTags
    position_rules: tuple[PositionRule, ...] = ()

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        return self.modifications.validate()

    def get_mass(self, monoisotopic: bool = True) -> float:
        return self.modifications.get_mass(monoisotopic=monoisotopic)

    def get_composition(self) -> Counter[ElementInfo]:
        return self.modifications.get_composition()

    def get_charge(self) -> int | None:
        """Get the charge of this modification, if any."""
        return self.modifications.get_charge()

    def __len__(self) -> int:
        """Get the number of tags in this modification"""
        return len(self.modifications)

    def find_indexes(self, sequence: str) -> list[int]:
        """Find all indexes in the sequence where this fixed modification applies."""
        # can be N [-1], C [-2] or internal (int)
        indexes: list[int] = []
        for rule in self.position_rules:
            match rule.terminal:
                case Terminal.N_TERM:
                    if rule.amino_acid is None or sequence[0] == rule.amino_acid:
                        indexes.append(-1)
                case Terminal.C_TERM:
                    if rule.amino_acid is None or sequence[-1] == rule.amino_acid:
                        indexes.append(-2)
                case Terminal.ANYWHERE:
                    for i, aa in enumerate(sequence):
                        if rule.amino_acid is None or aa == rule.amino_acid:
                            indexes.append(i)

        return indexes

    @staticmethod
    def from_string(s: str) -> FixedModification:
        from ..proforma_components.parsers import parse_fixed_modification

        return parse_fixed_modification(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_fixed_modification

        return serialize_fixed_modification(self)

    def __str__(self) -> str:
        return self.serialize()


# Type aliases for modification types (defined after all classes)
MODIFICATION_AMBIGUOUS_TYPE = ModificationAmbiguousPrimary | ModificationAmbiguousSecondary
MODIFICATION_TYPE = MODIFICATION_AMBIGUOUS_TYPE | ModificationCrossLinker | ModificationTags


@dataclass(frozen=True, slots=True)
class SequenceElement(MassPropertyMixin):
    """A single amino acid with optional modifications like 'M[Oxidation]' or 'K'"""

    amino_acid: AminoAcid
    modifications: tuple[MODIFICATION_TYPE, ...] = ()

    def get_mass(self, monoisotopic: bool = True) -> float:
        aa = AA_LOOKUP.one_letter(self.amino_acid)
        aa_mass = aa.monoisotopic_mass if monoisotopic else aa.average_mass
        if aa_mass is None:
            raise ValueError(f"Unknown mass for amino acid: {self.amino_acid}")
        mod_mass = sum_masses(self.modifications, monoisotopic=monoisotopic)
        return aa_mass + mod_mass

    def get_composition(self) -> Counter[ElementInfo]:
        aa = AA_LOOKUP.one_letter(self.amino_acid)
        composition = aa.composition
        if composition is None:
            raise ValueError(f"Unknown composition for amino acid: {self.amino_acid}")

        total_composition = Counter(composition)
        if self.modifications:
            add_composition(total_composition, merge_compositions(self.modifications))

        return total_composition

    @staticmethod
    def from_string(s: str) -> SequenceElement:
        from ..proforma_components.parsers import parse_sequence_element

        return parse_sequence_element(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_sequence_element

        return serialize_sequence_element(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class SequenceRegion(MassPropertyMixin):
    """A region of sequence with modifications"""

    sequence: tuple[SequenceElement, ...]
    modifications: tuple[MODIFICATION_TYPE, ...]
    ambiguous: bool

    def get_mass(self, monoisotopic: bool = True) -> float:
        return sum_masses(self.sequence + self.modifications, monoisotopic=monoisotopic)

    def get_composition(self) -> Counter[ElementInfo]:
        return merge_compositions(self.sequence + self.modifications)

    @staticmethod
    def from_string(s: str) -> SequenceRegion:
        from ..proforma_components.parsers import parse_sequence_region

        return parse_sequence_region(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_sequence_region

        return serialize_sequence_region(self)

    def __str__(self) -> str:
        return self.serialize()


SEQUENCE_TYPE = SequenceElement | SequenceRegion

# Type alias for global charge
GLOBAL_CHARGE_TYPE = int | tuple[GlobalChargeCarrier, ...]


@dataclass(frozen=True, slots=True)
class Peptidoform(MassPropertyMixin):
    """A Peptidoform"""

    sequence: tuple[SEQUENCE_TYPE, ...]
    name: str | None = None
    n_term_modifications: tuple[MODIFICATION_TYPE, ...] = ()
    c_term_modifications: tuple[MODIFICATION_TYPE, ...] = ()
    labile_modifications: tuple[ModificationTags, ...] = ()
    unlocalised_modifications: tuple[MODIFICATION_AMBIGUOUS_TYPE, ...] = ()

    def get_mass(self, monoisotopic: bool = True) -> float:
        return sum_masses(
            self.sequence + self.n_term_modifications + self.c_term_modifications + self.labile_modifications + self.unlocalised_modifications,
            monoisotopic=monoisotopic,
        )

    def get_composition(self) -> Counter[ElementInfo]:
        return merge_compositions(
            self.sequence + self.n_term_modifications + self.c_term_modifications + self.labile_modifications + self.unlocalised_modifications,
        )

    @staticmethod
    def from_string(s: str) -> Peptidoform:
        from ..proforma_components.parsers import parse_peptidoform

        return parse_peptidoform(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_peptidoform

        return serialize_peptidoform(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class PeptidoformIon(MassPropertyMixin):
    """A peptidoform ion"""

    peptidoforms: tuple[Peptidoform, ...]
    name: str | None = None
    charge: GLOBAL_CHARGE_TYPE | None = None

    def get_mass(self, monoisotopic: bool = True) -> float:
        """
        mass = sum_masses(self.peptidoforms, monoisotopic=monoisotopic)

        if isinstance(self.charge, int) and self.charge != 0:
            charge_mass = PROTON_MASS if self.charge > 0 else -ELECTRON_MASS
            mass += self.charge * charge_mass
        elif isinstance(self.charge, tuple):
            mass += sum_masses(self.charge, monoisotopic=monoisotopic)

        return mass
        """
        raise NotImplementedError()

    def get_composition(self) -> Counter[ElementInfo]:
        """
        comp: Counter[ElementInfo] = merge_compositions(self.peptidoforms)

        if isinstance(self.charge, int) and self.charge > 0:
            hydrogen_info = ELEMENT_LOOKUP["H"]
            comp[hydrogen_info] += self.charge  # Counter handles missing keys
        elif isinstance(self.charge, tuple):
            comp += merge_compositions(self.charge)

        return comp
        """
        raise NotImplementedError()

    @staticmethod
    def from_string(s: str) -> PeptidoformIon:
        from ..proforma_components.parsers import parse_peptidoform_ion

        return parse_peptidoform_ion(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_peptidoform_ion

        return serialize_peptidoform_ion(self)

    def __str__(self) -> str:
        return self.serialize()


@dataclass(frozen=True, slots=True)
class CompoundPeptidoformIon(MassPropertyMixin):
    """Describe compound peptidoform ion attributes"""

    peptidoform_ions: tuple[PeptidoformIon, ...]
    name: str | None = None
    fixed_modifications: tuple[FixedModification, ...] = ()
    isotope_replacement: tuple[IsotopeReplacement, ...] = ()

    def get_mass(self, monoisotopic: bool = True) -> float:
        raise NotImplementedError()
        total_mass: float = sum_masses(self.peptidoform_ions, monoisotopic=monoisotopic)

        # if fixed or isotope modifications raise NotImplementedError
        if self.fixed_modifications or self.isotope_replacement:
            raise NotImplementedError("Fixed modifications and isotope replacements are not yet supported in mass calculation.")

        return total_mass

    def get_composition(self) -> Counter[ElementInfo] | None:
        raise NotImplementedError()
        total_composition: Counter[ElementInfo] = merge_compositions(self.peptidoform_ions)

        # if fixed or isotope modifications raise NotImplementedError
        if self.fixed_modifications or self.isotope_replacement:
            raise NotImplementedError("Fixed modifications and isotope replacements are not yet supported in composition calculation.")

        return total_composition

    @staticmethod
    def from_string(s: str) -> CompoundPeptidoformIon:
        from ..proforma_components.parsers import parse_compound_peptidoform_ion

        return parse_compound_peptidoform_ion(s)

    def serialize(self) -> str:
        from ..proforma_components.serializers import serialize_compound_peptidoform_ion

        return serialize_compound_peptidoform_ion(self)

    def __str__(self) -> str:
        return self.serialize()


# H:z+1
PROTON_FORMULA = ChargedFormula(formula=(FormulaElement(element=Element.H, occurance=1),), charge=1)

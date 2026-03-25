import logging
import re
from collections import Counter
from collections.abc import Callable, Generator, Iterable, Mapping, Sequence
from enum import StrEnum
from itertools import product
from typing import (
    Any,
    Self,
    cast,
)

from tacular import (
    AA_LOOKUP,
    FRAGMENT_ION_LOOKUP,
    NEUTRAL_DELTA_LOOKUP,
    AminoAcid,
    Element,
    ElementInfo,
    FragmentIonInfo,
    IonType,
    IonTypeLiteral,
    NeutralDelta,
    NeutralDeltaInfo,
    NeutralDeltaLiteral,
)

from ..constants import PROTON_MASS, ModType, ModTypeLiteral, Terminal
from ..digestion.core import (
    EnzymeConfig,
    digest_annotation_by_aa,
    digest_annotation_by_regex,
    generate_regex,
    get_cleavage_sites,
    left_semi_spans,
    nonspecific_spans,
    right_semi_spans,
    semi_spans,
    sequential_digest_annotation,
)
from ..isotope import (
    IsotopicData,
    estimate_isotopic_distribution,
    isotopic_distribution,
)
from ..proforma_components import (
    MODIFICATION_TYPE,
    ChargedFormula,
    FixedModification,
    FormulaElement,
    GlobalChargeCarrier,
    IsotopeReplacement,
    ModificationTags,
    PositionRule,
    SequenceElement,
    SequenceRegion,
    TagMass,
)
from ..property.prop import AnnotationProperties
from ..spans import Span
from ..utils import get_mods
from .ambiguity import (
    annotate_ambiguity,
    condense_ambiguity_to_xnotation,
    group_by_ambiguity,
    unique_fragments,
)
from .cached_comps import DeltaInfo, IsotopeInfo
from .combinatorics import (
    generate_combinations,
    generate_combinations_with_replacement,
    generate_permutations,
    generate_product,
)
from .manipulation import (
    condense_mods_to_intervals,
    condense_static_mods,
    condense_to_peptidoform,
    count_residues,
    coverage,
    find_indices,
    is_subsequence,
    modification_coverage,
    percent_coverage,
    percent_residues,
)
from .mod import (
    Interval,
    Mod,
    Mods,
    convert_moddict_input,
    convert_single_mod_input,
)
from .mod_builder import modify
from .parser import ProFormaParser
from .positions import validate_position
from .randomizer import generate_random_proforma_annotation
from .serializer import serialize_annotation, serialize_charge
from .slicing import (
    generate_sliding_windows,
    join_annotations,
    reverse_annotation,
    shift_annotation,
    shuffle_annotation,
    slice_annotation,
    sort_annotation,
    split_annotation,
)
from .utils import (
    Fragment,
    adjust_comp,
    adjust_mass_mz,
    can_fragment_sequence,
)

logger = logging.getLogger(__name__)

fe = FormulaElement(element=Element.H, occurance=1)
H_CHARGE_FORMULA = ChargedFormula(formula=(fe,), charge=1)


ION_TYPE = IonTypeLiteral | IonType
CHARGE_TYPE = int | str | list[str] | Mods[GlobalChargeCarrier] | GlobalChargeCarrier | Mod[GlobalChargeCarrier]
ISOTOPE_TYPE = int | dict[str | ElementInfo, int]
LOSS_TYPE = NeutralDelta | NeutralDeltaLiteral | NeutralDeltaInfo | str
CUSTOM_LOSS_TYPE = str | ChargedFormula | float | dict[str | ChargedFormula | float, int]
POSITION_TYPE = int | tuple[int, int]

EMPTY_ISOTOPE_MODS = Mods[IsotopeReplacement](mod_type=ModType.ISOTOPE, _mods=None)
EMPTY_STATIC_MODS = Mods[FixedModification](mod_type=ModType.STATIC, _mods=None)
EMPTY_UNKNOWN_MODS = Mods[ModificationTags](mod_type=ModType.UNKNOWN, _mods=None)
EMPTY_LABILE_MODS = Mods[ModificationTags](mod_type=ModType.LABILE, _mods=None)
EMPTY_NTERM_MODS = Mods[ModificationTags](mod_type=ModType.NTERM, _mods=None)
EMPTY_CTERM_MODS = Mods[ModificationTags](mod_type=ModType.CTERM, _mods=None)
EMPTY_CHARGE_MODS = Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods=None)
EMPTY_INTERNAL_MODS = Mods[ModificationTags](mod_type=ModType.INTERNAL, _mods=None)


class ChargeType(StrEnum):
    INT = "int"
    ADDUCTS = "adducts"
    NONE = "none"


def get_loss_combinations(losses: dict[NeutralDeltaInfo, int], max_losses: int) -> list[DeltaInfo]:
    """Generate all combinations of losses up to max_losses."""
    if not losses:
        return [DeltaInfo.from_input(None)]

    # Generate all possible count combinations for each loss
    loss_items = list(losses.items())
    count_ranges = [range(count + 1) for _, count in loss_items]

    loss_combinations: list[dict[NeutralDeltaInfo, int] | None] = [None]

    for counts in product(*count_ranges):
        total_losses = sum(counts)
        # Skip if no losses or exceeds max
        if total_losses == 0 or total_losses > max_losses:
            continue

        # Build combination dict
        combo = {}
        for (loss, _), count in zip(loss_items, counts, strict=True):
            if count > 0:
                combo[loss] = count

        loss_combinations.append(combo)

    # Convert to DeltaInfo
    delta_combinations: list[DeltaInfo] = []
    for loss_combo in loss_combinations:
        if loss_combo is None:
            delta_combinations.append(DeltaInfo.from_input(None))
        else:
            formula_dict: dict[ChargedFormula, int] = {}
            for nd, count in loss_combo.items():
                nd_formula: ChargedFormula = ChargedFormula.from_composition(nd.composition)
                formula_dict[nd_formula] = count
            delta_combinations.append(DeltaInfo.from_input(formula_dict))  # type: ignore

    return delta_combinations


class ProFormaAnnotation:
    def __init__(
        self,
        sequence: str | None = None,
        compound_name: str | None = None,  # (>>>Name)
        ion_name: str | None = None,  # (>>Name)
        peptide_name: str | None = None,  # (>Name)
        isotope_mods: Any = None,
        static_mods: Any = None,
        labile_mods: Any = None,
        unknown_mods: Any = None,
        nterm_mods: Any = None,
        cterm_mods: Any = None,
        internal_mods: dict[int, Any] | None = None,
        intervals: list[Interval] | None = None,
        charge: Any = None,
        validate: bool = False,
    ) -> None:
        """Construct a ProFormaAnnotation.

        All modification parameters accept flexible input types and are
        normalised internally via :func:`convert_moddict_input`.  Pass
        ``validate=True`` to enforce structural validity on construction;
        leave it ``False`` (the default) for performance when building
        annotations programmatically.  Prefer the :meth:`parse` factory
        when constructing from a ProForma string.

        :param sequence: Bare amino-acid sequence (single-letter codes).
        :type sequence: str | None
        :param compound_name: Compound-level name encoded as ``(>>>Name)``.
        :type compound_name: str | None
        :param ion_name: Ion-level name encoded as ``(>>Name)``.
        :type ion_name: str | None
        :param peptide_name: Peptide-level name encoded as ``(>Name)``.
        :type peptide_name: str | None
        :param isotope_mods: Global isotope replacement modifications (``<13C>`` style).
        :type isotope_mods: Any
        :param static_mods: Global fixed modifications (``<[Mod]@AA>`` style).
        :type static_mods: Any
        :param labile_mods: Labile modifications that may be lost during fragmentation.
        :type labile_mods: Any
        :param unknown_mods: Unknown-localisation modifications (``?[Mod]`` style).
        :type unknown_mods: Any
        :param nterm_mods: N-terminal modifications.
        :type nterm_mods: Any
        :param cterm_mods: C-terminal modifications.
        :type cterm_mods: Any
        :param internal_mods: Per-position modifications keyed by 0-based index.
        :type internal_mods: dict[int, Any] | None
        :param intervals: Ambiguous sequence intervals.
        :type intervals: list[Interval] | None
        :param charge: Charge state as an integer or list of adduct strings.
        :type charge: Any
        :param validate: If ``True``, validate each field immediately after setting it.
        :type validate: bool
        """
        self._sequence: str | None = None
        self._compound_name: str | None = None
        self._ion_name: str | None = None
        self._peptide_name: str | None = None
        self._isotope_mods: dict[str, int] | None = None
        self._static_mods: dict[str, int] | None = None
        self._labile_mods: dict[str, int] | None = None
        self._unknown_mods: dict[str, int] | None = None
        self._nterm_mods: dict[str, int] | None = None
        self._cterm_mods: dict[str, int] | None = None
        self._internal_mods: dict[int, dict[str, int]] | None = None
        self._intervals: list[Interval] | None = None
        self._charge: int | list[str] | None = None
        self._validate = validate

        self.set_sequence(sequence, inplace=True, validate=validate)
        self.set_compound_name(compound_name, inplace=True, validate=validate)
        self.set_ion_name(ion_name, inplace=True, validate=validate)
        self.set_peptide_name(peptide_name, inplace=True, validate=validate)
        self.set_isotope_mods(isotope_mods, inplace=True, validate=validate)
        self.set_static_mods(static_mods, inplace=True, validate=validate)
        self.set_labile_mods(labile_mods, inplace=True, validate=validate)
        self.set_unknown_mods(unknown_mods, inplace=True, validate=validate)
        self.set_nterm_mods(nterm_mods, inplace=True, validate=validate)
        self.set_cterm_mods(cterm_mods, inplace=True, validate=validate)
        self.set_internal_mods(internal_mods, inplace=True, validate=validate)
        self.set_intervals(intervals, inplace=True, validate=validate)
        self.set_charge(charge, inplace=True, validate=validate)

    @property
    def charge_type(self) -> ChargeType:
        """Return the representation style of the stored charge (integer, adducts, or none).

        :return: The charge representation type.
        :rtype: ChargeType
        """
        if isinstance(self._charge, int):
            if self._charge == 0:
                return ChargeType.NONE
            return ChargeType.INT
        elif isinstance(self._charge, list):
            return ChargeType.ADDUCTS
        else:
            return ChargeType.NONE

    """
    Validators
    """

    def validate_sequence(self) -> None:
        """Check that every residue in the sequence is a recognised amino acid.

        :raises ValueError: If an unrecognised amino acid code is found.
        """
        for aa in self.sequence:
            if aa not in AA_LOOKUP:
                raise ValueError(f"Invalid amino acid '{aa}' in sequence '{self.sequence}'")

    def validate_isotope_mods(self) -> None:
        """Check that all isotope modifications are structurally valid.

        :raises ValueError: If any isotope modification is invalid.
        """
        if errors := self.isotope_mods.validate():
            raise ValueError(f"Invalid isotope modifications: {errors}")

    def validate_static_mods(self) -> None:
        """Check that all static (fixed) modifications are structurally valid.

        :raises ValueError: If any static modification is invalid.
        """
        if errors := self.static_mods.validate():
            raise ValueError(f"Invalid static modifications: {errors}")

    def validate_labile_mods(self) -> None:
        """Check that all labile modifications are structurally valid.

        :raises ValueError: If any labile modification is invalid.
        """
        if errors := self.labile_mods.validate():
            raise ValueError(f"Invalid labile modifications: {errors}")

    def validate_unknown_mods(self) -> None:
        """Check that all unknown-localisation modifications are structurally valid.

        :raises ValueError: If any unknown modification is invalid.
        """
        if errors := self.unknown_mods.validate():
            raise ValueError(f"Invalid unknown modifications: {errors}")

    def validate_nterm_mods(self) -> None:
        """Check that all N-terminal modifications are structurally valid.

        :raises ValueError: If any N-terminal modification is invalid.
        """
        if errors := self.nterm_mods.validate():
            raise ValueError(f"Invalid N-terminal modifications: {errors}")

    def validate_cterm_mods(self) -> None:
        """Check that all C-terminal modifications are structurally valid.

        :raises ValueError: If any C-terminal modification is invalid.
        """
        if errors := self.cterm_mods.validate():
            raise ValueError(f"Invalid C-terminal modifications: {errors}")

    def validate_internal_mods(self) -> None:
        """Check that all internal (per-position) modifications are structurally valid.

        :raises ValueError: If any internal modification at any position is invalid.
        """
        for pos, mods in self.internal_mods.items():
            if errors := mods.validate():
                raise ValueError(f"Invalid internal modifications at position {pos}: {errors}")

    def validate_intervals(self) -> None:
        """Check that all intervals are valid, non-overlapping, and within sequence bounds.

        :raises ValueError: If any interval is invalid, intervals overlap, or an interval
            falls outside the sequence length.
        """
        intervals = self.intervals
        for interval in intervals:
            if errors := interval.validate():
                raise ValueError(f"Invalid interval: {errors}")

        # ensure no overlapping intervals
        sorted_intervals = sorted(intervals, key=lambda x: x.start)
        for i in range(1, len(sorted_intervals)):
            if sorted_intervals[i].start < sorted_intervals[i - 1].end:
                raise ValueError(f"Overlapping intervals detected: {sorted_intervals[i - 1]} and {sorted_intervals[i]}")

        # ensure that intervals dont start/end out of bounds
        seq_len = len(self.sequence) if self._sequence is not None else 0
        for interval in intervals:
            if interval.start < 0 or interval.end > seq_len:
                raise ValueError(f"Interval {interval} is out of bounds for sequence length {seq_len}")

    def validate_charge(self) -> None:
        """Check that the charge value is structurally valid.

        :raises ValueError: If the charge adducts are invalid or the charge type is
            unrecognised.
        """
        charge_type = self.charge_type

        match charge_type:
            case ChargeType.INT:
                pass
            case ChargeType.ADDUCTS:
                if errors := self.charge_adducts.validate():
                    raise ValueError(f"Invalid charge adducts: {errors}")
            case ChargeType.NONE:
                pass
            case _:
                raise ValueError(f"Invalid charge type: {charge_type}")

    def validate_annotation(self) -> None:
        """Run all individual validators in order; raises on the first error found.

        :raises ValueError: If any component of the annotation is structurally invalid.
        """
        self.validate_sequence()
        self.validate_isotope_mods()
        self.validate_static_mods()
        self.validate_labile_mods()
        self.validate_unknown_mods()
        self.validate_nterm_mods()
        self.validate_cterm_mods()
        self.validate_internal_mods()
        self.validate_intervals()
        self.validate_charge()

    @property
    def start_aa(self) -> str | None:
        if self.has_sequence:
            return self.sequence[0]
        return None

    @property
    def end_aa(self) -> str | None:
        if self.has_sequence:
            return self.sequence[-1]
        return None

    # ============================================================================
    # Properties
    # ============================================================================

    @property
    def sequence_elements(self) -> tuple[SequenceElement, ...]:
        if self._sequence is None:
            return ()
        return tuple(
            SequenceElement.from_string(f"{aa}{self.get_internal_mods_str_at_index(i)}" if self.has_internal_mods_at_index(i) else aa)
            for i, aa in enumerate(self.sequence)
        )

    @property
    def sequence_regions(self) -> tuple[SequenceRegion, ...]:
        if self._sequence is None or self.has_intervals is False:
            return ()

        sequence_elements = self.sequence_elements
        sequence_regions: list[SequenceRegion] = []
        for interval in self.intervals:
            region_elements = sequence_elements[interval.start : interval.end]

            mods: list[MODIFICATION_TYPE] = []
            if interval.has_mods:
                for mod, count in interval.mods.parse_items():
                    for _ in range(count):
                        mods.append(mod)

            sequence_regions.append(
                SequenceRegion(
                    sequence=region_elements,
                    modifications=tuple(mods),
                    ambiguous=interval.ambiguous,
                )
            )

        return tuple(sequence_regions)

    @property
    def sequence_elements_and_regions(
        self,
    ) -> tuple[SequenceElement | SequenceRegion, ...]:
        if self._sequence is None:
            return ()

        if self.has_intervals is False:
            return self.sequence_elements

        elements_and_regions: list[SequenceElement | SequenceRegion] = []
        seq_index = 0
        for interval in self.intervals:
            # Add sequence elements before interval
            while seq_index < interval.start:
                elements_and_regions.append(self.sequence_elements[seq_index])
                seq_index += 1

            # Add sequence region for interval
            region_elements = self.sequence_elements[interval.start : interval.end]

            mods: list[MODIFICATION_TYPE] = []
            if interval.has_mods:
                for mod, count in interval.mods.parse_items():
                    for _ in range(count):
                        mods.append(mod)

            elements_and_regions.append(
                SequenceRegion(
                    sequence=region_elements,
                    modifications=tuple(mods),
                    ambiguous=interval.ambiguous,
                )
            )
            seq_index = interval.end

        # Add remaining sequence elements after last interval
        while seq_index < len(self.sequence_elements):
            elements_and_regions.append(self.sequence_elements[seq_index])
            seq_index += 1

        return tuple(elements_and_regions)

    @property
    def sequence(self) -> str:
        """The bare amino-acid sequence; returns an empty string when unset.

        :rtype: str
        """
        return self._sequence if self._sequence is not None else ""

    @sequence.setter
    def sequence(self, value: str | None) -> None:
        self.set_sequence(value, inplace=True, validate=self._validate)

    @property
    def compound_name(self) -> str:
        """Compound-level name (``>>>Name`` prefix); empty string when unset.

        :rtype: str
        """
        return self._compound_name if self._compound_name is not None else ""

    @compound_name.setter
    def compound_name(self, value: Any | None) -> None:
        self.set_compound_name(value, inplace=True, validate=self._validate)

    @property
    def compound_name_str(self) -> str:
        """Serialised compound name including the ``(>>>...)`` wrapper, or empty string.

        :rtype: str
        """
        if self._compound_name is None:
            return ""
        return f"(>>>{self._compound_name})"

    @property
    def ion_name(self) -> str:
        """Ion-level name (``>>Name`` prefix); empty string when unset.

        :rtype: str
        """
        return self._ion_name if self._ion_name is not None else ""

    @ion_name.setter
    def ion_name(self, value: Any | None) -> None:
        self.set_ion_name(value, inplace=True, validate=self._validate)

    @property
    def ion_name_str(self) -> str:
        """Serialised ion name including the ``(>>...)`` wrapper, or empty string.

        :rtype: str
        """
        if self._ion_name is None:
            return ""
        return f"(>>{self._ion_name})"

    @property
    def peptide_name(self) -> str:
        """Peptide-level name (``>Name`` prefix); empty string when unset.

        :rtype: str
        """
        return self._peptide_name if self._peptide_name is not None else ""

    @peptide_name.setter
    def peptide_name(self, value: Any | None) -> None:
        self.set_peptide_name(value, inplace=True, validate=self._validate)

    @property
    def peptide_name_str(self) -> str:
        """Serialised peptide name including the ``(>...)`` wrapper, or empty string.

        :rtype: str
        """
        if self._peptide_name is None:
            return ""
        return f"(>{self._peptide_name})"

    @property
    def isotope_mods(self) -> Mods[IsotopeReplacement]:
        """Global isotope replacement modifications; returns an empty ``Mods`` when unset.

        :rtype: Mods[IsotopeReplacement]
        """
        if self._isotope_mods is None:
            return EMPTY_ISOTOPE_MODS

        return Mods[IsotopeReplacement](mod_type=ModType.ISOTOPE, _mods=self._isotope_mods)

    @isotope_mods.setter
    def isotope_mods(self, value: Any) -> None:
        self.set_isotope_mods(value, inplace=True, validate=self._validate)

    @property
    def isotope_mods_str(self) -> str:
        """Serialised isotope modifications string, or empty string when unset.

        :rtype: str
        """
        if self._isotope_mods is None:
            return ""
        return self.isotope_mods.serialize()

    @property
    def static_mods(self) -> Mods[FixedModification]:
        """Global fixed modifications; returns an empty ``Mods`` when unset.

        :rtype: Mods[FixedModification]
        """
        if self._static_mods is None:
            return EMPTY_STATIC_MODS

        return Mods[FixedModification](mod_type=ModType.STATIC, _mods=self._static_mods)

    @static_mods.setter
    def static_mods(self, value: Any) -> None:
        self.set_static_mods(value, inplace=True, validate=self._validate)

    @property
    def static_mods_str(self) -> str:
        """Serialised static modifications string, or empty string when unset.

        :rtype: str
        """
        if self._static_mods is None:
            return ""
        return self.static_mods.serialize()

    @property
    def labile_mods(self) -> Mods[ModificationTags]:
        """Labile modifications that may be lost during fragmentation; empty ``Mods`` when unset.

        :rtype: Mods[ModificationTags]
        """
        if self._labile_mods is None:
            return EMPTY_LABILE_MODS
        return Mods[ModificationTags](mod_type=ModType.LABILE, _mods=self._labile_mods)

    @labile_mods.setter
    def labile_mods(self, value: Any) -> None:
        self.set_labile_mods(value, inplace=True, validate=self._validate)

    @property
    def labile_mods_str(self) -> str:
        """Serialised labile modifications string, or empty string when unset.

        :rtype: str
        """
        if self._labile_mods is None:
            return ""
        return self.labile_mods.serialize()

    @property
    def unknown_mods(self) -> Mods[ModificationTags]:
        """Unknown-localisation modifications; returns an empty ``Mods`` when unset.

        :rtype: Mods[ModificationTags]
        """
        if self._unknown_mods is None:
            return EMPTY_UNKNOWN_MODS
        return Mods[ModificationTags](mod_type=ModType.UNKNOWN, _mods=self._unknown_mods)

    @unknown_mods.setter
    def unknown_mods(self, value: Any) -> None:
        self.set_unknown_mods(value, inplace=True, validate=self._validate)

    @property
    def unknown_mods_str(self) -> str:
        """Serialised unknown modifications string, or empty string when unset.

        :rtype: str
        """
        if self._unknown_mods is None:
            return ""
        return self.unknown_mods.serialize()

    @property
    def nterm_mods(self) -> Mods[ModificationTags]:
        """N-terminal modifications; returns an empty ``Mods`` when unset.

        :rtype: Mods[ModificationTags]
        """
        if self._nterm_mods is None:
            return EMPTY_NTERM_MODS
        return Mods[ModificationTags](mod_type=ModType.NTERM, _mods=self._nterm_mods)

    @nterm_mods.setter
    def nterm_mods(self, value: Any) -> None:
        self.set_nterm_mods(value, inplace=True, validate=self._validate)

    @property
    def nterm_mods_str(self) -> str:
        """Serialised N-terminal modifications string, or empty string when unset.

        :rtype: str
        """
        if self._nterm_mods is None:
            return ""
        return self.nterm_mods.serialize()

    @property
    def cterm_mods(self) -> Mods[ModificationTags]:
        """C-terminal modifications; returns an empty ``Mods`` when unset.

        :rtype: Mods[ModificationTags]
        """
        if self._cterm_mods is None:
            return EMPTY_CTERM_MODS
        return Mods[ModificationTags](mod_type=ModType.CTERM, _mods=self._cterm_mods)

    @cterm_mods.setter
    def cterm_mods(self, value: Any) -> None:
        self.set_cterm_mods(value, inplace=True, validate=self._validate)

    @property
    def cterm_mods_str(self) -> str:
        """Serialised C-terminal modifications string, or empty string when unset.

        :rtype: str
        """
        if self._cterm_mods is None:
            return ""
        return self.cterm_mods.serialize()

    @property
    def internal_mods(self) -> dict[int, Mods[ModificationTags]]:
        """Per-position internal modifications keyed by 0-based index; empty dict when unset.

        :rtype: dict[int, Mods[ModificationTags]]
        """
        if self._internal_mods is None:
            return {}
        internal_mods_parsed: dict[int, Mods[ModificationTags]] = {}
        for pos, mods_dict in self._internal_mods.items():
            internal_mods_parsed[pos] = Mods[ModificationTags](mod_type=ModType.INTERNAL, _mods=mods_dict)
        return internal_mods_parsed

    @internal_mods.setter
    def internal_mods(self, value: dict[int, Any] | None) -> None:
        self.set_internal_mods(value, inplace=True, validate=self._validate)

    @property
    def validate(self) -> bool:
        return self._validate

    @validate.setter
    def validate(self, value: bool) -> None:
        self._validate = value
        if self.has_intervals:
            for interval in self.intervals:
                interval._validate = value

    def has_internal_mods_at_index(self, position: int) -> bool:
        """Check if there are any modifications at a specific position in the sequence."""
        if self._internal_mods is None:
            return False

        mods_dict = self._internal_mods.get(position, None)
        if mods_dict is None or len(mods_dict) == 0:
            return False

        return True

    def get_internal_mod_indexes(self) -> list[int]:
        """Get a list of all indexes that have internal modifications."""
        if self._internal_mods is None:
            return []
        return list(self._internal_mods.keys())

    def get_internal_mods_str_at_index(self, position: int) -> str:
        """Get the modification string at a specific position in the sequence."""
        if self.has_internal_mods_at_index(position) is False:
            return ""

        return self.get_internal_mods_at_index(position).serialize()

    def get_internal_mods_at_index(self, position: int) -> Mods[ModificationTags]:
        """Get all modifications at a specific position in the sequence."""
        if self.has_internal_mods_at_index(position) is False:
            return EMPTY_INTERNAL_MODS
        return Mods[ModificationTags](
            mod_type=ModType.INTERNAL,
            _mods=self._internal_mods[position],  # type: ignore
        )

    @property
    def intervals(self) -> tuple[Interval, ...]:
        """Ambiguous sequence intervals; empty tuple when unset.

        :rtype: tuple[Interval, ...]
        """
        return tuple(self._intervals) if self._intervals is not None else ()

    @intervals.setter
    def intervals(self, value: list[Interval] | None) -> None:
        self.set_intervals(value, inplace=True, validate=self._validate)

    @property
    def charge(self) -> int | Mods[GlobalChargeCarrier] | None:
        """Charge as an integer or adduct ``Mods``; ``None`` when uncharged or zero.

        :rtype: int | Mods[GlobalChargeCarrier] | None
        """
        if isinstance(self._charge, int):
            if self._charge == 0:
                return None
            return self._charge
        elif isinstance(self._charge, list):
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={mod_str: 1 for mod_str in self._charge})

        return None

    @charge.setter
    def charge(self, value: int | str | list[str] | Mods[GlobalChargeCarrier] | None) -> None:
        self.set_charge(value, inplace=True, validate=self._validate)

    @property
    def charge_state(self) -> int:
        """Numeric charge state derived from the stored charge; 0 when uncharged.

        :rtype: int
        :raises ValueError: If the stored charge value has an unexpected type.
        """
        charge = self.charge
        if isinstance(charge, int):
            return charge
        elif isinstance(charge, Mods):
            return sum(mod.get_charge() for mod in charge.mods)
        elif charge is None:
            return 0
        else:
            raise ValueError(f"Invalid charge type: {type(charge)}")

    @property
    def charge_adducts(self) -> Mods[GlobalChargeCarrier]:
        """Charge expressed as adduct ``Mods``; converts an integer charge to proton adducts.

        :rtype: Mods[GlobalChargeCarrier]
        """
        charge = self.charge
        if isinstance(charge, int):
            if charge == 0 or charge < 0:
                return EMPTY_CHARGE_MODS
            elif charge > 0:
                s = str(
                    GlobalChargeCarrier(
                        charged_formula=H_CHARGE_FORMULA,
                        occurance=charge,
                    )
                )
                return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={s: 1})
        elif isinstance(charge, Mods):
            return charge
        return EMPTY_CHARGE_MODS

    """
    Set Methods - Replace existing modifications
    """

    def set_sequence(self, sequence: str | None, inplace: bool = True, validate: bool | None = None) -> Self:
        """Set the amino-acid sequence.

        :param sequence: New sequence, or ``None`` to clear.
        :type sequence: str | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().set_sequence(sequence, inplace=True, validate=validate)
        self._sequence = sequence
        if validate:
            self.validate_sequence()
        return self

    def set_compound_name(self, name: str | None, inplace: bool = True, validate: bool | None = None) -> Self:
        """Set the compound-level name.

        :param name: New name, or ``None`` to clear.
        :type name: str | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._set_name_generic(name, "_compound_name", inplace, validate)

    def set_ion_name(self, name: str | None, inplace: bool = True, validate: bool | None = None) -> Self:
        """Set the ion-level name.

        :param name: New name, or ``None`` to clear.
        :type name: str | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._set_name_generic(name, "_ion_name", inplace, validate)

    def set_peptide_name(self, name: str | None, inplace: bool = True, validate: bool | None = None) -> Self:
        """Set the peptide-level name.

        :param name: New name, or ``None`` to clear.
        :type name: str | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._set_name_generic(name, "_peptide_name", inplace, validate)

    def set_isotope_mods(
        self,
        mods: dict[str, int] | Mods[IsotopeReplacement] | None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace all isotope modifications.

        :param mods: New isotope modifications, or ``None`` to clear.
        :type mods: dict[str, int] | Mods[IsotopeReplacement] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._set_mod_generic(mods, "_isotope_mods", "validate_isotope_mods", inplace, validate)

    def set_static_mods(
        self,
        mods: dict[str, int] | Mods[FixedModification] | None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace all static (fixed) modifications.

        :param mods: New static modifications, or ``None`` to clear.
        :type mods: dict[str, int] | Mods[FixedModification] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._set_mod_generic(mods, "_static_mods", "validate_static_mods", inplace, validate)

    def set_labile_mods(
        self,
        mods: dict[str, int] | Mods[ModificationTags] | None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace all labile modifications.

        :param mods: New labile modifications, or ``None`` to clear.
        :type mods: dict[str, int] | Mods[ModificationTags] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._set_mod_generic(mods, "_labile_mods", "validate_labile_mods", inplace, validate)

    def set_unknown_mods(
        self,
        mods: dict[str, int] | Mods[ModificationTags] | None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace all unknown-localisation modifications.

        :param mods: New unknown modifications, or ``None`` to clear.
        :type mods: dict[str, int] | Mods[ModificationTags] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._set_mod_generic(mods, "_unknown_mods", "validate_unknown_mods", inplace, validate)

    def set_nterm_mods(
        self,
        mods: dict[str, int] | Mods[ModificationTags] | None,
        inplace: bool = True,
        validate: bool | None = None,
        start_aa: str | None = None,
    ) -> Self:
        """Replace all N-terminal modifications.

        :param mods: New N-terminal modifications, or ``None`` to clear.
        :type mods: dict[str, int] | Mods[ModificationTags] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :param start_aa: Only apply if the sequence starts with this residue; no-op otherwise.
        :type start_aa: str | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if start_aa is not None:
            if self.start_aa != start_aa:
                return self if inplace else self.copy()
        return self._set_mod_generic(mods, "_nterm_mods", "validate_nterm_mods", inplace, validate)

    def set_cterm_mods(
        self,
        mods: dict[str, int] | Mods[ModificationTags] | None,
        inplace: bool = True,
        validate: bool | None = None,
        end_aa: str | None = None,
    ) -> Self:
        """Replace all C-terminal modifications.

        :param mods: New C-terminal modifications, or ``None`` to clear.
        :type mods: dict[str, int] | Mods[ModificationTags] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :param end_aa: Only apply if the sequence ends with this residue; no-op otherwise.
        :type end_aa: str | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if end_aa is not None:
            if self.end_aa != end_aa:
                return self if inplace else self.copy()
        return self._set_mod_generic(mods, "_cterm_mods", "validate_cterm_mods", inplace, validate)

    def set_internal_mods(
        self,
        mods: dict[int, dict[str, int] | Mods[ModificationTags] | None] | None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace all internal (per-position) modifications.

        :param mods: Mapping of 0-based residue index to modification dict, or ``None`` to clear.
        :type mods: dict[int, dict[str, int] | Mods[ModificationTags] | None] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate

        if not inplace:
            return self.copy().set_internal_mods(mods, inplace=True, validate=validate)

        if mods is None:
            self._internal_mods = None
            return self

        internal_mods: dict[int, dict[str, int]] = {}
        for pos, mods_dict in mods.items():
            internalmod = convert_moddict_input(mods_dict)
            if internalmod is None or len(internalmod) == 0:
                continue
            internal_mods[pos] = internalmod

        if len(internal_mods) == 0:
            self._internal_mods = None
            return self

        self._internal_mods = internal_mods
        if validate:
            self.validate_internal_mods()
        return self

    def set_intervals(
        self,
        intervals: list[Interval] | None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace all ambiguous sequence intervals.

        :param intervals: New list of intervals, or ``None`` to clear.
        :type intervals: list[Interval] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate

        if not inplace:
            return self.copy().set_intervals(intervals, inplace=True, validate=validate)

        if intervals is None:
            self._intervals = None
            return self

        if len(intervals) == 0:
            self._intervals = None
            return self

        self._intervals = intervals.copy()
        if validate:
            self.validate_intervals()
        return self

    def set_internal_mods_at_index(self, index: int, mods: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Replace internal modifications at a single 0-based sequence position.

        :param index: 0-based residue index.
        :type index: int
        :param mods: New modifications for this position, or ``None`` to remove them.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().set_internal_mods_at_index(index, mods, inplace=True, validate=validate)

        if mods is None:
            # remove mod at index
            if self._internal_mods is not None and index in self._internal_mods:
                del self._internal_mods[index]
            return self

        mods = convert_moddict_input(mods)

        if len(mods) == 0:
            # remove mod at index
            if self._internal_mods is not None and index in self._internal_mods:
                del self._internal_mods[index]
            return self

        if validate:
            if not Mods[ModificationTags](mod_type=ModType.INTERNAL, _mods=mods).is_valid:
                raise ValueError(f"Invalid internal modifications at position {index}")

        if self._internal_mods is None:
            self._internal_mods = {}

        self._internal_mods[index] = mods
        return self

    def set_charge(
        self,
        charge: int | str | list[str] | Mods[GlobalChargeCarrier] | GlobalChargeCarrier | Mod[GlobalChargeCarrier] | None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace the charge value.

        :param charge: New charge as an integer, adduct string(s), ``Mods``, or ``None`` to clear.
        :type charge: int | str | list[str] | Mods[GlobalChargeCarrier] | GlobalChargeCarrier | Mod[GlobalChargeCarrier] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        :raises ValueError: If the resolved charge value has an unsupported type.
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().set_charge(charge, inplace=True, validate=validate)

        set_value: None | int | list[str] = None
        if isinstance(charge, int):
            if charge == 0:
                set_value = None
            set_value = charge
        elif isinstance(charge, str):
            set_value = [charge]
        elif isinstance(charge, list):
            if len(charge) == 0:
                set_value = None
            else:
                set_value = [str(c) for c in charge]
                if len(set_value) == 0:
                    set_value = None
        elif charge is None:
            set_value = None
        elif isinstance(charge, Mods):
            set_value = [str(c) for c in charge._mods] if charge._mods is not None else None
        elif isinstance(charge, Mod):
            set_value = [str(charge)]
        elif isinstance(charge, GlobalChargeCarrier):
            set_value = [str(charge)]

        if not isinstance(set_value, (int, list)) and set_value is not None:
            raise ValueError(f"Invalid charge value: {set_value}")

        self._charge: None | int | list[str] = set_value

        if validate:
            self.validate_charge()

        return self

    def _set_name_generic(
        self,
        name: Any | None,
        attr_name: str,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy()._set_name_generic(name, attr_name, inplace=True, validate=validate)

        if name is not None and not isinstance(name, str):
            name = str(name)
        if name == "":
            name = None
        setattr(self, attr_name, name)
        return self

    def _set_mod_generic(
        self,
        mods: Any,
        attr_name: str,
        validator_method_name: str | None = None,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy()._set_mod_generic(mods, attr_name, validator_method_name, inplace=True, validate=validate)

        if mods is None:
            setattr(self, attr_name, None)
            return self

        converted_mods = convert_moddict_input(mods)
        if len(converted_mods) == 0:
            setattr(self, attr_name, None)
            return self

        setattr(self, attr_name, converted_mods)

        if validate and validator_method_name:
            getattr(self, validator_method_name)()

        return self

    def _set_mod_by_type(
        self,
        value: Any,
        mod_type: ModType,
    ) -> Self:
        match mod_type:
            case ModType.ISOTOPE:
                self.isotope_mods = value
            case ModType.STATIC:
                self.static_mods = value
            case ModType.LABILE:
                self.labile_mods = value
            case ModType.UNKNOWN:
                self.unknown_mods = value
            case ModType.NTERM:
                self.nterm_mods = value
            case ModType.CTERM:
                self.cterm_mods = value
            case ModType.INTERNAL:
                self.internal_mods = value
            case ModType.INTERVAL:
                self.intervals = value
            case ModType.CHARGE:
                self.charge = value
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

        return self

    def set_mods(
        self,
        mods: Mapping[ModType | ModTypeLiteral | int, Any] | None,
        inplace: bool = True,
    ) -> Self:
        """Set a modification by type, replacing any existing mods of that type"""

        if not inplace:
            return self.copy().set_mods(mods=mods, inplace=True)

        if mods is None:
            self.clear_mods(inplace=True)
            return self

        for mod_type, mod_value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise IndexError(f"Internal modification index out of range: {mod_type}")
                self.set_internal_mods_at_index(mod_type, mod_value, inplace=True)
                continue

            self._set_mod_by_type(mod_value, ModType(mod_type))

        return self

    """
    Append Methods
    """

    def _append_mod_generic(
        self,
        mod: Any,
        attr_name: str,
        validator: Callable[[str], Any],
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        if validate is None:
            validate = self._validate

        if not inplace:
            return self.copy()._append_mod_generic(mod, attr_name, validator, inplace=True, validate=validate)

        mod_str, count = convert_single_mod_input(mod)

        if validate:
            if not validator(mod_str).is_valid:
                raise ValueError(f"Invalid modification: {mod_str}")

        mod_dict = getattr(self, attr_name)
        if mod_dict is None:
            setattr(self, attr_name, {})
            mod_dict = getattr(self, attr_name)

        if mod_str in mod_dict:
            mod_dict[mod_str] += count
        else:
            mod_dict[mod_str] = count

        return self

    def append_isotope_mod(self, mod: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Append an isotope modification.

        :param mod: Modification to append.
        :type mod: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._append_mod_generic(mod, "_isotope_mods", IsotopeReplacement.from_string, inplace, validate)

    def append_static_mod(self, mod: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Append a static (fixed) modification.

        :param mod: Modification to append.
        :type mod: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._append_mod_generic(mod, "_static_mods", FixedModification.from_string, inplace, validate)

    def append_labile_mod(self, mod: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Append a labile modification.

        :param mod: Modification to append.
        :type mod: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._append_mod_generic(mod, "_labile_mods", ModificationTags.from_string, inplace, validate)

    def append_unknown_mod(self, mod: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Append an unknown-localisation modification.

        :param mod: Modification to append.
        :type mod: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._append_mod_generic(mod, "_unknown_mods", ModificationTags.from_string, inplace, validate)

    def append_nterm_mod(
        self,
        mod: Any,
        inplace: bool = True,
        validate: bool | None = None,
        start_aa: str | None = None,
    ) -> Self:
        """Append an N-terminal modification.

        :param mod: Modification to append.
        :type mod: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :param start_aa: Only apply if the sequence starts with this residue; no-op otherwise.
        :type start_aa: str | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if start_aa is not None:
            if self.start_aa != start_aa:
                return self if inplace else self.copy()
        return self._append_mod_generic(mod, "_nterm_mods", ModificationTags.from_string, inplace, validate)

    def append_cterm_mod(
        self,
        mod: Any,
        inplace: bool = True,
        validate: bool | None = None,
        end_aa: str | None = None,
    ) -> Self:
        """Append a C-terminal modification.

        :param mod: Modification to append.
        :type mod: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :param end_aa: Only apply if the sequence ends with this residue; no-op otherwise.
        :type end_aa: str | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if end_aa is not None:
            if self.end_aa != end_aa:
                return self if inplace else self.copy()
        return self._append_mod_generic(mod, "_cterm_mods", ModificationTags.from_string, inplace, validate)

    def append_internal_mod_at_index(self, index: int, mod: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Append an internal modification at a specific 0-based sequence position.

        :param index: 0-based residue index.
        :type index: int
        :param mod: Modification to append.
        :type mod: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate

        if not inplace:
            return self.copy().append_internal_mod_at_index(index, mod, inplace=True, validate=validate)

        mod_str, count = convert_single_mod_input(mod)

        if validate:
            if not ModificationTags.from_string(mod_str).is_valid:
                raise ValueError(f"Invalid modification: {mod_str}")

        if self._internal_mods is None:
            self._internal_mods = {}

        if index not in self._internal_mods:
            self._internal_mods[index] = {}

        if mod_str in self._internal_mods[index]:
            self._internal_mods[index][mod_str] += count
        else:
            self._internal_mods[index][mod_str] = count

        return self

    def append_interval(
        self,
        interval: Interval | tuple[int, int, bool, Any],
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Append an ambiguous sequence interval.

        :param interval: ``Interval`` object or a ``(start, end, ambiguous, mods)`` tuple.
        :type interval: Interval | tuple[int, int, bool, Any]
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().append_interval(interval, inplace=True, validate=validate)

        if isinstance(interval, tuple):
            start, end, ambiguous, mods_input = interval
            mods_converted = convert_moddict_input(mods_input)
            interval = Interval(
                start=start,
                end=end,
                ambiguous=ambiguous,
                mods=mods_converted,
                validate=validate,
            )
        else:
            interval = interval.copy()
            interval._validate = validate

        if validate:
            if not isinstance(interval, Interval):
                raise TypeError(f"Expected Interval object, got {type(interval)}")
            if not interval.is_valid:
                raise ValueError(f"Invalid interval: {interval}")

        if self._intervals is None:
            self._intervals = []

        self._intervals.append(interval)
        return self

    def _append_by_type(
        self,
        value: Any,
        mod_type: ModType,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        if not inplace:
            return self.copy()._append_by_type(value, mod_type, inplace=True, validate=validate)

        match mod_type:
            case ModType.ISOTOPE:
                self.append_isotope_mod(value, inplace=True, validate=validate)
            case ModType.STATIC:
                self.append_static_mod(value, inplace=True, validate=validate)
            case ModType.LABILE:
                self.append_labile_mod(value, inplace=True, validate=validate)
            case ModType.UNKNOWN:
                self.append_unknown_mod(value, inplace=True, validate=validate)
            case ModType.NTERM:
                self.append_nterm_mod(value, inplace=True, validate=validate)
            case ModType.CTERM:
                self.append_cterm_mod(value, inplace=True, validate=validate)
            case ModType.INTERNAL:
                for key, val in value.items():
                    self.append_internal_mod_at_index(key, val, inplace=True, validate=validate)
            case ModType.INTERVAL:
                self.append_interval(value, inplace=True, validate=validate)
            case ModType.CHARGE:
                self.set_charge(value, inplace=True, validate=validate)
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

        return self

    def append_mods(
        self,
        mods: Mapping[ModType | ModTypeLiteral | int, Any],
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Append modifications of multiple types from a mapping of mod-type to value.

        :param mods: Mapping of :class:`ModType` (or literal/index) to modification value.
        :type mods: Mapping[ModType | ModTypeLiteral | int, Any]
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        :raises IndexError: If an integer key is out of range for the current sequence.
        """
        if not inplace:
            return self.copy().append_mods(mods, inplace=True, validate=validate)

        for mod_type, value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise IndexError(f"Internal modification index out of range: {mod_type}")
                self.append_internal_mod_at_index(mod_type, value, inplace=True, validate=validate)
                continue

            self._append_by_type(value, ModType(mod_type), inplace=True, validate=validate)

        return self

    """
    Extend Methods - Add multiple modifications
    """

    def _extend_generic(
        self,
        mods: Any,
        append_method: Callable[[Any, bool, bool | None], Self],
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy()._extend_generic(mods, append_method, inplace=True, validate=validate)
        if mods is not None:
            for mod in mods:
                append_method(mod, inplace=True, validate=validate)  # type: ignore
        return self

    def extend_isotope_mods(self, mods: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Extend isotope modifications by appending each item in *mods*.

        :param mods: Iterable of modifications to append.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._extend_generic(mods, self.append_isotope_mod, inplace, validate)

    def extend_static_mods(self, mods: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Extend static modifications by appending each item in *mods*.

        :param mods: Iterable of modifications to append.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._extend_generic(mods, self.append_static_mod, inplace, validate)

    def extend_labile_mods(self, mods: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Extend labile modifications by appending each item in *mods*.

        :param mods: Iterable of modifications to append.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._extend_generic(mods, self.append_labile_mod, inplace, validate)

    def extend_unknown_mods(self, mods: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Extend unknown-localisation modifications by appending each item in *mods*.

        :param mods: Iterable of modifications to append.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._extend_generic(mods, self.append_unknown_mod, inplace, validate)

    def extend_nterm_mods(
        self,
        mods: Any,
        inplace: bool = True,
        validate: bool | None = None,
        start_aa: str | None = None,
    ) -> Self:
        """Extend N-terminal modifications by appending each item in *mods*.

        :param mods: Iterable of modifications to append.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :param start_aa: Only apply if the sequence starts with this residue; no-op otherwise.
        :type start_aa: str | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().extend_nterm_mods(mods, inplace=True, validate=validate, start_aa=start_aa)
        if start_aa is not None:
            if self.start_aa != start_aa:
                return self
        if mods is not None:
            for mod in mods:
                self.append_nterm_mod(mod, inplace=True, validate=validate, start_aa=start_aa)
        return self

    def extend_cterm_mods(
        self,
        mods: Any,
        inplace: bool = True,
        validate: bool | None = None,
        end_aa: str | None = None,
    ) -> Self:
        """Extend C-terminal modifications by appending each item in *mods*.

        :param mods: Iterable of modifications to append.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :param end_aa: Only apply if the sequence ends with this residue; no-op otherwise.
        :type end_aa: str | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().extend_cterm_mods(mods, inplace=True, validate=validate, end_aa=end_aa)
        if end_aa is not None:
            if self.end_aa != end_aa:
                return self
        if mods is not None:
            for mod in mods:
                self.append_cterm_mod(mod, inplace=True, validate=validate, end_aa=end_aa)
        return self

    def extend_internal_mods_at_index(self, index: int, mods: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Extend internal modifications at a single position by appending each item in *mods*.

        :param index: 0-based residue index.
        :type index: int
        :param mods: Iterable of modifications to append at this position.
        :type mods: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().extend_internal_mods_at_index(index, mods, inplace=True, validate=validate)
        if mods is not None:
            for mod in mods:
                self.append_internal_mod_at_index(index, mod, inplace=True, validate=validate)
        return self

    def extend_intervals(self, intervals: Any, inplace: bool = True, validate: bool | None = None) -> Self:
        """Extend ambiguous sequence intervals by appending each item in *intervals*.

        :param intervals: Iterable of intervals to append.
        :type intervals: Any
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._extend_generic(intervals, self.append_interval, inplace, validate)

    def _extend_by_type(
        self,
        value: Any,
        mod_type: ModType,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy()._extend_by_type(value, mod_type, inplace=True, validate=validate)

        match mod_type:
            case ModType.ISOTOPE:
                self.extend_isotope_mods(value, inplace=True, validate=validate)
            case ModType.STATIC:
                self.extend_static_mods(value, inplace=True, validate=validate)
            case ModType.LABILE:
                self.extend_labile_mods(value, inplace=True, validate=validate)
            case ModType.UNKNOWN:
                self.extend_unknown_mods(value, inplace=True, validate=validate)
            case ModType.NTERM:
                self.extend_nterm_mods(value, inplace=True, validate=validate)
            case ModType.CTERM:
                self.extend_cterm_mods(value, inplace=True, validate=validate)
            case ModType.INTERNAL:
                for index, mod in value.items():
                    self.extend_internal_mods_at_index(index, mod, inplace=True, validate=validate)
            case ModType.INTERVAL:
                self.extend_intervals(value, inplace=True, validate=validate)
            case ModType.CHARGE:
                raise NotImplementedError("Extending charge not supported.")
            case _:
                raise NotImplementedError(f"Appending {mod_type} not supported.")

        return self

    def extend_mods(
        self,
        mods: Mapping[ModType | ModTypeLiteral | int, Any],
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Extend modifications of multiple types by iterating through each mapped iterable.

        :param mods: Mapping of :class:`ModType` (or literal/index) to iterable of modification values.
        :type mods: Mapping[ModType | ModTypeLiteral | int, Any]
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        :raises IndexError: If an integer key is out of range for the current sequence.
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().extend_mods(mods, inplace=True, validate=validate)

        for mod_type, value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise IndexError(f"Internal modification index out of range: {mod_type}")
                self.extend_internal_mods_at_index(mod_type, value, inplace=True, validate=validate)
                continue
            self._extend_by_type(value, ModType(mod_type), inplace=True, validate=validate)

        return self

    """
    REMOVE Methods
    """

    def remove_mods(self, mods: Mapping[ModType | ModTypeLiteral | int, Any], inplace: bool = True) -> Self:
        """Remove modifications by decrementing their counts."""
        if not inplace:
            return self.copy().remove_mods(mods, inplace=True)

        for mod_type, mod_value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise IndexError(f"Internal modification index out of range: {mod_type}")
                self.remove_internal_mod_at_index(mod_type, mod_value, inplace=True)
                continue

            match ModType(mod_type):
                case ModType.ISOTOPE:
                    self.remove_isotope_mod(mod_value, inplace=True)
                case ModType.STATIC:
                    self.remove_static_mod(mod_value, inplace=True)
                case ModType.LABILE:
                    self.remove_labile_mod(mod_value, inplace=True)
                case ModType.UNKNOWN:
                    self.remove_unknown_mod(mod_value, inplace=True)
                case ModType.NTERM:
                    self.remove_nterm_mod(mod_value, inplace=True)
                case ModType.CTERM:
                    self.remove_cterm_mod(mod_value, inplace=True)
                case ModType.INTERNAL:
                    raise NotImplementedError("Use remove_internal_mod_at_index for internal modifications.")
                case ModType.INTERVAL:
                    self.remove_interval(mod_value, inplace=True)
                case ModType.CHARGE:
                    raise NotImplementedError("Removing charge modifications not supported.")
                case _:
                    raise TypeError(f"Unknown mod type: {mod_type}")

        return self

    def _remove_mod_generic(
        self,
        mod: Any,
        attr_name: str,
        inplace: bool = True,
    ) -> Self:
        """Generic method to remove a modification by decrementing its count."""
        if not inplace:
            return self.copy()._remove_mod_generic(mod, attr_name, inplace=True)

        mod_dict = getattr(self, attr_name)
        if mod_dict is None:
            return self

        mod_str, count = convert_single_mod_input(mod)

        if mod_str not in mod_dict:
            return self

        # Decrement count, ensuring it doesn't go below 0
        mod_dict[mod_str] = max(0, mod_dict[mod_str] - count)

        # Remove if count reaches 0
        if mod_dict[mod_str] == 0:
            del mod_dict[mod_str]

        # Clean up if dict is now empty
        if len(mod_dict) == 0:
            setattr(self, attr_name, None)

        return self

    def remove_isotope_mod(self, mod: Any, inplace: bool = True) -> Self:
        """Remove a specific isotope modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_isotope_mods", inplace)

    def remove_static_mod(self, mod: Any, inplace: bool = True) -> Self:
        """Remove a specific static modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_static_mods", inplace)

    def remove_labile_mod(self, mod: Any, inplace: bool = True) -> Self:
        """Remove a specific labile modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_labile_mods", inplace)

    def remove_unknown_mod(self, mod: Any, inplace: bool = True) -> Self:
        """Remove a specific unknown modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_unknown_mods", inplace)

    def remove_nterm_mod(self, mod: Any, inplace: bool = True, start_aa: str | None = None) -> Self:
        """Remove a specific N-terminal modification by decrementing its count."""
        if start_aa is not None and self.start_aa != start_aa:
            return self if inplace else self.copy()
        return self._remove_mod_generic(mod, "_nterm_mods", inplace)

    def remove_cterm_mod(self, mod: Any, inplace: bool = True, end_aa: str | None = None) -> Self:
        """Remove a specific C-terminal modification by decrementing its count."""
        if end_aa is not None and self.end_aa != end_aa:
            return self if inplace else self.copy()
        return self._remove_mod_generic(mod, "_cterm_mods", inplace)

    def remove_internal_mod_at_index(self, index: int, mod: Any, inplace: bool = True) -> Self:
        """Remove a specific internal modification at a position by decrementing its count."""
        if not inplace:
            return self.copy().remove_internal_mod_at_index(index, mod, inplace=True)

        if self._internal_mods is None or index not in self._internal_mods:
            return self

        mod_str, count = convert_single_mod_input(mod)

        if mod_str not in self._internal_mods[index]:
            return self

        # Decrement count, ensuring it doesn't go below 0
        self._internal_mods[index][mod_str] = max(0, self._internal_mods[index][mod_str] - count)

        # Remove if count reaches 0
        if self._internal_mods[index][mod_str] == 0:
            del self._internal_mods[index][mod_str]

        # Remove position if no mods left
        if len(self._internal_mods[index]) == 0:
            del self._internal_mods[index]

        # Clean up if internal_mods is now empty
        if len(self._internal_mods) == 0:
            self._internal_mods = None

        return self

    def remove_interval(self, interval: Interval, inplace: bool = True) -> Self:
        """Remove a specific interval from the intervals list."""
        if not inplace:
            return self.copy().remove_interval(interval, inplace=True)

        if self._intervals is None:
            return self

        try:
            self._intervals.remove(interval)
        except ValueError:
            # Interval not found, just return
            pass

        if len(self._intervals) == 0:
            self._intervals = None

        return self

    """
    Magic Methods
    """

    def compare(self, other: Self) -> bool:
        """Compare this annotation to *other* field-by-field, logging differences at DEBUG level.

        :param other: The annotation to compare against.
        :type other: Self
        :return: ``True`` if all fields are equal, ``False`` otherwise.
        :rtype: bool
        """
        # check each attribute for equality
        diffs = []  # hold the string values of differing attributes (ACTUALLY SHOW THE ATTRIBUTES)
        for attr in [
            "_sequence",
            "_isotope_mods",
            "_static_mods",
            "_labile_mods",
            "_unknown_mods",
            "_nterm_mods",
            "_cterm_mods",
            "_internal_mods",
            "_intervals",
            "_charge",
        ]:
            if getattr(self, attr) != getattr(other, attr):
                self_attribute = getattr(self, attr)
                other_attribute = getattr(other, attr)
                diffs.append(f"{attr} (self: {self_attribute}, other: {other_attribute})")
        if diffs:
            logger.debug("Differences found in attributes: %s", ", ".join(diffs))
            return False
        return True

    def __len__(self) -> int:
        return len(self.stripped_sequence)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ProFormaAnnotation):
            raise NotImplementedError(f"Cannot compare ProFormaAnnotationBase with {type(other)}")

        return (
            self.sequence == other.sequence
            and self._isotope_mods == other._isotope_mods
            and self._static_mods == other._static_mods
            and self._labile_mods == other._labile_mods
            and self._unknown_mods == other._unknown_mods
            and self._nterm_mods == other._nterm_mods
            and self._cterm_mods == other._cterm_mods
            and self._internal_mods == other._internal_mods
            and self._intervals == other._intervals
            and self._charge == other._charge
        )

    def __repr__(self) -> str:
        seq = f"ProFormaAnnot(sequence={self.sequence}"

        if self.has_isotope_mods:
            seq += f", {ModType.ISOTOPE.value}={self.isotope_mods}"
        if self.has_static_mods:
            seq += f", {ModType.STATIC.value}={self.static_mods}"
        if self.has_labile_mods:
            seq += f", {ModType.LABILE.value}={self.labile_mods}"
        if self.has_unknown_mods:
            seq += f", {ModType.UNKNOWN.value}={self.unknown_mods}"
        if self.has_nterm_mods:
            seq += f", {ModType.NTERM.value}={self.nterm_mods}"
        if self.has_cterm_mods:
            seq += f", {ModType.CTERM.value}={self.cterm_mods}"
        if self.has_internal_mods:
            internal_mod_items = ", ".join(f"{pos}: {mods}" for pos, mods in sorted(self.internal_mods.items()))
            seq += f", {ModType.INTERNAL.value}={{{internal_mod_items}}}"
        if self.has_intervals:
            seq += f", {ModType.INTERVAL.value}={self.intervals}"
        if self.has_charge:
            seq += f", {ModType.CHARGE.value}={self.charge}"
        seq += ")"

        return seq

    def __str__(self) -> str:
        return self.serialize()

    def __hash__(self):
        return hash(
            (
                self._sequence,
                tuple(self._isotope_mods.items()) if self._isotope_mods else None,
                tuple(self._static_mods.items()) if self._static_mods else None,
                tuple(self._labile_mods.items()) if self._labile_mods else None,
                tuple(self._unknown_mods.items()) if self._unknown_mods else None,
                tuple(self._nterm_mods.items()) if self._nterm_mods else None,
                tuple(self._cterm_mods.items()) if self._cterm_mods else None,
                tuple((pos, tuple(sorted(mods.items()))) for pos, mods in sorted(self._internal_mods.items())) if self._internal_mods else None,
                tuple(self._intervals) if self._intervals else None,
                self._charge,
            )
        )

    def copy(self) -> Self:
        """Return a deep copy of this annotation.

        :return: A new independent annotation with the same field values.
        :rtype: Self
        """
        return self.__class__(
            sequence=self._sequence,
            compound_name=self._compound_name,
            ion_name=self._ion_name,
            peptide_name=self._peptide_name,
            isotope_mods=self._isotope_mods.copy() if self._isotope_mods is not None else None,
            static_mods=self._static_mods.copy() if self._static_mods is not None else None,
            labile_mods=self._labile_mods.copy() if self._labile_mods is not None else None,
            unknown_mods=self._unknown_mods.copy() if self._unknown_mods is not None else None,
            nterm_mods=self._nterm_mods.copy() if self._nterm_mods is not None else None,
            cterm_mods=self._cterm_mods.copy() if self._cterm_mods is not None else None,
            internal_mods={pos: mods.copy() for pos, mods in self._internal_mods.items()} if self._internal_mods is not None else None,
            intervals=self._intervals.copy() if self._intervals is not None else None,
            charge=self._charge,
            validate=self._validate,
        )

    def update(self, other: Self) -> None:
        """In-place update: replace all fields of this annotation with copies from *other*.

        :param other: Source annotation whose field values will be copied.
        :type other: Self
        """
        self._sequence = other._sequence
        self._compound_name = other._compound_name
        self._ion_name = other._ion_name
        self._peptide_name = other._peptide_name
        self._isotope_mods = other._isotope_mods.copy() if other._isotope_mods is not None else None
        self._static_mods = other._static_mods.copy() if other._static_mods is not None else None
        self._labile_mods = other._labile_mods.copy() if other._labile_mods is not None else None
        self._unknown_mods = other._unknown_mods.copy() if other._unknown_mods is not None else None
        self._nterm_mods = other._nterm_mods.copy() if other._nterm_mods is not None else None
        self._cterm_mods = other._cterm_mods.copy() if other._cterm_mods is not None else None
        self._internal_mods = {pos: mods.copy() for pos, mods in other._internal_mods.items()} if other._internal_mods is not None else None
        self._intervals = other._intervals.copy() if other._intervals is not None else None
        self._charge = other._charge

    def __getitem__(self, key: int | slice | Span | tuple[int, int, int]) -> Self:
        if isinstance(key, (tuple, Span)):
            return self.slice_by_span(key, inplace=False)
        if isinstance(key, slice):
            start, stop, step = key.start, key.stop, key.step
            if step is not None and step != 1:
                raise ValueError("Step slicing not supported")
            return self.slice(start, stop, inplace=False)
        elif isinstance(key, int):
            raise NotImplementedError("Single index access not supported for ProFormaAnnotation")

    def sort_mods(self, inplace: bool = True) -> Self:
        """Sort all modification dictionaries and the intervals list deterministically.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation with sorted modifications.
        :rtype: Self
        """
        if not inplace:
            return self.copy().sort_mods(inplace=True)

        if self._isotope_mods is not None:
            self._isotope_mods = dict(sorted(self._isotope_mods.items()))
        if self._static_mods is not None:
            self._static_mods = dict(sorted(self._static_mods.items()))
        if self._labile_mods is not None:
            self._labile_mods = dict(sorted(self._labile_mods.items()))
        if self._unknown_mods is not None:
            self._unknown_mods = dict(sorted(self._unknown_mods.items()))
        if self._nterm_mods is not None:
            self._nterm_mods = dict(sorted(self._nterm_mods.items()))
        if self._cterm_mods is not None:
            self._cterm_mods = dict(sorted(self._cterm_mods.items()))
        if self._internal_mods is not None:
            for pos in self._internal_mods:
                self._internal_mods[pos] = dict(sorted(self._internal_mods[pos].items()))
            self._internal_mods = dict(sorted(self._internal_mods.items()))
        if self._intervals is not None:
            self._intervals.sort(key=lambda x: (x.start, x.end))

        if self._charge is not None and isinstance(self._charge, list):
            self._charge = list(sorted(self._charge))

        return self

    @property
    def has_sequence(self) -> bool:
        return bool(self._sequence)

    @property
    def has_compound_name(self) -> bool:
        return bool(self._compound_name)

    @property
    def has_ion_name(self) -> bool:
        return bool(self._ion_name)

    @property
    def has_peptide_name(self) -> bool:
        return bool(self._peptide_name)

    @property
    def has_isotope_mods(self) -> bool:
        return bool(self._isotope_mods)

    @property
    def has_static_mods(self) -> bool:
        return bool(self._static_mods)

    @property
    def has_labile_mods(self) -> bool:
        return bool(self._labile_mods)

    @property
    def has_unknown_mods(self) -> bool:
        return bool(self._unknown_mods)

    @property
    def has_nterm_mods(self) -> bool:
        return bool(self._nterm_mods)

    @property
    def has_cterm_mods(self) -> bool:
        return bool(self._cterm_mods)

    @property
    def has_internal_mods(self) -> bool:
        return bool(self._internal_mods)

    @property
    def has_intervals(self) -> bool:
        return bool(self._intervals)

    @property
    def has_charge(self) -> bool:
        if isinstance(self._charge, list):
            return len(self._charge) > 0
        elif isinstance(self._charge, int):
            return self._charge != 0
        return self._charge is not None

    def _has_mods_by_type(self, mod_type: ModType) -> bool:
        match mod_type:
            case ModType.ISOTOPE:
                return self.has_isotope_mods
            case ModType.STATIC:
                return self.has_static_mods
            case ModType.LABILE:
                return self.has_labile_mods
            case ModType.UNKNOWN:
                return self.has_unknown_mods
            case ModType.NTERM:
                return self.has_nterm_mods
            case ModType.CTERM:
                return self.has_cterm_mods
            case ModType.INTERNAL:
                return self.has_internal_mods
            case ModType.INTERVAL:
                return self.has_intervals
            case ModType.CHARGE:
                return self.has_charge
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def has_mods(
        self,
        mod_types: (Iterable[ModTypeLiteral] | Iterable[ModType] | ModType | ModTypeLiteral | None) = None,
    ) -> bool:
        """Return ``True`` if any of the specified modification types are present.

        :param mod_types: Types to check; all types when ``None``.
        :type mod_types: Iterable[ModTypeLiteral] | Iterable[ModType] | ModType | ModTypeLiteral | None
        :return: ``True`` if at least one matching modification exists.
        :rtype: bool
        """
        mod_enums = get_mods(mod_types)
        return any(self._has_mods_by_type(mod_enum) for mod_enum in mod_enums)

    @property
    def has_sequence_ambiguity(self) -> bool:
        return self.has_intervals or self.has_unknown_mods

    @property
    def has_residue_ambiguity(self) -> bool:
        return len(self.ambiguous_residues) > 0

    @property
    def ambiguous_residues(self) -> tuple[str, ...]:
        return tuple(aa for aa in self.stripped_sequence if AA_LOOKUP.is_ambiguous(aa))

    @property
    def has_mass_ambiguity(self) -> bool:
        return len(self.mass_ambiguous_residues) > 0

    @property
    def mass_ambiguous_residues(self) -> tuple[str, ...]:
        return tuple(aa for aa in self.stripped_sequence if AA_LOOKUP.is_mass_ambiguous(aa))

    def _get_mods_by_type(self, mod_type: ModType) -> Any:
        match mod_type:
            case ModType.ISOTOPE:
                return self.isotope_mods
            case ModType.STATIC:
                return self.static_mods
            case ModType.LABILE:
                return self.labile_mods
            case ModType.UNKNOWN:
                return self.unknown_mods
            case ModType.NTERM:
                return self.nterm_mods
            case ModType.CTERM:
                return self.cterm_mods
            case ModType.INTERNAL:
                return self.internal_mods
            case ModType.INTERVAL:
                return self.intervals
            case ModType.CHARGE:
                return self.charge
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def get_mods(
        self,
        mod_types: (Iterable[ModTypeLiteral] | Iterable[ModType] | ModType | ModTypeLiteral | None) = None,
    ) -> dict[ModType | ModTypeLiteral, Any]:
        """Return a dict of present modification types mapped to their values.

        Only types that currently have modifications are included in the result.

        :param mod_types: Types to include; all types when ``None``.
        :type mod_types: Iterable[ModTypeLiteral] | Iterable[ModType] | ModType | ModTypeLiteral | None
        :return: Mapping of mod type to modification value.
        :rtype: dict[ModType | ModTypeLiteral, Any]
        """
        mod_enums = get_mods(mod_types)
        return {mod_enum: self._get_mods_by_type(mod_enum) for mod_enum in mod_enums if self._has_mods_by_type(mod_enum)}

    @classmethod
    def parse_chimeric(cls, sequence: str, validate: bool | None = None) -> Generator["ProFormaAnnotation", None, None]:
        """Parse a ProForma string into multiple ProFormaAnnotation objects"""
        if validate is None:
            validate = False
        # Initialize the Generator
        parser_gen = ProFormaParser(sequence).parse()

        for prof_parser, _ in parser_gen:
            # Split Global Mods into Static (Fixed) vs Isotope (Global)
            # The parser groups all <...> tags together; we separate them by the '@' symbol.
            static_mods: dict[str, int] | None = None  # e.g., <[Oxidation]@C>
            isotope_mods: dict[str, int] | None = None  # e.g., <13C>

            if prof_parser.global_mods:
                for mod, count in prof_parser.global_mods.items():
                    if "@" in mod:
                        if static_mods is None:
                            static_mods = {}
                        static_mods[mod] = count
                    else:
                        if isotope_mods is None:
                            isotope_mods = {}
                        isotope_mods[mod] = count

            charge = None
            if prof_parser.charge is not None:
                charge = prof_parser.charge
            elif prof_parser.charge_adducts is not None:
                charge = prof_parser.charge_adducts

            # Construct the object
            # We cast defaultdicts to standard dicts to prevent side effects
            annot = ProFormaAnnotation(
                sequence="".join(prof_parser.amino_acids),
                compound_name=prof_parser.compound_name,
                ion_name=prof_parser.ion_name,
                peptide_name=prof_parser.peptide_name,
                isotope_mods=isotope_mods,
                static_mods=static_mods,
                labile_mods=prof_parser.labile_mods,
                unknown_mods=prof_parser.unknown_mods,
                nterm_mods=prof_parser.nterm_mods,
                cterm_mods=prof_parser.cterm_mods,
                internal_mods=prof_parser.internal_mods,
                intervals=prof_parser.intervals,
                charge=charge,
                validate=validate,
            )

            yield annot

    @classmethod
    def parse(cls, sequence: str, validate: bool | None = None) -> "ProFormaAnnotation":
        """Parse a ProForma string into a ProFormaAnnotation object"""
        if validate is None:
            validate = False
        # Initialize the Generator
        parser_gen = ProFormaParser(sequence).parse()

        # Get first annotation segment
        try:
            prof_parser, connection = next(parser_gen)
        except StopIteration as e:
            raise ValueError(f"Invalid ProForma sequence: {sequence}") from e

        # Validate that this is a single peptide (not chimeric/crosslinked)
        if connection is not None:
            raise ValueError(f"Chimeric and crosslinked peptides not supported in single annotation: {sequence}")

        # Ensure there are no subsequent segments waiting in the generator
        try:
            next(parser_gen)
            raise ValueError(f"Multiple peptide segments found in sequence: {sequence}")
        except StopIteration:
            pass  # This is expected for a single annotation

        # Split Global Mods into Static (Fixed) vs Isotope (Global)
        # The parser groups all <...> tags together; we separate them by the '@' symbol.
        static_mods: dict[str, int] | None = None  # e.g., <[Oxidation]@C>
        isotope_mods: dict[str, int] | None = None  # e.g., <13C>

        if prof_parser.global_mods:
            for mod, count in prof_parser.global_mods.items():
                if "@" in mod:
                    if static_mods is None:
                        static_mods = {}
                    static_mods[mod] = count
                else:
                    if isotope_mods is None:
                        isotope_mods = {}
                    isotope_mods[mod] = count

        charge = None
        if prof_parser.charge is not None:
            charge = prof_parser.charge
        elif prof_parser.charge_adducts is not None:
            charge = prof_parser.charge_adducts

            def convert_charge_count(cnt: int) -> str:
                if cnt <= 0:
                    raise ValueError("Charge count cannot be less than or equal to zero.")
                elif cnt == 1:
                    return ""
                else:
                    return f"^{int(cnt)}"

            charge = [f"{adduct}{convert_charge_count(count)}" for adduct, count in charge.items()]

        # Construct the object
        # We cast defaultdicts to standard dicts to prevent side effects
        annot = ProFormaAnnotation(
            sequence="".join(prof_parser.amino_acids),
            compound_name=prof_parser.compound_name,
            ion_name=prof_parser.ion_name,
            peptide_name=prof_parser.peptide_name,
            isotope_mods=isotope_mods,
            static_mods=static_mods,
            labile_mods=prof_parser.labile_mods,
            unknown_mods=prof_parser.unknown_mods,
            nterm_mods=prof_parser.nterm_mods,
            cterm_mods=prof_parser.cterm_mods,
            internal_mods=prof_parser.internal_mods,
            intervals=prof_parser.intervals,
            charge=charge,
            validate=validate,
        )

        return annot

    def serialize(self, exclude_charge: bool = False) -> str:
        """Serialise this annotation to a ProForma string.

        :param exclude_charge: If ``True``, omit the charge suffix from the output.
        :type exclude_charge: bool
        :return: A ProForma-compliant string representation.
        :rtype: str
        """
        return serialize_annotation(self, exclude_charge=exclude_charge)

    def serialize_charge(self) -> str:
        return serialize_charge(self)

    def get_sequence_composition(self) -> Counter[ElementInfo]:
        sequence_composition: Counter[ElementInfo] = Counter()
        for aa in self.stripped_sequence:
            aa_info = AA_LOOKUP.one_letter(aa)
            if aa_info.composition is None:
                raise ValueError(f"Composition not available for amino acid: {aa}")
            for element, count in aa_info.composition.items():
                sequence_composition[element] += count
        return sequence_composition

    @property
    def stripped_sequence(self) -> str:
        """Get the unmodified amino acid sequence without any modifications"""
        return self._sequence if self._sequence is not None else ""

    def map_static_mods_to_indexes(self) -> dict[int, list[Mod[ModificationTags]]]:
        # can be N, C or internal
        if self._static_mods is None:
            return {}

        mapped_mods: dict[int, list[Mod[ModificationTags]]] = {}
        for static_mod in self.static_mods:
            valid_indexes = static_mod.value.find_indexes(self.sequence)

            if not valid_indexes:
                continue

            for index in valid_indexes:
                if index not in mapped_mods:
                    mapped_mods[index] = []
                mapped_mods[index].append(Mod(static_mod.value.modifications, count=static_mod.count))
        return mapped_mods

    def _map_isotopes(self) -> dict[ElementInfo, ElementInfo]:
        """Map isotope modifications to element replacements."""
        isotope_map: dict[ElementInfo, ElementInfo] = {}
        if not self.has_isotope_mods:
            return isotope_map

        for isotope_mod in self._isotope_mods.keys():  # type: ignore
            template, replaced = IsotopeReplacement.from_string(isotope_mod).get_isotope_replacements()
            isotope_map[template] = replaced

        return isotope_map

    def _base_comp(self, skip_labile: bool = False, monoisotopic: bool = True) -> tuple[Counter[ElementInfo], int, float]:
        total_composition: Counter[ElementInfo] = self.get_sequence_composition()
        total_charge = 0  # results from internal formula mods
        total_delta_mass = 0.0  # only from MassTags

        if self.has_unknown_mods:
            unknown_mods = self.unknown_mods
            composition, delta_mass, charge = unknown_mods.get_composition_with_delta_mass_charge(monoisotopic)
            total_composition += composition
            total_delta_mass += delta_mass
            total_charge += charge

        if not skip_labile and self.has_labile_mods:
            labile_mods = self.labile_mods
            composition, delta_mass, charge = labile_mods.get_composition_with_delta_mass_charge(monoisotopic)
            total_composition += composition
            total_delta_mass += delta_mass
            total_charge += charge

        if self.has_nterm_mods:
            nterm_mods = self.nterm_mods
            composition, delta_mass, charge = nterm_mods.get_composition_with_delta_mass_charge(monoisotopic)
            total_composition += composition
            total_delta_mass += delta_mass
            total_charge += charge

        if self.has_cterm_mods:
            cterm_mods = self.cterm_mods
            composition, delta_mass, charge = cterm_mods.get_composition_with_delta_mass_charge(monoisotopic)
            total_composition += composition
            total_delta_mass += delta_mass
            total_charge += charge

        if self.has_static_mods:
            static_mod_map = self.map_static_mods_to_indexes()
            for _, mods in static_mod_map.items():
                for mod in mods:
                    try:
                        total_composition += mod.get_composition()
                    except ValueError as e:
                        if isinstance(mod.value, ModificationTags) and isinstance(mod.value.first_tag, TagMass):
                            # MassTag does not have composition, only delta mass
                            total_delta_mass += mod.get_mass(monoisotopic)
                        else:
                            raise e
                    total_charge += mod.get_charge()

        # Internal mods
        if self.has_internal_mods:
            for mods in self.internal_mods.values():
                composition, delta_mass, charge = mods.get_composition_with_delta_mass_charge(monoisotopic)
                total_composition += composition
                total_delta_mass += delta_mass
                total_charge += charge

        # Intervals
        if self.has_intervals:
            for interval in self.intervals:
                composition, delta_mass, charge = interval.mods.get_composition_with_delta_mass_charge(monoisotopic)
                total_composition += composition
                total_delta_mass += delta_mass
                total_charge += charge

        return total_composition, total_charge, total_delta_mass

    def comp(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        charge: CHARGE_TYPE | None = None,
        *,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
    ) -> Counter[ElementInfo]:
        """Calculate composition, preferring user charge over annotation charge."""

        frag = self.frag(
            ion_type=ion_type,
            charge=charge,
            monoisotopic=True,
            isotopes=isotopes,
            deltas=deltas,
            calculate_composition=True,
            _include_sequence=False,
        )

        if frag.composition is None:
            raise ValueError("Fragment composition could not be calculated.")

        return frag.composition

    def _base_mass(self, monoisotopic: bool = True, skip_labile: bool = False) -> tuple[float, int]:
        """Optimized mass calculation with minimal overhead."""
        total_mass = 0.0
        total_charge = 0  # results from internal formula mods

        # Inline mass lookup to avoid function call overhead
        # Amino acids - hot path, optimize heavily
        aa_lookup = AA_LOOKUP.one_letter_to_info
        for aa in self.stripped_sequence:
            mass = aa_lookup[aa].get_mass(monoisotopic=monoisotopic)
            if mass is None:
                raise ValueError(f"Mass not available for amino acid: {aa}")
            total_mass += mass

        # Unknown mods
        if self.has_unknown_mods:
            m, c = self.unknown_mods.get_mass_charge(monoisotopic=monoisotopic)
            total_mass += m
            total_charge += c

        # Labile mods
        if not skip_labile and self.has_labile_mods:
            m, c = self.labile_mods.get_mass_charge(monoisotopic=monoisotopic)
            total_mass += m
            total_charge += c

        # N-terminal mods
        if self.has_nterm_mods:
            m, c = self.nterm_mods.get_mass_charge(monoisotopic=monoisotopic)
            total_mass += m
            total_charge += c

        # Internal mods
        if self.has_internal_mods:
            for mods in self.internal_mods.values():
                m, c = mods.get_mass_charge(monoisotopic=monoisotopic)
                total_mass += m
                total_charge += c

        # Interval mods
        if self.has_intervals:
            for interval in self.intervals:
                m, c = interval.mods.get_mass_charge(monoisotopic=monoisotopic)
                total_mass += m
                total_charge += c

        # C-terminal mods
        if self.has_cterm_mods:
            m, c = self.cterm_mods.get_mass_charge(monoisotopic=monoisotopic)
            total_mass += m
            total_charge += c

        # Static mods
        if self.has_static_mods:
            static_mod_map = self.map_static_mods_to_indexes()
            for mods in static_mod_map.values():
                for mod in mods:
                    total_mass += mod.get_mass(monoisotopic=monoisotopic)
                    total_charge += mod.get_charge()

        return total_mass, total_charge

    def _get_mass_vector(self, monoisotopic: bool = True) -> list[float]:
        # slice sequence into single aa slcices
        vec: list[float] = []
        for i in range(len(self)):
            sub_annot = self.slice(i, i + 1, inplace=False)
            vec.append(
                sub_annot.mass(
                    ion_type=IonType.NEUTRAL,
                    charge=None,
                    monoisotopic=monoisotopic,
                    isotopes=None,
                )
            )
        return vec

    def _build_mass_vector(self, monoisotopic: bool = True) -> list[float]:
        """Build a per-residue mass array without sequence slicing.

        Each element is the residue mass plus any modifications localised to that
        position (internal mods, N-/C-terminal mods at the respective ends, and
        static mods mapped to their target residues).

        :param monoisotopic: Use monoisotopic masses when ``True``, average masses when ``False``.
        :type monoisotopic: bool
        :return: List of per-residue masses, length == len(self).
        :rtype: list[float]
        :raises ValueError: If the annotation contains unknown mods or interval mods.
        """
        if self.has_unknown_mods or self.has_intervals:
            raise ValueError(f"fragment_masses not supported for sequences with unknown modifications or intervals: {str(self)}")

        aa_lookup = AA_LOOKUP.one_letter_to_info
        masses: list[float] = []
        for aa in self.stripped_sequence:
            m = aa_lookup[aa].get_mass(monoisotopic=monoisotopic)
            if m is None:
                raise ValueError(f"Mass not available for amino acid: {aa}")
            masses.append(m)

        if self.has_nterm_mods:
            m, _ = self.nterm_mods.get_mass_charge(monoisotopic=monoisotopic)
            masses[0] += m

        if self.has_internal_mods:
            for pos, mods in self.internal_mods.items():
                m, _ = mods.get_mass_charge(monoisotopic=monoisotopic)
                masses[pos] += m

        if self.has_cterm_mods:
            m, _ = self.cterm_mods.get_mass_charge(monoisotopic=monoisotopic)
            masses[-1] += m

        if self.has_static_mods:
            static_mod_map = self.map_static_mods_to_indexes()
            for pos, mods_list in static_mod_map.items():
                for mod in mods_list:
                    masses[pos] += mod.get_mass(monoisotopic=monoisotopic)

        return masses

    def _get_comp_vector(self) -> list[Counter[ElementInfo]]:
        # slice sequence into single aa slcices
        vec: list[Counter[ElementInfo]] = []
        for i in range(len(self)):
            sub_annot = self.slice(i, i + 1, inplace=False)
            vec.append(
                sub_annot.comp(
                    ion_type=IonType.NEUTRAL,
                    charge=None,
                    isotopes=None,
                )
            )
        return vec

    @property
    def monoisotopic_base_mass(self) -> float:
        """Calculate monoisotopic mass of the unmodified sequence."""
        return self._base_mass(monoisotopic=True)[0]

    @property
    def average_base_mass(self) -> float:
        """Calculate average mass of the unmodified sequence."""
        return self._base_mass(monoisotopic=False)[0]

    def mass(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        charge: CHARGE_TYPE | None = None,
        monoisotopic: bool = True,
        *,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        calculate_with_composition: bool = False,
    ) -> float:
        """Calculate mass, preferring user charge over annotation charge."""

        f = self.frag(
            ion_type=ion_type,
            charge=charge,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_composition=calculate_with_composition,
            _include_sequence=False,
        )
        return f.mass

    def neutral_mass(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        monoisotopic: bool = True,
        *,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        calculate_with_composition: bool = False,
    ) -> float:
        """Calculate the neutral (uncharged) mass for the given ion type.

        :param ion_type: Fragment ion type to use for the calculation.
        :type ion_type: ION_TYPE
        :param monoisotopic: Use monoisotopic masses when ``True``, average masses when ``False``.
        :type monoisotopic: bool
        :param isotopes: Isotope offsets or element-count overrides.
        :type isotopes: ISOTOPE_TYPE | None
        :param deltas: Custom neutral-loss or gain formula(s).
        :type deltas: CUSTOM_LOSS_TYPE | None
        :param calculate_with_composition: Derive mass from elemental composition instead of
            the fast mass-lookup path.
        :type calculate_with_composition: bool
        :return: Neutral mass in daltons.
        :rtype: float
        """
        f = self.frag(
            ion_type=ion_type,
            charge=0,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_composition=calculate_with_composition,
        )
        return f.mass

    def _frag(
        self,
        ion_type: IonType,
        monoisotopic: bool,
        isotope: IsotopeInfo,
        delta: DeltaInfo,
        calculate_composition: bool,
        parent_sequence: str,
        parent_sequence_length: int,
        position: int | tuple[int, int] | None,
    ) -> Fragment:
        # Dont include labile mods for fragment ions
        skip_labile = True
        if ion_type == IonType.NEUTRAL or ion_type == IonType.PRECURSOR:
            skip_labile = False

        if self.has_isotope_mods or calculate_composition:
            base_comp, base_charge, delta_mass = self._base_comp(skip_labile=skip_labile, monoisotopic=monoisotopic)

            if calculate_composition and delta_mass != 0.0:
                raise ValueError("Cannot calculate composition with delta mass changes.")
            elif calculate_composition and delta_mass == 0.0:
                return adjust_comp(
                    base_comp=base_comp,
                    charge=self.charge_adducts,
                    ion_type=ion_type,
                    monoisotopic=monoisotopic,
                    isotope=isotope,
                    delta=delta,
                    inplace=True,
                    isotope_map=self._map_isotopes() if self.has_isotope_mods else None,
                    position=position,
                    parent_sequence=parent_sequence,
                    parent_sequence_length=parent_sequence_length,
                    internal_charge=base_charge,
                )
            else:
                base_mass = sum(element.get_mass(monoisotopic=monoisotopic) * count for element, count in base_comp.items())
                return adjust_mass_mz(
                    base=base_mass + delta_mass,
                    charge=self.charge_adducts,
                    monoisotopic=monoisotopic,
                    ion_type=ion_type,
                    isotope=isotope,
                    delta=delta,
                    position=position,
                    parent_sequence=parent_sequence,
                    parent_sequence_length=parent_sequence_length,
                    internal_charge=base_charge,
                )

        base_mass, base_charge = self._base_mass(monoisotopic=monoisotopic, skip_labile=skip_labile)

        return adjust_mass_mz(
            base=base_mass,
            charge=self.charge_adducts,
            monoisotopic=monoisotopic,
            ion_type=ion_type,
            isotope=isotope,
            delta=delta,
            position=position,
            parent_sequence=parent_sequence,
            parent_sequence_length=parent_sequence_length,
            internal_charge=base_charge,
        )

    def frag(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        charge: CHARGE_TYPE | None = None,
        monoisotopic: bool = True,
        *,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        calculate_composition: bool = False,
        position: int | tuple[int, int] | None = None,
        _include_sequence: bool = True,
    ) -> Fragment:
        """Calculate mass, preferring user charge over annotation charge."""
        ion_type = IonType(ion_type)
        can_fragment_sequence(self.sequence, ion_type)
        delta_info = DeltaInfo.from_input(deltas)

        inplace = False

        frag_annot = self
        if charge is not None:  # update charge
            frag_annot = frag_annot.set_charge(charge, inplace=inplace)
            inplace = True

        parent_sequence = None
        if _include_sequence:
            parent_sequence = frag_annot.serialize(exclude_charge=False)
        else:
            parent_sequence = ""

        if position is None:
            ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]
            if ion_info.is_forward | ion_info.is_backward:
                position = len(self)  # default to full length for terminal ions
            if ion_info.is_internal:
                if ion_info.ion_type == IonType.IMMONIUM and len(frag_annot) != 1:
                    raise ValueError("Immonium ions must be single amino acids, or the position must be specified.")
                position = (
                    1,
                    len(self),
                )  # use whole sequence for internal ions by default
            if ion_info.is_intact:
                position = None  # use whole sequence for intact ions

        _pos: tuple[int, int] | None = validate_position(ion_type, position, len(self))

        iso_info = IsotopeInfo.from_input(isotopes)

        # get the appropriate fragment annotation based on the position parameter
        match _pos:
            case None:
                pass
            case tuple() as pios_tuple:
                pos_start, pos_end = pios_tuple
                frag_annot = frag_annot.slice(pos_start, pos_end, inplace=inplace)
            case _:
                raise ValueError(f"Invalid position type: {type(position)}")

        return frag_annot._frag(
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=iso_info,
            delta=delta_info,
            calculate_composition=calculate_composition,
            parent_sequence=parent_sequence,
            parent_sequence_length=len(self),
            position=position,
        )

    def mz(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        charge: CHARGE_TYPE | None = None,
        monoisotopic: bool = True,
        *,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        calculate_with_composition: bool = False,
    ) -> float:
        """Calculate m/z, preferring user charge over annotation charge."""

        f = self.frag(
            ion_type=ion_type,
            charge=charge,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_composition=calculate_with_composition,
            _include_sequence=False,
        )
        return f.mz

    def _fragment(
        self,
        ion_type: IonType,
        monoisotopic: bool = True,
        *,
        isotopes: list[IsotopeInfo],
        deltas: list[DeltaInfo],
        neutral_deltas: list[NeutralDeltaInfo],
        calculate_composition: bool,
        parent_sequence: str,
        parent_sequence_length: int,
        max_deltas: int,
        min_length: int | None,
        max_length: int | None,
    ) -> Generator[Fragment, None, None]:
        if self.has_unknown_mods or self.has_intervals:
            raise ValueError(f"Fragmentation not supported for sequences with unknown modifications or intervals: {str(self)}")

        # check for fragments wa/wb and da/db, then return both
        match ion_type:
            case IonType.W:
                # yeild both wa and wb
                yield from self._fragment(
                    IonType.WA,
                    monoisotopic,
                    isotopes=isotopes,
                    deltas=deltas,
                    neutral_deltas=neutral_deltas,
                    calculate_composition=calculate_composition,
                    parent_sequence=parent_sequence,
                    parent_sequence_length=parent_sequence_length,
                    max_deltas=max_deltas,
                    min_length=min_length,
                    max_length=max_length,
                )
                yield from self._fragment(
                    IonType.WB,
                    monoisotopic,
                    isotopes=isotopes,
                    deltas=deltas,
                    neutral_deltas=neutral_deltas,
                    calculate_composition=calculate_composition,
                    parent_sequence=parent_sequence,
                    parent_sequence_length=parent_sequence_length,
                    max_deltas=max_deltas,
                    min_length=min_length,
                    max_length=max_length,
                )
                return
            case IonType.D:
                # yeild both da and db
                yield from self._fragment(
                    IonType.DA,
                    monoisotopic,
                    isotopes=isotopes,
                    deltas=deltas,
                    neutral_deltas=neutral_deltas,
                    calculate_composition=calculate_composition,
                    parent_sequence=parent_sequence,
                    parent_sequence_length=parent_sequence_length,
                    max_deltas=max_deltas,
                    min_length=min_length,
                    max_length=max_length,
                )
                yield from self._fragment(
                    IonType.DB,
                    monoisotopic,
                    isotopes=isotopes,
                    deltas=deltas,
                    neutral_deltas=neutral_deltas,
                    calculate_composition=calculate_composition,
                    parent_sequence=parent_sequence,
                    parent_sequence_length=parent_sequence_length,
                    max_deltas=max_deltas,
                    min_length=min_length,
                    max_length=max_length,
                )
                return
            case _:
                pass

        ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]

        loss_dict: dict[NeutralDeltaInfo, int] = {}
        # Forward ions: b1, b2, b3, ... (cumulative from N-terminus)
        if ion_info.is_forward:
            for i in range(1, len(self) + 1):
                if min_length is not None and i < min_length:
                    continue

                if max_length is not None and i > max_length:
                    break

                sub_annot = self.slice(0, i, inplace=False)

                try:
                    ion_type = can_fragment_sequence(sub_annot.sequence, ion_type)
                except ValueError:
                    continue

                if neutral_deltas:
                    loss_dict.clear()
                    for nd in neutral_deltas:
                        loss_dict[nd] = min(nd.calculate_loss_sites(sub_annot.sequence), max_deltas)

                neutral_delta_combinations: list[DeltaInfo] = get_loss_combinations(loss_dict, max_deltas)

                for isotope in isotopes:
                    for delta in deltas:
                        for ndelta in neutral_delta_combinations:
                            combined_delta = delta + ndelta
                            yield sub_annot._frag(
                                ion_type=ion_type,
                                monoisotopic=monoisotopic,
                                isotope=isotope,
                                delta=combined_delta,
                                calculate_composition=calculate_composition,
                                parent_sequence=parent_sequence,
                                parent_sequence_length=parent_sequence_length,
                                position=i,  # Changed: position is the cleavage site
                            )

        # Backward ions: y1, y2, y3, ... (cumulative from C-terminus)
        elif ion_info.is_backward:
            for i in range(1, len(self) + 1):
                if min_length is not None and i < min_length:
                    continue

                if max_length is not None and i > max_length:
                    break

                sub_annot = self[len(self) - i : len(self)]

                try:
                    ion_type = can_fragment_sequence(sub_annot.sequence, ion_type)
                except ValueError:
                    continue

                if neutral_deltas:
                    loss_dict.clear()
                    for nd in neutral_deltas:
                        loss_dict[nd] = min(nd.calculate_loss_sites(sub_annot.sequence), max_deltas)

                neutral_delta_combinations = get_loss_combinations(loss_dict, max_deltas)

                for isotope in isotopes:
                    for delta in deltas:
                        for ndelta in neutral_delta_combinations:
                            combined_delta = delta + ndelta

                            yield sub_annot._frag(
                                ion_type=ion_type,
                                monoisotopic=monoisotopic,
                                isotope=isotope,
                                delta=combined_delta,
                                calculate_composition=calculate_composition,
                                parent_sequence=parent_sequence,
                                parent_sequence_length=parent_sequence_length,
                                position=i,
                            )

        elif ion_info.is_intact:
            if min_length is not None and len(self) < min_length:
                return

            if max_length is not None and len(self) > max_length:
                return

            if neutral_deltas:
                loss_dict.clear()
                for nd in neutral_deltas:
                    loss_dict[nd] = min(nd.calculate_loss_sites(self.sequence), max_deltas)

            neutral_delta_combinations = get_loss_combinations(loss_dict, max_deltas)

            for isotope in isotopes:
                for delta in deltas:
                    for ndelta in neutral_delta_combinations:
                        combined_delta = delta + ndelta
                        yield self._frag(
                            ion_type=ion_type,
                            monoisotopic=monoisotopic,
                            isotope=isotope,
                            delta=combined_delta,
                            calculate_composition=calculate_composition,
                            parent_sequence=parent_sequence,
                            parent_sequence_length=parent_sequence_length,
                            position=len(self),  # Precursor position
                        )
        elif ion_info.is_internal:
            if ion_info.ion_type == IonType.IMMONIUM:
                # Immonium ions are single residue fragments
                for i in range(1, len(self) + 1):
                    sub_annot = self.slice(i - 1, i, inplace=False)

                    if neutral_deltas:
                        loss_dict.clear()
                        for nd in neutral_deltas:
                            loss_dict[nd] = min(nd.calculate_loss_sites(self.sequence), max_deltas)

                    neutral_delta_combinations = get_loss_combinations(loss_dict, max_deltas)

                    for isotope in isotopes:
                        for delta in deltas:
                            for ndelta in neutral_delta_combinations:
                                combined_delta = delta + ndelta
                                yield sub_annot._frag(
                                    ion_type=ion_type,
                                    monoisotopic=monoisotopic,
                                    isotope=isotope,
                                    delta=combined_delta,
                                    calculate_composition=calculate_composition,
                                    parent_sequence=parent_sequence,
                                    parent_sequence_length=parent_sequence_length,
                                    position=i,  # Position is the residue index
                                )
            else:
                # gen all internal fragmetns from 1 to n-1
                for start in range(2, len(self)):  # Start from position 1 to len-1
                    for end in range(start, len(self)):  # End before C-terminus
                        # Apply length filters
                        if min_length is not None and (end - start + 1) < min_length:
                            continue
                        if max_length is not None and (end - start + 1) > max_length:
                            continue
                        sub_annot = self.slice(start - 1, end, inplace=False)

                        if neutral_deltas:
                            loss_dict.clear()
                            for nd in neutral_deltas:
                                loss_dict[nd] = min(nd.calculate_loss_sites(self.sequence), max_deltas)

                        neutral_delta_combinations = get_loss_combinations(loss_dict, max_deltas)

                        for isotope in isotopes:
                            for delta in deltas:
                                for ndelta in neutral_delta_combinations:
                                    combined_delta = delta + ndelta
                                    yield sub_annot._frag(
                                        ion_type=ion_type,
                                        monoisotopic=monoisotopic,
                                        isotope=isotope,
                                        delta=combined_delta,
                                        calculate_composition=calculate_composition,
                                        parent_sequence=parent_sequence,
                                        parent_sequence_length=parent_sequence_length,
                                        position=(
                                            start,
                                            end,
                                        ),
                                    )

    def fragment(
        self,
        ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
        charges: Sequence[CHARGE_TYPE] | None = None,
        monoisotopic: bool = True,
        *,
        isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
        deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
        neutral_deltas: Sequence[LOSS_TYPE] = (),
        calculate_composition: bool = False,
        max_ndeltas: int = 1,
        min_length: int | None = None,
        max_length: int | None = None,
    ) -> list[Fragment]:
        """Generate fragment annotation for given ion type."""

        if charges is None:
            cstate = self.charge_state
            if cstate != 0:
                charges = (cstate,)
            else:
                charges = (1,)

        # charge_infos: list[ChargeCarrierInfo] = [ChargeCarrierInfo.from_input(charge) for charge in charges]

        isotope_infos: list[IsotopeInfo] = [IsotopeInfo.from_input(isotope) for isotope in isotopes]

        neutral_deltas_infos: list[NeutralDeltaInfo] = []
        if neutral_deltas_infos is not None:
            for loss in neutral_deltas:
                if loss is None:
                    continue
                if isinstance(loss, NeutralDeltaInfo):
                    neutral_deltas_infos.append(loss)
                else:
                    nd: NeutralDeltaInfo = NEUTRAL_DELTA_LOOKUP[loss]
                    neutral_deltas_infos.append(nd)

        delta_infos = [DeltaInfo.from_input(loss) for loss in deltas]

        fragments: list[Fragment] = []
        for charge in charges:
            charged_annot = self.set_charge(charge, inplace=False)
            sequence = charged_annot.serialize()
            for ion in ion_types:
                fragments.extend(
                    list(
                        charged_annot._fragment(
                            ion_type=IonType(ion),
                            monoisotopic=monoisotopic,
                            isotopes=isotope_infos,
                            deltas=delta_infos,
                            neutral_deltas=neutral_deltas_infos,
                            calculate_composition=calculate_composition,
                            parent_sequence=sequence,
                            parent_sequence_length=len(charged_annot),
                            max_deltas=max_ndeltas,
                            min_length=min_length,
                            max_length=max_length,
                        )
                    )
                )
        return fragments

    _UNSUPPORTED_META_ION_TYPES: frozenset[IonType] = frozenset({IonType.W, IonType.D})

    def fragment_masses(
        self,
        ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
        charges: Sequence[int] = (1,),
        monoisotopic: bool = True,
    ) -> dict[tuple[IonType, int], list[float]]:
        """Compute fragment ion m/z values using a fast prefix/suffix-sum approach.

        Returns a mapping of ``(ion_type, charge)`` to a list of m/z values of
        length ``len(self)``, ordered from fragment position 1 to N. No neutral
        losses, isotope shifts, or custom deltas are applied.

        :param ion_types: Ion series to generate (e.g. ``IonType.B``, ``IonType.Y``).
            Meta-types such as ``IonType.W`` and ``IonType.D`` are not supported;
            pass their concrete sub-types (``WA``/``WB``, ``DA``/``DB``) directly.
        :type ion_types: Sequence[ION_TYPE]
        :param charges: Proton charge states to compute m/z for.
        :type charges: Sequence[int]
        :param monoisotopic: Use monoisotopic masses when ``True``, average masses when ``False``.
        :type monoisotopic: bool
        :return: Dict mapping ``(IonType, charge)`` to a list of m/z values.
        :rtype: dict[tuple[IonType, int], list[float]]
        :raises ValueError: If a meta ion type (``W``, ``D``) is requested, or if the
            annotation contains unknown mods or interval mods.
        """
        for ion_type_input in ion_types:
            if IonType(ion_type_input) in self._UNSUPPORTED_META_ION_TYPES:
                raise ValueError(
                    f"Meta ion type {ion_type_input!r} is not supported in fragment_masses(). "
                    "Pass concrete sub-types (e.g. IonType.WA, IonType.WB) directly."
                )

        n = len(self)
        mass_vec = self._build_mass_vector(monoisotopic=monoisotopic)
        result: dict[tuple[IonType, int], list[float]] = {}

        for charge in charges:
            charge_offset = charge * PROTON_MASS
            for ion_type_input in ion_types:
                ion_type = IonType(ion_type_input)
                ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]
                ion_offset = ion_info.get_mass(monoisotopic=monoisotopic)

                if ion_info.is_forward:
                    prefix = 0.0
                    masses_out: list[float] = []
                    for m in mass_vec:
                        prefix += m
                        masses_out.append((prefix + ion_offset + charge_offset) / charge)
                    result[(ion_type, charge)] = masses_out

                elif ion_info.is_backward:
                    prefix = 0.0
                    masses_out = []
                    for i in range(n - 1, -1, -1):
                        prefix += mass_vec[i]
                        masses_out.append((prefix + ion_offset + charge_offset) / charge)
                    result[(ion_type, charge)] = masses_out

                elif ion_info.is_intact:
                    total = sum(mass_vec)
                    mz = (total + ion_offset + charge_offset) / charge
                    result[(ion_type, charge)] = [mz] * n

                elif ion_info.is_internal:
                    result[(ion_type, charge)] = [
                        (mass_vec[i] + ion_offset + charge_offset) / charge
                        for i in range(n)
                    ]

        return result

    """
    Pop Methods
    """

    def pop_isotope_mods(self, inplace: bool = True) -> Mods[IsotopeReplacement]:
        """Pop and return isotope modifications, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed isotope modifications (empty ``Mods`` if none were present).
        :rtype: Mods[IsotopeReplacement]
        """
        if not self.has_isotope_mods:
            return EMPTY_ISOTOPE_MODS

        if not inplace:
            return self.copy().pop_isotope_mods(inplace=True)

        value = self.isotope_mods
        self._isotope_mods = None
        return value

    def pop_static_mods(self, inplace: bool = True) -> Mods[FixedModification]:
        """Pop and return static modifications, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed static modifications (empty ``Mods`` if none were present).
        :rtype: Mods[FixedModification]
        """
        if not self.has_static_mods:
            return EMPTY_STATIC_MODS

        if not inplace:
            return self.copy().pop_static_mods(inplace=True)

        value = self.static_mods
        self._static_mods = None
        return value

    def pop_labile_mods(self, inplace: bool = True) -> Mods[ModificationTags]:
        """Pop and return labile modifications, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed labile modifications (empty ``Mods`` if none were present).
        :rtype: Mods[ModificationTags]
        """
        if not self.has_labile_mods:
            return EMPTY_LABILE_MODS

        if not inplace:
            return self.copy().pop_labile_mods(inplace=True)

        value = self.labile_mods
        self._labile_mods = None
        return value

    def pop_unknown_mods(self, inplace: bool = True) -> Mods[ModificationTags]:
        """Pop and return unknown-localisation modifications, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed unknown modifications (empty ``Mods`` if none were present).
        :rtype: Mods[ModificationTags]
        """
        if not self.has_unknown_mods:
            return EMPTY_UNKNOWN_MODS

        if not inplace:
            return self.copy().pop_unknown_mods(inplace=True)

        value = self.unknown_mods
        self._unknown_mods = None
        return value

    def pop_nterm_mods(self, inplace: bool = True) -> Mods[ModificationTags]:
        """Pop and return N-terminal modifications, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed N-terminal modifications (empty ``Mods`` if none were present).
        :rtype: Mods[ModificationTags]
        """
        if not self.has_nterm_mods:
            return EMPTY_NTERM_MODS

        if not inplace:
            return self.copy().pop_nterm_mods(inplace=True)

        value = self.nterm_mods
        self._nterm_mods = None
        return value

    def pop_cterm_mods(self, inplace: bool = True) -> Mods[ModificationTags]:
        """Pop and return C-terminal modifications, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed C-terminal modifications (empty ``Mods`` if none were present).
        :rtype: Mods[ModificationTags]
        """
        if not self.has_cterm_mods:
            return EMPTY_CTERM_MODS

        if not inplace:
            return self.copy().pop_cterm_mods(inplace=True)

        value = self.cterm_mods
        self._cterm_mods = None
        return value

    def pop_internal_mods(self, inplace: bool = True) -> dict[int, Mods[ModificationTags]]:
        """Pop and return all internal modifications, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed internal modifications (empty dict if none were present).
        :rtype: dict[int, Mods[ModificationTags]]
        """
        if not self.has_internal_mods:
            return {}

        if not inplace:
            return self.copy().pop_internal_mods(inplace=True)

        value = self.internal_mods
        self._internal_mods = None
        return value

    def pop_intervals(self, inplace: bool = True) -> list[Interval]:
        """Pop and return all intervals, clearing them from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed intervals (empty list if none were present).
        :rtype: list[Interval]
        """
        if not self.has_intervals:
            return []

        if not inplace:
            return self.copy().pop_intervals(inplace=True)

        value = self._intervals.copy() if self._intervals else []
        self._intervals = None
        return value

    def pop_charge(self, inplace: bool = True) -> int | Mods[GlobalChargeCarrier] | None:
        """Pop and return the charge, clearing it from the annotation.

        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: The removed charge value, or ``None`` if no charge was set.
        :rtype: int | Mods[GlobalChargeCarrier] | None
        """
        if not self.has_charge:
            return None

        if not inplace:
            return self.copy().pop_charge(inplace=True)

        value = self.charge
        self._charge = None
        return value

    def pop_internal_mod_at_index(self, index: int, inplace: bool = True) -> tuple[tuple[MODIFICATION_TYPE, int], ...]:
        """Pop and return internal modifications at a single 0-based position.

        :param index: 0-based residue index.
        :type index: int
        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: Tuple of ``(modification, count)`` pairs; empty tuple if none were present.
        :rtype: tuple[tuple[MODIFICATION_TYPE, int], ...]
        """
        if self._internal_mods is None:
            return ()

        if index not in self._internal_mods:
            return ()

        if not inplace:
            return self.copy().pop_internal_mod_at_index(index, inplace=True)

        # Parse the modifications at this index before removing
        mods_dict = self._internal_mods[index]
        mods = tuple((ModificationTags.from_string(mod_str), count) for mod_str, count in mods_dict.items())

        # Remove the mod dict at this index
        del self._internal_mods[index]

        # Clean up if internal_mods is now empty
        if len(self._internal_mods) == 0:
            self._internal_mods = None

        return mods

    def _pop_mod_by_type(self, mod_type: ModType) -> Any:
        match mod_type:
            case ModType.ISOTOPE:
                return self.pop_isotope_mods(inplace=True)
            case ModType.STATIC:
                return self.pop_static_mods(inplace=True)
            case ModType.LABILE:
                return self.pop_labile_mods(inplace=True)
            case ModType.UNKNOWN:
                return self.pop_unknown_mods(inplace=True)
            case ModType.NTERM:
                return self.pop_nterm_mods(inplace=True)
            case ModType.CTERM:
                return self.pop_cterm_mods(inplace=True)
            case ModType.INTERNAL:
                return self.pop_internal_mods(inplace=True)
            case ModType.INTERVAL:
                return self.pop_intervals(inplace=True)
            case ModType.CHARGE:
                return self.pop_charge(inplace=True)
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def pop_mods(
        self,
        mod_types: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None = None,
        inplace: bool = True,
    ) -> dict[ModType, Any]:
        """Pop and return modifications of the specified types, clearing them from the annotation.

        :param mod_types: Types to pop; all types when ``None``.
        :type mod_types: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None
        :param inplace: Clear from this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: Mapping of :class:`ModType` to the removed modification values.
        :rtype: dict[ModType, Any]
        """
        if inplace is False:
            return self.copy().pop_mods(mod_types=mod_types, inplace=True)

        mod_enums: list[ModType] = get_mods(mod_types)

        d: dict[ModType, Any] = {}
        for mod_enum in mod_enums:
            d[mod_enum] = self._pop_mod_by_type(mod_enum)

        return d

    def filter_mods(
        self,
        mods: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None = None,
        inplace: bool = True,
        keep: bool = True,
    ) -> Self:
        """Filter modifications by type, either keeping or removing the specified types.

        :param mods: Modification types to keep or remove; all types when ``None``.
        :type mods: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param keep: When ``True`` keep only the specified types; when ``False`` remove them.
        :type keep: bool
        :return: The (possibly new) annotation with filtered modifications.
        :rtype: Self
        """
        if inplace is False:
            return self.copy().filter_mods(mods=mods, inplace=True)

        if keep:
            # Keep only specified mods
            mod_types_to_keep = set(get_mods(mods))

            all_mod_types = {mod_type for mod_type in ModType}
            mod_types_to_remove = all_mod_types - mod_types_to_keep
        else:
            # Remove only specified mods
            mod_types_to_remove = get_mods(mods)

        if len(mod_types_to_remove) == 0:
            # If no mods to remove, return the annotation as is
            return self
        self.pop_mods(mod_types_to_remove)
        return self

    """
    Remove Methods
    """

    def _clear_mod_dict(self, attr_name: str, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy()._clear_mod_dict(attr_name, inplace=True)
        setattr(self, attr_name, None)
        return self

    def clear_isotope_mods(self, inplace: bool = True) -> Self:
        """Clear all isotope modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_isotope_mods", inplace)

    def clear_static_mods(self, inplace: bool = True) -> Self:
        """Clear all static modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_static_mods", inplace)

    def clear_nterm_mods(self, inplace: bool = True) -> Self:
        """Clear all N-terminal modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_nterm_mods", inplace)

    def clear_cterm_mods(self, inplace: bool = True) -> Self:
        """Clear all C-terminal modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_cterm_mods", inplace)

    def clear_labile_mods(self, inplace: bool = True) -> Self:
        """Clear all labile modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_labile_mods", inplace)

    def clear_unknown_mods(self, inplace: bool = True) -> Self:
        """Clear all unknown-localisation modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_unknown_mods", inplace)

    def clear_internal_mods(self, inplace: bool = True) -> Self:
        """Clear all internal (per-position) modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_internal_mods", inplace)

    def clear_internal_mod_at_index(self, index: int, inplace: bool = True) -> Self:
        """Clear internal modifications at a single 0-based sequence position.

        :param index: 0-based residue index.
        :type index: int
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        if not inplace:
            return self.copy().clear_internal_mod_at_index(index, inplace=True)
        if self._internal_mods is None or index not in self._internal_mods:
            return self
        del self._internal_mods[index]
        if len(self._internal_mods) == 0:
            self._internal_mods = None
        return self

    def clear_intervals(self, inplace: bool = True) -> Self:
        """Clear all ambiguous sequence intervals.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_intervals", inplace)

    def clear_charge(self, inplace: bool = True) -> Self:
        """Clear the charge value.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_charge", inplace)

    def _clear_mod_by_type(self, mod_type: ModType) -> None:
        match mod_type:
            case ModType.ISOTOPE:
                self.clear_isotope_mods(inplace=True)
            case ModType.STATIC:
                self.clear_static_mods(inplace=True)
            case ModType.LABILE:
                self.clear_labile_mods(inplace=True)
            case ModType.UNKNOWN:
                self.clear_unknown_mods(inplace=True)
            case ModType.NTERM:
                self.clear_nterm_mods(inplace=True)
            case ModType.CTERM:
                self.clear_cterm_mods(inplace=True)
            case ModType.INTERNAL:
                self.clear_internal_mods(inplace=True)
            case ModType.INTERVAL:
                self.clear_intervals(inplace=True)
            case ModType.CHARGE:
                self.clear_charge(inplace=True)
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def clear_mods(
        self,
        mods: (ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None) = None,
        inplace: bool = True,
    ) -> Self:
        """Clear modifications of the specified types (all types when ``None``).

        :param mods: Types to clear; all types when ``None``.
        :type mods: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation with selected modifications cleared.
        :rtype: Self
        """
        if inplace is False:
            return self.copy().clear_mods(mods=mods, inplace=True)
        mod_enums = get_mods(mods)
        for mod_enum in mod_enums:
            self._clear_mod_by_type(mod_enum)
        return self

    def strip_mods(self, inplace: bool = False) -> Self:
        """Remove all modifications of every type, leaving only the bare sequence.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) bare annotation.
        :rtype: Self
        """
        return self.clear_mods(None, inplace=inplace)

    """
    Slicing Methods
    """

    def slice_by_span(
        self,
        span: Span | tuple[int, int, int],
        inplace: bool = False,
    ) -> Self:
        return self.slice(span[0], span[1], inplace=inplace)

    def slice(
        self,
        start: int | None,
        stop: int | None,
        inplace: bool = False,
    ) -> Self:
        """Return a sub-annotation spanning ``sequence[start:stop]``, carrying over applicable mods.

        :param start: 0-based start index (inclusive); ``None`` means the beginning.
        :type start: int | None
        :param stop: 0-based stop index (exclusive); ``None`` means the end.
        :type stop: int | None
        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :return: The sliced annotation.
        :rtype: Self
        """
        return cast(
            Self,
            slice_annotation(
                self,
                start=start,
                stop=stop,
                inplace=inplace,
            ),
        )

    def split(self) -> list[Self]:
        """Split this annotation into a list of single-residue annotations.

        :return: List of single-residue annotations in sequence order.
        :rtype: list[Self]
        """
        return [cast(Self, a) for a in split_annotation(self)]

    @staticmethod
    def join(annotations: Sequence["ProFormaAnnotation"]) -> "ProFormaAnnotation":
        return join_annotations(annotations)

    def shift(
        self,
        n: int,
        keep_nterm: int = 0,
        keep_cterm: int = 0,
        inplace: bool = False,
    ) -> Self:
        """Cyclically shift the sequence by *n* positions, optionally anchoring termini.

        :param n: Number of positions to shift (positive = rightward).
        :type n: int
        :param keep_nterm: Number of N-terminal residues to keep in place.
        :type keep_nterm: int
        :param keep_cterm: Number of C-terminal residues to keep in place.
        :type keep_cterm: int
        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :return: The (possibly new) shifted annotation.
        :rtype: Self
        """
        return cast(Self, shift_annotation(self, n, keep_nterm, keep_cterm, inplace))

    def shuffle(
        self,
        seed: Any = None,
        keep_nterm: int = 0,
        keep_cterm: int = 0,
        inplace: bool = False,
    ) -> Self:
        """Randomly shuffle the sequence residues, optionally anchoring termini.

        :param seed: Random seed for reproducibility; ``None`` for a random shuffle.
        :type seed: Any
        :param keep_nterm: Number of N-terminal residues to keep in place.
        :type keep_nterm: int
        :param keep_cterm: Number of C-terminal residues to keep in place.
        :type keep_cterm: int
        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :return: The (possibly new) shuffled annotation.
        :rtype: Self
        """
        return cast(Self, shuffle_annotation(self, seed, keep_nterm, keep_cterm, inplace))

    def reverse(
        self,
        keep_nterm: int = 0,
        keep_cterm: int = 0,
        inplace: bool = False,
    ) -> Self:
        """Reverse the sequence residues, optionally anchoring termini.

        :param keep_nterm: Number of N-terminal residues to keep in place.
        :type keep_nterm: int
        :param keep_cterm: Number of C-terminal residues to keep in place.
        :type keep_cterm: int
        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :return: The (possibly new) reversed annotation.
        :rtype: Self
        """
        return cast(Self, reverse_annotation(self, keep_nterm, keep_cterm, inplace))

    def sort(
        self,
        inplace: bool = False,
        key: Callable[[str], Any] | None = None,
        reverse: bool = False,
    ) -> Self:
        """Sort the sequence residues, optionally with a custom key.

        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :param key: Key function applied to each residue for sorting; default alphabetical.
        :type key: Callable[[str], Any] | None
        :param reverse: If ``True``, sort in descending order.
        :type reverse: bool
        :return: The (possibly new) sorted annotation.
        :rtype: Self
        """
        return cast(Self, sort_annotation(self, inplace, key, reverse))

    def sliding_windows(
        self,
        window_size: int,
        reverse: bool = False,
    ) -> Generator[Self, None, None]:
        """Yield overlapping sub-annotations of a fixed window size.

        :param window_size: Number of residues in each window.
        :type window_size: int
        :param reverse: Iterate from C-terminus to N-terminus when ``True``.
        :type reverse: bool
        :return: Generator of window annotations.
        :rtype: Generator[Self, None, None]
        """
        for window in generate_sliding_windows(self, window_size, reverse):
            yield cast(Self, window)

    """
    Modification Methods
    """

    def condense_static_mods(self, inplace: bool = True) -> Self:
        return cast(Self, condense_static_mods(self, inplace=inplace))

    def condense_to_peptidoform(self, inplace: bool = True) -> Self:
        return cast(Self, condense_to_peptidoform(self, inplace=inplace))

    def count_residues(self, include_mods: bool = True) -> dict[str, int]:
        return count_residues(self, include_mods=include_mods)

    def percent_residues(self, include_mods: bool = True) -> dict[str, float]:
        return percent_residues(self, include_mods=include_mods)

    def is_subsequence(
        self,
        other: Self,
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> bool:
        return is_subsequence(self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals)

    def find_indices(
        self,
        other: Self,
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> list[int]:
        return find_indices(self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals)

    def condense_mods_to_intervals(self, inplace: bool = True) -> Self:
        return cast(Self, condense_mods_to_intervals(self, inplace=inplace))

    def coverage(
        self,
        annotations: Iterable[Self],
        accumulate: bool = False,
        ignore_mods: bool = False,
        ignore_ambiguity: bool = False,
    ) -> list[int]:
        return coverage(
            annotation=self,
            annotations=annotations,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )

    def percent_coverage(
        self,
        annotations: Iterable[Self],
        accumulate: bool = False,
        ignore_mods: bool = False,
        ignore_ambiguity: bool = False,
    ) -> float:
        return percent_coverage(
            annotation=self,
            annotations=annotations,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )

    def modification_coverage(
        self,
        annotations: Iterable[Self],
        ignore_ambiguity: bool = False,
        accumulate: bool = False,
    ) -> dict[int, int]:
        return modification_coverage(
            annotation=self,
            annotations=annotations,
            ignore_ambiguity=ignore_ambiguity,
            accumulate=accumulate,
        )

    def permutations(self, size: int | None = None) -> Generator[Self, None, None]:
        for item in generate_permutations(self, size):
            yield cast(Self, item)

    def product(self, repeat: int | None = None) -> Generator[Self, None, None]:
        for item in generate_product(self, repeat):
            yield cast(Self, item)

    def combinations(self, r: int | None = None) -> Generator[Self, None, None]:
        for item in generate_combinations(self, r):
            yield cast(Self, item)

    def combinations_with_replacement(self, r: int | None = None) -> Generator[Self, None, None]:
        for item in generate_combinations_with_replacement(self, r):
            yield cast(Self, item)

    def modify(
        self,
        *,
        nterm_static: Mapping[str | None, Iterable[Any]] | None = None,
        cterm_static: Mapping[str | None, Iterable[Any]] | None = None,
        internal_static: Mapping[str | None, Iterable[Any]] | None = None,
        labile_static: Mapping[str | None, Iterable[Any]] | None = None,
        nterm_variable: Mapping[str | None, Iterable[Any]] | None = None,
        cterm_variable: Mapping[str | None, Iterable[Any]] | None = None,
        internal_variable: (Mapping[str | None, Iterable[Any]] | None) = None,
        labile_variable: Mapping[str | None, Iterable[Any]] | None = None,
        max_variable_mods: int = 2,
        use_regex: bool = False,
        inplace: bool = False,
        use_static_notation: bool = False,
        unique_peptidoforms: bool = True,
    ) -> Generator[Self, None, None]:
        """
        Build all modifications from intervals and mass shifts.
        """
        for annot in modify(
            self,
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
            inplace=inplace,
            use_static_notation=use_static_notation,
            unique_peptidoforms=unique_peptidoforms,
        ):
            yield cast(Self, annot)

    def add_static_mod_by_residue(
        self,
        residue: str | Iterable[str],
        mod: Any,
        inplace: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_static_mod_by_residue(residue, mod, inplace=True)

        residues = list(residue)

        mod_str, count = convert_single_mod_input(mod)

        if count != 1:
            raise ValueError("Fixed modifications added by residue must have a count of 1.")

        # filter residues to only those in the sequence
        residues = [aa for aa in residues if aa in self.stripped_sequence]

        if len(residues) == 0:
            return self

        rules: list[PositionRule] = []
        for aa in residues:
            rules.append(PositionRule(terminal=Terminal.ANYWHERE, amino_acid=AminoAcid.from_str(aa)))
        fixed_mod = FixedModification(
            modifications=ModificationTags.from_string(mod_str),
            position_rules=tuple(rules),
        )

        self.append_static_mod(fixed_mod, inplace=True)
        return self

    def get_interval(self, start: int, end: int) -> Interval | None:
        """Get the interval modification that spans the given start and end positions.

        Args:
            start (int): The start position of the interval (1-based).
            end (int): The end position of the interval (1-based).
        Returns:
            Interval | None: The Interval object if found, otherwise None.
        """
        if not self.has_intervals:
            return None

        for interval in self.intervals:
            if interval.start == start and interval.end == end:
                return interval

        return None

    def annotate_ambiguity(
        self,
        forward_coverage: list[int],
        reverse_coverage: list[int],
        mass_shift: Any | None = None,
        add_mods_to_intervals: bool = False,
        sort_mods: bool = True,
        inplace: bool = False,
    ) -> Self:
        """Annotate modification-site ambiguity using forward and reverse fragment coverage vectors.

        :param forward_coverage: Per-position coverage counts from N-terminal fragments.
        :type forward_coverage: list[int]
        :param reverse_coverage: Per-position coverage counts from C-terminal fragments.
        :type reverse_coverage: list[int]
        :param mass_shift: Optional mass shift to associate with ambiguous intervals.
        :type mass_shift: Any | None
        :param add_mods_to_intervals: If ``True``, include modifications in the interval objects.
        :type add_mods_to_intervals: bool
        :param sort_mods: If ``True``, sort modifications after annotation.
        :type sort_mods: bool
        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotated annotation.
        :rtype: Self
        """
        return cast(
            Self,
            annotate_ambiguity(
                self,
                forward_coverage=forward_coverage,
                reverse_coverage=reverse_coverage,
                mass_shift=mass_shift,
                add_mods_to_intervals=add_mods_to_intervals,
                sort_mods=sort_mods,
                inplace=inplace,
            ),
        )

    def condense_ambiguity_to_xnotation(self, inplace: bool = True) -> Self:
        """Condense ambiguous interval regions to X-notation placeholders.

        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation with intervals condensed to X notation.
        :rtype: Self
        """
        return cast(Self, condense_ambiguity_to_xnotation(self, inplace=inplace))

    @staticmethod
    def group_by_ambiguity(annotations: Iterable["ProFormaAnnotation"], precision: int = 5) -> list[tuple["ProFormaAnnotation", ...]]:
        """Group annotations that are ambiguous equivalents of each other.

        :param annotations: Annotations to group.
        :type annotations: Iterable[ProFormaAnnotation]
        :param precision: Decimal precision for mass comparison when grouping.
        :type precision: int
        :return: List of groups, each group being a tuple of equivalent annotations.
        :rtype: list[tuple[ProFormaAnnotation, ...]]
        """
        return group_by_ambiguity(annotations, precision=precision)

    @staticmethod
    def unique_fragments(annotations: Iterable["ProFormaAnnotation"], precision: int = 4) -> list[int]:
        """Return the indices of annotations that produce unique fragment masses.

        :param annotations: Annotations to compare.
        :type annotations: Iterable[ProFormaAnnotation]
        :param precision: Decimal precision for fragment mass comparison.
        :type precision: int
        :return: List of indices into *annotations* whose fragment masses are unique.
        :rtype: list[int]
        """
        return unique_fragments(annotations, precision=precision)

    def to_ip2(self) -> str:
        """Convert this annotation to IP2 format string.

        :raises NotImplementedError: This conversion is not yet implemented.
        """
        raise NotImplementedError("Conversion to IP2 format is not yet implemented.")

    @staticmethod
    def from_ip2_sequence(sequence: str) -> "ProFormaAnnotation":
        """Internal function for converting a single IP2 sequence."""
        if re.match(r"^([A-Z]|-)\..*\.([A-Z]|-)$", sequence):
            sequence = sequence[2:-2]
        sequence = re.sub(r"\(([^)]+)\)", r"[\1]", sequence)
        sequence = re.sub(r"^\[([^\]]+)\]", r"[\1]-", sequence)
        sequence = re.sub(r"\]\[", r"]-[", sequence)

        return ProFormaAnnotation.parse(sequence)

    def to_diann(self) -> str:
        """Convert this annotation to DIA-NN format string.

        :raises NotImplementedError: This conversion is not yet implemented.
        """
        raise NotImplementedError("Conversion to DIANN format is not yet implemented.")

    @staticmethod
    def from_diann(sequence: str) -> "ProFormaAnnotation":
        """Internal function for converting a single DIANN sequence."""
        if sequence.startswith("_"):
            sequence = sequence[1:]
            if re.match(r"^\[[^\]]+\]", sequence):
                sequence = re.sub(r"^\[([^\]]+)\]", r"[\1]-", sequence)

        if sequence.endswith("_"):
            sequence = sequence[:-1]

        elif re.search(r"_\[[^\]]+\]$", sequence):
            sequence = re.sub(r"_\[([^\]]+)\]$", r"-[\1]", sequence)

        return ProFormaAnnotation.parse(sequence)

    def to_casanovo(self) -> str:
        """Convert this annotation to Casanovo format string.

        :raises NotImplementedError: This conversion is not yet implemented.
        """
        raise NotImplementedError("Conversion to Casanovo format is not yet implemented.")

    @staticmethod
    def from_casanovo(sequence: str) -> "ProFormaAnnotation":
        """Internal function for converting a single Casanovo sequence."""
        new_sequence_comps: list[str] = []
        in_mod = False  # Tracks if we are within a modification
        is_nterm = False  # Tracks if the current modification is at the N-terminus

        for _, char in enumerate(sequence):
            if char in {"+", "-"}:
                # Check if it's at the start (N-terminal)
                is_nterm = len(new_sequence_comps) == 0

                # Start a new modification block
                new_sequence_comps.append("[")
                new_sequence_comps.append(char)
                in_mod = True
            elif in_mod and char.isalpha():
                # End the modification block
                new_sequence_comps.append("]")

                if is_nterm:
                    # Add a dash if it's an N-terminal modification
                    new_sequence_comps.append("-")
                    is_nterm = False

                # Add the current character and close modification
                in_mod = False
                new_sequence_comps.append(char)
            else:
                # Add regular characters
                new_sequence_comps.append(char)

        # Close any unclosed modification at the end of the sequence
        if in_mod:
            new_sequence_comps.append("]")

        sequence = "".join(new_sequence_comps)
        return ProFormaAnnotation.parse(sequence)

    def to_ms2_pip(
        self,
        inplace: bool = False,
    ) -> tuple[str, str]:
        """Convert a single peptide sequence to MS2PIP format

        Returns:
            tuple[str, str]: (unmodified_sequence, modification_string)
                where modification_string is in format "loc1|name1|loc2|name2|..."
        """

        if self.has_mods(
            (
                ModType.ISOTOPE,
                ModType.LABILE,
                ModType.UNKNOWN,
                ModType.INTERVAL,
                ModType.CHARGE,
            )
        ):
            raise ValueError("MS2PIP format does not support isotope, labile, unknown, interval, charge, or charge adduct modifications.")

        if not inplace:
            # Create a copy to condense
            annot_copy = self.copy()
            annot_copy.condense_static_mods(inplace=True)
        else:
            self.condense_static_mods(inplace=True)
            annot_copy = self

        mod_tuples: list[tuple[int, str]] = []

        # Process N-terminal modifications
        if annot_copy._nterm_mods is not None:
            for mod_name, count in annot_copy._nterm_mods.items():
                if count != 1:
                    raise ValueError("MS2PIP format does not support modification multipliers.")
                mod_tuples.append((0, mod_name))

        # Process C-terminal modifications
        if annot_copy._cterm_mods is not None:
            for mod_name, count in annot_copy._cterm_mods.items():
                if count != 1:
                    raise ValueError("MS2PIP format does not support modification multipliers.")
                mod_tuples.append((-1, mod_name))

        # Process internal modifications
        if annot_copy._internal_mods is not None:
            for index, mods_dict in annot_copy._internal_mods.items():
                if len(mods_dict) > 1:
                    raise ValueError("MS2PIP format does not support multiple modifications at the same site.")
                for mod_name, count in mods_dict.items():
                    if count != 1:
                        raise ValueError("MS2PIP format does not support modification multipliers.")
                    # MS2PIP uses 1-indexed positions
                    mod_tuples.append((index + 1, mod_name))

        unmod_sequence = annot_copy.stripped_sequence

        mod_str = "|".join(f"{loc}|{name}" for loc, name in mod_tuples)

        return unmod_sequence, mod_str

    @staticmethod
    def from_ms2_pip(
        sequence: str,
        mod_str: str,
        static_mods: Mapping[str, float | int | str] | None = None,
    ) -> "ProFormaAnnotation":
        """Create ProFormaAnnotation from MS2PIP format"""

        # Create annotation with just the sequence
        annot = ProFormaAnnotation(sequence=sequence)

        if mod_str.strip() == "":
            # No modifications - just add static mods if provided
            for static_aa, mass in (static_mods or {}).items():
                annot.add_static_mod_by_residue(static_aa, mass, inplace=True)
            return annot

        # Parse modification string: format is "loc1|name1|loc2|name2|..."
        mod_parts = mod_str.split("|")

        if len(mod_parts) % 2 != 0:
            raise ValueError(f"Invalid MS2PIP modification string format: {mod_str}")

        # Process modifications in pairs (location, name)
        for i in range(0, len(mod_parts), 2):
            loc_str = mod_parts[i]
            mod_name = mod_parts[i + 1]

            # Parse location
            loc = int(loc_str)

            # Add to appropriate location
            if loc == 0:
                # N-terminal modification
                annot.append_nterm_mod(mod_name, inplace=True)
            elif loc == -1:
                # C-terminal modification
                annot.append_cterm_mod(mod_name, inplace=True)
            else:
                # Internal modification (1-indexed in MS2PIP, convert to 0-indexed)
                internal_index = loc - 1
                if internal_index < 0 or internal_index >= len(sequence):
                    raise ValueError(f"Modification location {loc} is out of range for sequence of length {len(sequence)}")
                annot.append_internal_mod_at_index(internal_index, mod_name, inplace=True)

        # Add static modifications
        for static_aa, mass in (static_mods or {}).items():
            annot.add_static_mod_by_residue(static_aa, mass, inplace=True)

        return annot

    def isotopic_distribution(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        charge: CHARGE_TYPE | None = None,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        max_isotopes: int | None = 10,
        min_abundance_threshold: float = 0.001,  # based on the most abundant peak
        distribution_resolution: int | None = 5,
        use_neutron_count: bool = False,
        conv_min_abundance_threshold: float = 10e-15,
    ) -> list[IsotopicData]:
        """Calculate the exact isotopic distribution from elemental composition.

        :param ion_type: Fragment ion type to use.
        :type ion_type: ION_TYPE
        :param charge: Charge override; uses the annotation charge when ``None``.
        :type charge: CHARGE_TYPE | None
        :param isotopes: Isotope offsets or element-count overrides.
        :type isotopes: ISOTOPE_TYPE | None
        :param deltas: Custom neutral-loss or gain formula(s).
        :type deltas: CUSTOM_LOSS_TYPE | None
        :param max_isotopes: Maximum number of isotope peaks to return.
        :type max_isotopes: int | None
        :param min_abundance_threshold: Minimum relative abundance (vs. the most abundant peak).
        :type min_abundance_threshold: float
        :param distribution_resolution: Decimal places used when binning m/z values.
        :type distribution_resolution: int | None
        :param use_neutron_count: Use neutron count instead of exact mass offsets.
        :type use_neutron_count: bool
        :param conv_min_abundance_threshold: Minimum absolute abundance during convolution.
        :type conv_min_abundance_threshold: float
        :return: List of isotopic data points sorted by m/z.
        :rtype: list[IsotopicData]
        """
        # check if any deltas provided are float?

        frag_annot = self
        if charge is not None:  # update charge
            frag_annot = frag_annot.set_charge(charge, inplace=False)

        composition: Counter[ElementInfo] = frag_annot.comp(ion_type=ion_type, isotopes=isotopes, deltas=deltas)

        return isotopic_distribution(
            chemical_formula=cast(Mapping[str | ElementInfo, int | float], composition),
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
            distribution_resolution=distribution_resolution,
            use_neutron_count=use_neutron_count,
            conv_min_abundance_threshold=conv_min_abundance_threshold,
            charge_state=frag_annot.charge_state,
        )

    def estimate_isotopic_distribution(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        charge: CHARGE_TYPE | None = None,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        max_isotopes: int | None = 10,
        min_abundance_threshold: float = 0.001,
        distribution_resolution: int | None = 5,
        use_neutron_count: bool = False,
        conv_min_abundance_threshold: float = 10e-15,
    ) -> list[IsotopicData]:
        """Estimate isotopic distribution based on mass."""

        mass = self.mass(ion_type=ion_type, charge=charge, isotopes=isotopes, deltas=deltas)

        return estimate_isotopic_distribution(
            neutral_mass=mass,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
            distribution_resolution=distribution_resolution,
            use_neutron_count=use_neutron_count,
            conv_min_abundance_threshold=conv_min_abundance_threshold,
        )

    @staticmethod
    def random(
        min_length: int = 6,
        max_length: int = 20,
        mod_probability: float = 0.05,
        include_internal_mods: bool = True,
        include_nterm_mods: bool = True,
        include_cterm_mods: bool = True,
        include_labile_mods: bool = True,
        include_unknown_mods: bool = True,
        include_isotopic_mods: bool = True,
        include_static_mods: bool = True,
        include_intervals: bool = True,
        include_charge: bool = True,
        require_composition: bool = True,
    ) -> "ProFormaAnnotation":
        """Generate a random ProFormaAnnotation for testing purposes."""
        return generate_random_proforma_annotation(
            min_length=min_length,
            max_length=max_length,
            mod_probability=mod_probability,
            include_internal_mods=include_internal_mods,
            include_nterm_mods=include_nterm_mods,
            include_cterm_mods=include_cterm_mods,
            include_labile_mods=include_labile_mods,
            include_unknown_mods=include_unknown_mods,
            include_isotopic_mods=include_isotopic_mods,
            include_static_mods=include_static_mods,
            include_intervals=include_intervals,
            include_charge=include_charge,
            require_composition=require_composition,
        )

    def left_semi_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get left semi-enzymatic sequences (N-terminus fixed)."""
        return left_semi_spans(self, min_len, max_len)

    def right_semi_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get right semi-enzymatic sequences (C-terminus fixed)."""
        return right_semi_spans(self, min_len, max_len)

    def semi_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get all semi-enzymatic sequences."""
        return semi_spans(self, min_len, max_len)

    def nonspecific_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get all non-enzymatic sequences (all possible subsequences)."""
        return nonspecific_spans(self, min_len, max_len)

    def cleavage_sites(
        self,
        enzyme: str | re.Pattern[str],
    ) -> Generator[int, None, None]:
        """Yield 0-based cleavage positions matching the given enzyme regex.

        :param enzyme: Regex pattern (or pre-compiled pattern) defining cleavage sites.
        :type enzyme: str | re.Pattern[str]
        :return: Generator of 0-based indices where cleavage occurs.
        :rtype: Generator[int, None, None]
        """
        # Call the underlying function
        return get_cleavage_sites(self, enzyme)

    def simple_cleavage_sites(
        self,
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
    ) -> Generator[int, None, None]:
        """Get cleavage sites using simple amino acid rules."""
        enzyme_regex = generate_regex(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        )
        return self.cleavage_sites(enzyme_regex)

    def digest(
        self,
        enzyme: str,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Digest this annotation using a regex pattern."""
        return digest_annotation_by_regex(
            annotation=self,
            enzyme_regex=enzyme,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )

    def simple_digest(
        self,
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Digest this annotation with specified enzyme parameters."""
        return digest_annotation_by_aa(
            annotation=self,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )

    def sequential_digest(
        self,
        enzyme_configs: list[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Perform sequential digestion with multiple enzymes."""
        return sequential_digest_annotation(self, enzyme_configs, min_len, max_len)

    @property
    def prop(self) -> AnnotationProperties:
        """Get the properties of this annotation."""
        return AnnotationProperties(self.stripped_sequence)

import logging
import re
from collections import Counter
from collections.abc import Callable, Generator, Iterable, Mapping, Sequence
from types import FunctionType
from typing import (
    TYPE_CHECKING,
    Any,
    Self,
    cast,
)

from tacular import (
    AA_LOOKUP,
    AminoAcid,
    ElementInfo,
    IonType,
    IonTypeLiteral,
    NeutralDelta,
    NeutralDeltaInfo,
    NeutralDeltaLiteral,
)

from ..constants import ModType, Terminal
from ..diagnostics import (
    CompositionError,
    InvalidAdjustmentError,
    PeptacularError,
    ProFormaFormatError,
    UnknownModificationError,
    UnsupportedOperationError,
)
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
    brain_isotopic_distribution,
    estimate_isotopic_distribution,
)
from ..proforma_components import (
    MODIFICATION_TYPE,
    ChargedFormula,
    FixedModification,
    GlobalChargeCarrier,
    IsotopeReplacement,
    ModificationTags,
    PositionRule,
    SequenceElement,
    SequenceRegion,
)
from ..property.prop import AnnotationProperties
from ..spans import Span
from . import _frag_engine, _mass
from ._frag_engine import get_loss_combinations
from ._mass import (
    EMPTY_CHARGE_MODS,
    H_CHARGE_FORMULA,
    H_DECHARGE_FORMULA,
)
from ._mod_access import (
    EMPTY_CTERM_MODS,
    EMPTY_INTERNAL_MODS,
    EMPTY_ISOTOPE_MODS,
    EMPTY_LABILE_MODS,
    EMPTY_NTERM_MODS,
    EMPTY_STATIC_MODS,
    EMPTY_UNKNOWN_MODS,
    ChargeType,
    _ModAccessMixin,
)
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
from .frag import proton_binding_offset
from .localization import DEFAULT_MAX_ISOMERS, candidate_sites, localization_isomers
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
    convert_single_mod_input,
)
from .mod_builder import modify
from .parser import ProFormaParser
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
from .utils import Fragment

if TYPE_CHECKING:
    import numpy as np

# Names that were module attributes here before the frag engine, mass functions and
# modification accessors moved to their own modules. Re-exported so imports from
# ``peptacular.annotation.annotation`` keep working.
from ._frag_engine import (  # noqa: E402, F401
    _ION_TYPE_TO_MZPAF_SERIES,
    ELECTRON_MASS,
    FRAGMENT_ION_LOOKUP,
    NEUTRAL_DELTA_LOOKUP,
    PROTON_MASS,
    SATELLITE_TRIM_END,
    SATELLITE_TRIM_START,
    FragmentIonInfo,
    _as_options,
    _carrier_mass,
    _ion_mass,
    _unless_impossible_loss,
    adjust_comp,
    adjust_mass_mz,
    validate_mass,
    validate_position,
)
from ._mass import (  # noqa: E402, F401
    _AA_COMPOSITIONS,
    _AVERAGE_AA_MASSES,
    _MONOISOTOPIC_AA_MASSES,
    H_ELEMENT_INFO,
    HYDROGEN_BINDING_MASS,
    Element,
    FormulaElement,
    TagMass,
    _adjust_mass_value,
    can_fragment_sequence,
    fe,
    to_ion_type,
)
from ._mod_access import (  # noqa: E402, F401
    InvalidPositionError,
    ModTypeLiteral,
    PositionScore,
    _concrete_position_labels,
    _resolve_mod_types,
    as_mod_iterable,
    convert_moddict_input,
    is_mod_collection,
)

__all__ = [
    "H_CHARGE_FORMULA",
    "H_DECHARGE_FORMULA",
    "ION_TYPE",
    "CHARGE_TYPE",
    "ISOTOPE_TYPE",
    "LOSS_TYPE",
    "CUSTOM_LOSS_TYPE",
    "POSITION_TYPE",
    "EMPTY_ISOTOPE_MODS",
    "EMPTY_STATIC_MODS",
    "EMPTY_UNKNOWN_MODS",
    "EMPTY_LABILE_MODS",
    "EMPTY_NTERM_MODS",
    "EMPTY_CTERM_MODS",
    "EMPTY_CHARGE_MODS",
    "EMPTY_INTERNAL_MODS",
    "ChargeType",
    "get_loss_combinations",
    "ProFormaAnnotation",
]

logger = logging.getLogger(__name__)

# Errors that already say what went wrong; any other ValueError from the parser is a
# notation error and is re-raised as ProFormaFormatError.
_TYPED_ERRORS = (ProFormaFormatError, UnsupportedOperationError, UnknownModificationError, CompositionError, InvalidAdjustmentError)

ION_TYPE = IonTypeLiteral | IonType
CHARGE_TYPE = int | str | list[str] | Mods[GlobalChargeCarrier] | GlobalChargeCarrier | Mod[GlobalChargeCarrier]
ISOTOPE_TYPE = int | dict[str | ElementInfo, int]
LOSS_TYPE = NeutralDelta | NeutralDeltaLiteral | NeutralDeltaInfo | str
CUSTOM_LOSS_TYPE = str | ChargedFormula | float | dict[str | ChargedFormula | float, int]
POSITION_TYPE = int | tuple[int, int]

_SIGNED_DELTA_RE = re.compile(r"^[+-](?:\d+(?:\.\d*)?|\.\d+)$")
_SKYLINE_BARE_NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")
# Skyline's own mass regex is `^[+-]?[0-9]*[.,]?[0-9]*$` (MassModification.cs:108), so a comma
# decimal is a mass to Skyline. Its value depends on the writer's culture, so this matches it
# in order to refuse it, not to read it.
_SKYLINE_COMMA_NUMBER = re.compile(r"[+-]?[0-9]*,[0-9]*")
_SKYLINE_UNIMOD = re.compile(r"unimod:(\d+)", re.IGNORECASE)


def _skyline_require_residue(tokens: list[str], sequence: str) -> None:
    """Raise if no residue has been read yet, so a Skyline modification has nothing to sit on."""
    if not any(not token.startswith("[") for token in tokens):
        raise ProFormaFormatError(f"Modification before the first residue in Skyline sequence: {sequence!r}.")


# `[0-9]` rather than `\\d`: `\\d` matches a full-width or Arabic-Indic digit and `int()`
# accepts it, so a position written with one was read as a number. Philosopher writes the
# position with strconv.Itoa (lib/rep/rep.go:594), which is ASCII.
_MSFRAGGER_ENTRY_RE = re.compile(r"^(?:(N-term)|(C-term)|([0-9]+)([a-zA-Z]))\(([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+))\)$")
# `A-Z` and not `str.isalpha()`, for the reason above: `isalpha()` is true of a non-ASCII
# letter, and the FragPipe `Peptide` column holds bare residues.
_MSFRAGGER_PEPTIDE_RE = re.compile(r"^[A-Z]+$")


class ProFormaAnnotation(_ModAccessMixin):
    """A single ProForma 2.1 peptidoform: a sequence plus its modifications, intervals, names and charge.

    Create one with :meth:`parse` (or :func:`peptacular.parse`) from a ProForma string, or with the
    constructor from parts. Modifications are stored as strings and resolved against tacular only
    when a mass, m/z or composition is requested, so parsing an unknown modification succeeds and
    ``mass()`` raises :class:`UnknownModificationError`.

    The ``set_*`` / ``append_*`` / ``extend_*`` / ``remove_*`` methods take ``inplace`` (default
    ``True``) and return the annotation, so calls can be chained. Because annotations are mutable
    they are not hashable; use :meth:`serialize` as a dict key or set member. ``==`` compares
    contents and ignores the order of modifications at one site.

    Methods ending in ``_spans`` (:meth:`digest_spans`, :meth:`semi_spans`, ...) yield
    :class:`~peptacular.spans.Span` objects; slice with ``annot[span]`` to get the peptide.

    >>> import peptacular as pt
    >>> annot = pt.ProFormaAnnotation.parse("PEM[Oxidation]TIDE/2")
    >>> annot.sequence
    'PEMTIDE'
    >>> annot.charge
    2
    >>> annot.serialize()
    'PEM[Oxidation]TIDE/2'
    """

    def __init__(
        self,
        sequence: str | None = None,
        *,
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
            # Tally identical adduct strings into occurrence counts: two ``'Na:z+1'``
            # entries must become ``{'Na:z+1': 2}`` (Mod scales charge/mass/composition
            # by count), not collapse to a single carrier via a hardcoded count of 1.
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods=dict(Counter(self._charge)))

        return None

    @charge.setter
    def charge(self, value: int | str | list[str] | Mods[GlobalChargeCarrier] | None) -> None:
        self.set_charge(value, inplace=True, validate=self._validate)

    @property
    def charge_state(self) -> int:
        """Numeric charge state derived from the stored charge; 0 when uncharged.

        This is the charge written after the slash, and it does not include any charge carried
        inside a ``Formula`` modification. ProForma 2.1 section 11.5 says that "if both localised
        charges and charge carriers are used all of these should be accounted for in the total
        charge of the peptidoform ion", and gives ``PEPT[Formula:Zn:z+2]IDE/2`` as an example whose
        total charge is +4. For that peptidoform ``charge_state`` is 2, while :meth:`mz` divides by
        the total charge of 4.

        .. code-block:: python

            >>> import peptacular as pt
            >>> annot = pt.ProFormaAnnotation.parse("PEPT[Formula:Zn:z+2]IDE/2")
            >>> annot.charge_state
            2
            >>> round(annot.mass() / annot.mz())
            4

        :rtype: int
        :raises PeptacularError: If the stored charge value has an unexpected type.
        """
        return _mass.charge_state(self)

    @property
    def charge_adducts(self) -> Mods[GlobalChargeCarrier]:
        """Charge expressed as adduct ``Mods``; converts an integer charge to proton/deproton adducts.

        Positive integer charges are converted to proton additions (``[M+nH]n+``);
        negative integer charges are converted to deprotonations (``[M-nH]n-``).

        :rtype: Mods[GlobalChargeCarrier]
        """
        return _mass.charge_adducts(self)

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
            return NotImplemented

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

    # ProFormaAnnotation is mutable (set_charge, append_mods, inplace=True edits ...),
    # and __eq__ compares that mutable state. A hash that changes while the object
    # sits in a set or dict silently corrupts the container, so annotations are
    # unhashable. Use ``annotation.serialize()`` (a str) as a set member or dict key.
    __hash__ = None  # type: ignore[assignment]

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
            intervals=[iv.copy() for iv in self._intervals] if self._intervals is not None else None,
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
        self._intervals = [iv.copy() for iv in other._intervals] if other._intervals is not None else None
        self._charge = other._charge

    def __getitem__(self, key: slice | Span | tuple[int, int, int]) -> Self:
        """Slice by ``annot[start:stop]`` or by a :class:`Span`; returns a new annotation.

        An integer index is not supported, because a residue carries modifications that a
        plain letter cannot. Use ``annot[i:i + 1]`` for a one-residue annotation or
        ``annot.stripped_sequence[i]`` for the letter.
        """
        if isinstance(key, (tuple, Span)):
            return self.slice_by_span(key, inplace=False)
        if isinstance(key, slice):
            start, stop, step = key.start, key.stop, key.step
            if step is not None and step != 1:
                raise PeptacularError("Step slicing not supported")
            return self.slice(start, stop, inplace=False)
        if isinstance(key, int):
            raise UnsupportedOperationError(
                f"ProFormaAnnotation does not support integer indexing (annot[{key}]); "
                f"use annot[{key}:{key + 1}] for a one-residue annotation or annot.stripped_sequence[{key}] for the letter"
            )
        raise TypeError(f"ProFormaAnnotation indices must be slices or Spans, not {type(key).__name__}")

    def sort_mods(self, *, inplace: bool = True) -> Self:
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

    @classmethod
    def parse_chimeric(cls, sequence: str, *, validate: bool | None = None) -> Generator["ProFormaAnnotation", None, None]:
        """Parse a ProForma string into multiple ProFormaAnnotation objects.

        :raises ProFormaFormatError: The string is not valid ProForma.
        :raises UnsupportedOperationError: Cross-linked (``//``) peptidoforms.
        """
        try:
            yield from cls._parse_chimeric(sequence, validate)
        except _TYPED_ERRORS:
            raise
        except ValueError as e:
            raise ProFormaFormatError(str(e)) from e

    @classmethod
    def _parse_chimeric(cls, sequence: str, validate: bool | None = None) -> Generator["ProFormaAnnotation", None, None]:
        if validate is None:
            validate = False
        # Initialize the Generator
        parser_gen = ProFormaParser(sequence).parse()

        for prof_parser, connection in parser_gen:
            if connection is True:
                raise UnsupportedOperationError(f"Cross-linked peptidoforms joined by '//' are not supported (use '+' for chimeric ions): {sequence}")
            yield cls._from_parser(prof_parser, validate)

    @classmethod
    def parse(cls, sequence: str, *, validate: bool | None = None) -> "ProFormaAnnotation":
        """Parse a ProForma string into a ProFormaAnnotation object.

        :raises ProFormaFormatError: The string is not valid ProForma.
        :raises UnsupportedOperationError: Valid ProForma that one annotation cannot hold
            (chimeric ``+`` or cross-linked ``//`` peptidoforms).
        """
        try:
            return cls._parse_single(sequence, validate)
        except _TYPED_ERRORS:
            raise
        except ValueError as e:
            raise ProFormaFormatError(str(e)) from e

    @classmethod
    def _parse_single(cls, sequence: str, validate: bool | None = None) -> "ProFormaAnnotation":
        if validate is None:
            validate = False
        # Initialize the Generator
        parser_gen = ProFormaParser(sequence).parse()

        # Get first annotation segment
        try:
            prof_parser, connection = next(parser_gen)
        except StopIteration as e:
            raise PeptacularError(f"Invalid ProForma sequence: {sequence}") from e

        # Validate that this is a single peptide (not chimeric/crosslinked)
        if connection is not None:
            raise UnsupportedOperationError(f"Chimeric and crosslinked peptides not supported in single annotation: {sequence}")

        # Ensure there are no subsequent segments waiting in the generator
        try:
            next(parser_gen)
            raise PeptacularError(f"Multiple peptide segments found in sequence: {sequence}")
        except StopIteration:
            pass  # This is expected for a single annotation

        return cls._from_parser(prof_parser, validate)

    @classmethod
    def _from_parser(cls, prof_parser: ProFormaParser, validate: bool) -> "ProFormaAnnotation":
        """Build one annotation from one parsed peptidoform (shared by ``parse`` and ``parse_chimeric``)."""
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
                    raise PeptacularError("Charge count cannot be less than or equal to zero.")
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

    def serialize(self, *, exclude_charge: bool = False) -> str:
        """Serialise this annotation to a ProForma string.

        :param exclude_charge: If ``True``, omit the charge suffix from the output.
        :type exclude_charge: bool
        :return: A ProForma-compliant string representation.
        :rtype: str
        """
        return serialize_annotation(self, exclude_charge=exclude_charge)

    def to_dict(self) -> dict[str, Any]:
        """Return a lossless, versioned JSON-compatible representation."""
        from ..proforma_json import to_proforma_dict

        return to_proforma_dict(self)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Self:
        """Restore an annotation from its versioned JSON-compatible representation."""
        from ..proforma_json import from_proforma_dict

        return from_proforma_dict(data, expected_type=cls)

    def to_json(self, *, indent: int | None = None) -> str:
        """Return deterministic JSON text for this annotation."""
        from ..proforma_json import to_proforma_json

        return to_proforma_json(self, indent=indent)

    @classmethod
    def from_json(cls, data: str | bytes | bytearray) -> Self:
        """Restore an annotation from versioned JSON text."""
        from ..proforma_json import from_proforma_json

        return from_proforma_json(data, expected_type=cls)

    def serialize_charge(self) -> str:
        return serialize_charge(self)

    def get_sequence_composition(self) -> Counter[ElementInfo]:
        """Elemental composition of the residues alone (no modifications, no terminal water).

        :rtype: Counter[ElementInfo]
        :raises CompositionError: If a residue has no defined composition: ``B`` or ``Z``, each of
            which stands for two amino acids of different composition, or a bare ``X``, which has
            no composition of its own (write ``X[Formula:...]`` for one, or ``X[+0]`` for zero
            mass, a deliberate departure from ProForma 2.1 section 7.1).
        """
        return _mass.sequence_composition(self)

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

    def map_isotopes(self) -> dict[ElementInfo, ElementInfo]:
        """Map each global isotope label to the element it replaces and its replacement.

        ``<13C>PEPTIDE`` gives ``{C: 13C}``: every carbon atom in the composition is counted
        as carbon-13. The keys and values are tacular ``ElementInfo`` objects (the key is the
        element with no mass number). An annotation without global isotope labels gives ``{}``.

        :return: ``{element: isotope}`` for each global isotope label.
        :rtype: dict[ElementInfo, ElementInfo]

        .. code-block:: python

            >>> import peptacular as pt
            >>> [(str(k.symbol), str(v.symbol), v.mass_number) for k, v in pt.parse("<13C>PEP").map_isotopes().items()]
            [('C', 'C', 13)]
            >>> pt.parse("PEP").map_isotopes()
            {}
        """
        isotope_map: dict[ElementInfo, ElementInfo] = {}
        if not self.has_isotope_mods:
            return isotope_map

        for isotope_mod in self._isotope_mods.keys():  # type: ignore
            template, replaced = IsotopeReplacement.from_string(isotope_mod).get_isotope_replacements()
            isotope_map[template] = replaced

        return isotope_map

    def _base_comp(self, skip_labile: bool = False, monoisotopic: bool = True) -> tuple[Counter[ElementInfo], int, float]:
        """Composition of the residues and modifications, with the internal charge and mass-only delta (see :func:`._mass.base_comp`)."""
        return _mass.base_comp(self, skip_labile=skip_labile, monoisotopic=monoisotopic)

    def comp(
        self,
        charge: CHARGE_TYPE | None = None,
        *,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
    ) -> Counter[ElementInfo]:
        """Calculate composition, preferring user charge over annotation charge."""
        return _mass.comp(self, charge, ion_type=ion_type, isotopes=isotopes, deltas=deltas)

    def _base_mass(self, monoisotopic: bool = True, skip_labile: bool = False) -> tuple[float, int]:
        """Optimized mass calculation with minimal overhead."""
        return _mass.base_mass(self, monoisotopic=monoisotopic, skip_labile=skip_labile)

    def _get_mass_vector(self, monoisotopic: bool = True) -> list[float]:
        """Neutral mass of each one-residue slice (see :func:`._mass.residue_mass_vector`)."""
        return _mass.residue_mass_vector(self, monoisotopic=monoisotopic)

    def _build_mass_vector(self, monoisotopic: bool = True) -> list[float]:
        """Build a per-residue mass array without sequence slicing.

        Each element is the residue mass plus any modifications localised to that
        position (internal mods, N-/C-terminal mods at the respective ends, and
        static mods mapped to their target residues).

        :param monoisotopic: Use monoisotopic masses when ``True``, average masses when ``False``.
        :type monoisotopic: bool
        :return: List of per-residue masses, length == len(self).
        :rtype: list[float]
        :raises PeptacularError: If the annotation contains unknown mods or interval mods.
        """
        return _frag_engine.build_mass_vector(self, monoisotopic=monoisotopic)

    def _get_comp_vector(self) -> list[Counter[ElementInfo]]:
        """Neutral composition of each one-residue slice (see :func:`._mass.residue_comp_vector`)."""
        return _mass.residue_comp_vector(self)

    @property
    def monoisotopic_base_mass(self) -> float:
        """Calculate monoisotopic mass of the unmodified sequence."""
        return _mass.base_mass(self, monoisotopic=True)[0]

    @property
    def average_base_mass(self) -> float:
        """Calculate average mass of the unmodified sequence."""
        return _mass.base_mass(self, monoisotopic=False)[0]

    def mass(
        self,
        charge: CHARGE_TYPE | None = None,
        *,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        monoisotopic: bool = True,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        calculate_with_composition: bool = False,
    ) -> float:
        """Calculate mass, preferring user charge over annotation charge.

        A negative charge is deprotonation, so ``PEPTIDE/-2`` is :math:`[M-2H]^{2-}`. See
        :func:`peptacular.mass` for how that differs from ProForma 2.1 section 11.5, which
        words a negative charge as the addition of electrons.
        """
        return _mass.mass_and_charge(self, ion_type, charge, monoisotopic, isotopes, deltas, calculate_with_composition)[0]

    def _mass_and_charge(
        self,
        ion_type: ION_TYPE,
        charge: CHARGE_TYPE | None,
        monoisotopic: bool,
        isotopes: ISOTOPE_TYPE | None,
        deltas: CUSTOM_LOSS_TYPE | None,
        calculate_with_composition: bool,
    ) -> tuple[float, int]:
        """Avoid fragment allocation and annotation copies for ordinary intact ions."""
        return _mass.mass_and_charge(self, ion_type, charge, monoisotopic, isotopes, deltas, calculate_with_composition)

    def neutral_mass(
        self,
        *,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        monoisotopic: bool = True,
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
        return self.mass(
            ion_type=ion_type,
            charge=0,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_with_composition=calculate_with_composition,
        )

    def _frag(
        self,
        ion_type: IonType,
        monoisotopic: bool,
        isotope: IsotopeInfo,
        delta: DeltaInfo,
        calculate_with_composition: bool,
        parent_sequence: str,
        parent_sequence_length: int,
        position: int | tuple[int, int] | None,
    ) -> Fragment:
        """Build one ion; satellite ions drop the cleaved residue first (see :func:`._frag_engine.frag_one`)."""
        return _frag_engine.frag_one(
            self,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            delta=delta,
            calculate_with_composition=calculate_with_composition,
            parent_sequence=parent_sequence,
            parent_sequence_length=parent_sequence_length,
            position=position,
        )

    def _satellite_mod_error(self, ion_type: IonType) -> str | None:
        """Why a d or w ion of this (sub)sequence is undefined, or None when it is defined.

        A d or w ion keeps part of the cleaved residue's side chain (its beta substituent), so
        it is not defined when that residue carries a modification, explicit or from a global
        fixed modification (as in paftacular). A v ion loses the whole side chain, and its
        modification with it, so v ions are always defined.
        """
        return _frag_engine.satellite_mod_error(self, ion_type)

    def _frag_impl(
        self,
        ion_type: IonType,
        monoisotopic: bool,
        isotope: IsotopeInfo,
        delta: DeltaInfo,
        calculate_with_composition: bool,
        parent_sequence: str,
        parent_sequence_length: int,
        position: int | tuple[int, int] | None,
    ) -> Fragment:
        """Build one ion of this (sub)sequence by mass or composition (see :func:`._frag_engine.frag_impl`)."""
        return _frag_engine.frag_impl(
            self,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            delta=delta,
            calculate_with_composition=calculate_with_composition,
            parent_sequence=parent_sequence,
            parent_sequence_length=parent_sequence_length,
            position=position,
        )

    def frag(
        self,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        charge: CHARGE_TYPE | None = None,
        *,
        monoisotopic: bool = True,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        calculate_with_composition: bool = False,
        position: int | tuple[int, int] | None = None,
        _include_sequence: bool = True,
    ) -> Fragment:
        """Calculate mass, preferring user charge over annotation charge.

        A cross-linked or branched peptidoform, a ``#XL`` or ``#BRANCH`` label on two residues,
        raises :class:`UnsupportedOperationError`: the two residues are bonded to each other and
        that bond is not modelled, so the ions of the linear peptide would be wrong.
        """
        return _frag_engine.frag(
            self,
            ion_type,
            charge,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_with_composition=calculate_with_composition,
            position=position,
            _include_sequence=_include_sequence,
        )

    def mz(
        self,
        charge: CHARGE_TYPE | None = None,
        *,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        monoisotopic: bool = True,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        calculate_with_composition: bool = False,
    ) -> float:
        """Calculate m/z, preferring user charge over annotation charge.

        The denominator is the absolute charge. A negative charge is deprotonation, so
        ``PEPTIDE/-2`` is :math:`[M-2H]^{2-}`; see :func:`peptacular.mz`.
        """
        mass, total_charge = _mass.mass_and_charge(self, ion_type, charge, monoisotopic, isotopes, deltas, calculate_with_composition)
        return mass / abs(total_charge) if total_charge else mass

    def _series_mass_vector(self, monoisotopic: bool, calculate_with_composition: bool) -> list[float] | None:
        """Per-residue masses for the terminal-series fast path, or None when it does not apply.

        Terminal mods sit on the first/last residue, so a prefix sum of length ``i`` equals the
        mass of ``self.slice(0, i)`` (the C-terminal mods only join at ``i == len``) and a suffix
        sum equals ``self[len - i:]``. Anything that needs the composition path, or that slicing
        treats specially, returns None so the caller slices instead.
        """
        return _frag_engine.series_mass_vector(self, monoisotopic, calculate_with_composition)

    def _fragment_series(
        self,
        ion_type: IonType,
        *,
        forward: bool,
        monoisotopic: bool,
        isotopes: list[IsotopeInfo],
        deltas: list[DeltaInfo],
        neutral_deltas: list[NeutralDeltaInfo],
        calculate_with_composition: bool,
        parent_sequence: str,
        parent_sequence_length: int,
        max_deltas: int,
        min_length: int | None,
        max_length: int | None,
        _fast: bool = True,
    ) -> Generator[Fragment, None, None]:
        """Yield one terminal ion series (forward: b1..bn, backward: y1..yn).

        Ions with plain mass adjustments are built from prefix sums without slicing the
        annotation; everything else (composition mode, formula deltas, isotope swaps,
        satellite ions, ...) slices and goes through :meth:`_frag`.
        """
        return _frag_engine.fragment_series(
            self,
            ion_type,
            forward=forward,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            calculate_with_composition=calculate_with_composition,
            parent_sequence=parent_sequence,
            parent_sequence_length=parent_sequence_length,
            max_deltas=max_deltas,
            min_length=min_length,
            max_length=max_length,
            _fast=_fast,
        )

    def _fragment(
        self,
        ion_type: IonType,
        monoisotopic: bool = True,
        *,
        isotopes: list[IsotopeInfo],
        deltas: list[DeltaInfo],
        neutral_deltas: list[NeutralDeltaInfo],
        calculate_with_composition: bool,
        parent_sequence: str,
        parent_sequence_length: int,
        max_deltas: int,
        min_length: int | None,
        max_length: int | None,
        _expand: bool = True,
    ) -> Generator[Fragment, None, None]:
        """Yield every ion of one ion type (see :func:`._frag_engine.fragment_ions`)."""
        return _frag_engine.fragment_ions(
            self,
            ion_type,
            monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            calculate_with_composition=calculate_with_composition,
            parent_sequence=parent_sequence,
            parent_sequence_length=parent_sequence_length,
            max_deltas=max_deltas,
            min_length=min_length,
            max_length=max_length,
            _expand=_expand,
        )

    @staticmethod
    def _default_fragment_charges(charge_state: int) -> tuple[int, ...]:
        """Return default fragment charge states derived from a precursor charge state.

        For a positive precursor charge ``c``, returns ``1, 2, …, c-1``.
        For a negative precursor charge ``c``, returns ``-1, -2, …, c+1``.
        Falls back to ``(1,)`` when ``charge_state`` is 0 (unannotated) or ±1.
        """
        return _frag_engine.default_fragment_charges(charge_state)

    def fragment(
        self,
        ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
        charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
        *,
        monoisotopic: bool = True,
        isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
        deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
        neutral_deltas: Sequence[LOSS_TYPE | None] = (),
        calculate_with_composition: bool = False,
        max_ndeltas: int = 1,
        min_length: int | None = None,
        max_length: int | None = None,
    ) -> list[Fragment]:
        """Generate fragment ions for each ion type and charge.

        A single ion type, charge, isotope or neutral delta may be passed without a list
        (``fragment("by", 2)`` is ``fragment(["by"], [2])``).

        A cross-linked or branched peptidoform, a ``#XL`` or ``#BRANCH`` label on two residues,
        raises :class:`UnsupportedOperationError`: the two residues are bonded to each other and
        that bond is not modelled, so the ions of the linear peptide would be wrong.
        """
        return _frag_engine.fragment(
            self,
            ion_types,
            charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            calculate_with_composition=calculate_with_composition,
            max_ndeltas=max_ndeltas,
            min_length=min_length,
            max_length=max_length,
        )

    def fragment_arrays(
        self,
        ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
        charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
        *,
        monoisotopic: bool = True,
        isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
        deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
        neutral_deltas: Sequence[LOSS_TYPE | None] = (),
        calculate_with_composition: bool = False,
        max_ndeltas: int = 1,
        min_length: int | None = None,
        max_length: int | None = None,
    ) -> "dict[str, np.ndarray]":
        """The ions of :meth:`fragment` as numpy columns, one row per ion (needs ``peptacular[numpy]``).

        Takes the same arguments as :meth:`fragment` and returns the same ions in the same
        order, as a dict of equal-length arrays keyed by :data:`~peptacular.FRAGMENT_ARRAY_KEYS`
        (see :func:`peptacular.fragment_arrays`). ``pl.DataFrame(result)``,
        ``pa.table(result)`` and ``pd.DataFrame(result)`` accept it directly.

        :raises MissingOptionalDependencyError: If numpy is not installed.
        :raises UnsupportedOperationError: If a ``#XL`` or ``#BRANCH`` label is on two residues
            (a cross-linked or branched peptidoform).
        """
        from .frag_arrays import fragment_arrays

        return fragment_arrays(
            [self],
            ion_types,
            charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            calculate_with_composition=calculate_with_composition,
            max_ndeltas=max_ndeltas,
            min_length=min_length,
            max_length=max_length,
        )

    def fast_fragment(
        self, ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y), charges: int | Sequence[int] | None = None, *, monoisotopic: bool = True
    ) -> dict[tuple[IonType, int], list[float]]:
        """Compute fragment ion m/z values using a fast prefix/suffix-sum approach.

        Returns a mapping of ``(ion_type, charge)`` to a list of m/z values of
        length ``len(self)``, ordered from fragment position 1 to N. Annotation
        isotope labels, intrinsic charges, and labile modifications use the
        regular calculation path. No additional losses or deltas are applied.

        :param ion_types: Ion series to generate (e.g. ``IonType.B``, ``IonType.Y``).
            Supports a, b, c, x, y, z, p, and n. Use ``fragment()`` for other series.
        :type ion_types: Sequence[ION_TYPE]
        :param charges: Proton charge states to compute m/z for.  When ``None``,
            defaults to ``1`` through ``precursor_charge - 1`` if the annotation
            carries a positive charge state, otherwise falls back to ``(1,)``.
        :type charges: Sequence[int] or None
        :param monoisotopic: Use monoisotopic masses when ``True``, average masses when ``False``.
        :type monoisotopic: bool
        :return: Dict mapping ``(IonType, charge)`` to a list of m/z values.
        :rtype: dict[tuple[IonType, int], list[float]]
        :raises PeptacularError: If the ion type is unsupported, a charge is invalid,
            or the annotation contains unknown mods or interval mods.
        :raises UnsupportedOperationError: If a ``#XL`` or ``#BRANCH`` label is on two residues
            (a cross-linked or branched peptidoform).
        """
        return _frag_engine.fast_fragment(self, ion_types, charges, monoisotopic=monoisotopic)

    """
    Slicing Methods
    """

    def slice_by_span(self, span: Span | tuple[int, int, int], *, inplace: bool = False) -> Self:
        return self.slice(span[0], span[1], inplace=inplace)

    def slice(
        self,
        start: int | None,
        stop: int | None,
        *,
        inplace: bool = False,
        _keep_cut_localization_groups: bool = False,
    ) -> Self:
        """Return a sub-annotation spanning ``sequence[start:stop]``, carrying over applicable mods.

        :param start: 0-based start index (inclusive); ``None`` means the beginning.
        :type start: int | None
        :param stop: 0-based stop index (exclusive); ``None`` means the end.
        :type stop: int | None
        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :param _keep_cut_localization_groups: Private. Skip the check that refuses a slice
            cutting a localization group. Only for callers that take every residue of the
            annotation apart, such as :meth:`split` and the fragment engine, where a group
            is divided between the pieces rather than lost.
        :type _keep_cut_localization_groups: bool
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
                _keep_cut_localization_groups=_keep_cut_localization_groups,
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

    def shift(self, n: int, *, keep_nterm: int = 0, keep_cterm: int = 0, inplace: bool = False) -> Self:
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

    def shuffle(self, *, seed: Any = None, keep_nterm: int = 0, keep_cterm: int = 0, inplace: bool = False) -> Self:
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

    def reverse(self, *, keep_nterm: int = 0, keep_cterm: int = 0, inplace: bool = False) -> Self:
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

    def sort(self, *, inplace: bool = False, key: Callable[[str], Any] | None = None, reverse: bool = False) -> Self:
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

    def sliding_windows(self, window_size: int, *, reverse: bool = False) -> Generator[Self, None, None]:
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

    def condense_static_mods(self, *, inplace: bool = True) -> Self:
        return cast(Self, condense_static_mods(self, inplace=inplace))

    def condense_to_peptidoform(self, *, inplace: bool = True) -> Self:
        return cast(Self, condense_to_peptidoform(self, inplace=inplace))

    def count_residues(self, *, include_mods: bool = True) -> dict[str, int]:
        return count_residues(self, include_mods=include_mods)

    def percent_residues(self, *, include_mods: bool = True) -> dict[str, float]:
        return percent_residues(self, include_mods=include_mods)

    def is_subsequence(self, other: Self, *, ignore_mods: bool = False, ignore_intervals: bool = True) -> bool:
        return is_subsequence(self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals)

    def find_indices(self, other: Self, *, ignore_mods: bool = False, ignore_intervals: bool = True) -> list[int]:
        return find_indices(self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals)

    def condense_mods_to_intervals(self, *, inplace: bool = True) -> Self:
        return cast(Self, condense_mods_to_intervals(self, inplace=inplace))

    def coverage(self, subsequences: Iterable[Self], *, accumulate: bool = False, ignore_mods: bool = False, ignore_ambiguity: bool = False) -> list[int]:
        return coverage(
            annotation=self,
            subsequences=subsequences,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )

    def percent_coverage(self, subsequences: Iterable[Self], *, accumulate: bool = False, ignore_mods: bool = False, ignore_ambiguity: bool = False) -> float:
        return percent_coverage(
            annotation=self,
            subsequences=subsequences,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )

    def modification_coverage(self, subsequences: Iterable[Self], *, ignore_ambiguity: bool = False, accumulate: bool = False) -> dict[int, int]:
        return modification_coverage(
            annotation=self,
            subsequences=subsequences,
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
        unique_peptidoforms: bool = False,
    ) -> Generator[Self, None, None]:
        """
        Build all modifications from intervals and mass shifts.

        By default every positional isomer is yielded (``unique_peptidoforms=False``). Set
        ``unique_peptidoforms=True`` to collapse isomers that share a modification
        composition to a single representative.
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

    def add_static_mod_by_residue(self, residue: str | Iterable[str], mod: Any, *, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().add_static_mod_by_residue(residue, mod, inplace=True)

        residues = list(residue)

        mod_str, count = convert_single_mod_input(mod)

        if count != 1:
            raise PeptacularError("Fixed modifications added by residue must have a count of 1.")

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
        *,
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

    def condense_ambiguity_to_xnotation(self, *, inplace: bool = True) -> Self:
        """Condense ambiguous interval regions to X-notation placeholders.

        :param inplace: Modify this object when ``True``; return a new annotation when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation with intervals condensed to X notation.
        :rtype: Self
        """
        return cast(Self, condense_ambiguity_to_xnotation(self, inplace=inplace))

    def localization_isomers(self, *, max_isomers: int | None = DEFAULT_MAX_ISOMERS) -> list[Self]:
        """Expand every ambiguous modification position into its concrete placements.

        Expands ``#label`` groups, ranges (``PEP(ST)[Phospho]IDE``) and unknown-position mods
        (``[Phospho]?PEPTIDE``; these can go on any residue), one mod per residue. A group's label and the chosen
        residue's score stay on the placed mod (``S[Phospho#g1(0.8)]``). See
        :func:`peptacular.localization_isomers` for the full rules and ordering.

        :param max_isomers: Raise :class:`PeptacularError` if there would be more isomers than
            this. Defaults to 10,000; pass ``None`` for no limit.
        :type max_isomers: int | None
        :return: Deduplicated isomers in a fixed order. This annotation is not modified.
        :rtype: list[Self]
        """
        return cast(list[Self], localization_isomers(self, max_isomers=max_isomers))

    def candidate_sites(self, mod: Any, *, residues: str) -> list[tuple[int, Self]]:
        """Place ``mod`` on each unmodified residue whose letter is in ``residues``.

        :param mod: The modification to place (``"Phospho"``, a mass, ...).
        :type mod: Any
        :param residues: One-letter codes that can carry ``mod`` (required; there is no default site table).
        :type residues: str
        :return: ``(position, isomer)`` pairs in sequence order, position 0-based.
        :rtype: list[tuple[int, Self]]
        """
        return cast(list[tuple[int, Self]], candidate_sites(self, mod, residues=residues))

    @staticmethod
    def group_by_ambiguity(annotations: Iterable["ProFormaAnnotation"], *, precision: int = 5) -> list[tuple["ProFormaAnnotation", ...]]:
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
    def unique_fragments(annotations: Iterable["ProFormaAnnotation"], *, precision: int = 4) -> list[int]:
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
        """Convert this annotation to an IP2 sequence string.

        IP2 writes each modification in parentheses directly after the residue it sits on.
        N-terminal modifications are leading parenthesis groups, and C-terminal modifications
        are trailing groups after the last residue's own modification. The modification string
        is copied through unchanged, and no flanking residues are written. The output is the
        inverse of :meth:`from_ip2_sequence`, so reading it back gives an equal annotation.
        Static modifications are written onto the residues they apply to.

        >>> import peptacular as pt
        >>> pt.parse("PEP[Phospho]TIDE").to_ip2()
        'PEP(Phospho)TIDE'
        >>> pt.parse("PEP[+79.9663]TIDE").to_ip2()
        'PEP(+79.9663)TIDE'
        >>> pt.parse("[Acetyl]-PEPTIDE").to_ip2()
        '(Acetyl)PEPTIDE'
        >>> pt.parse("PEPTIDE[Phospho]-[Amidated]").to_ip2()
        'PEPTIDE(Phospho)(Amidated)'

        :return: The annotation as an IP2 sequence without flanking residues.
        :rtype: str
        :raises PeptacularError: If the annotation has isotope labels, labile modifications,
            unknown-position modifications, interval modifications, a charge or charge adducts,
            a C-terminal modification on a last residue with no modification of its own, more
            than one modification on the last residue, a modification multiplier, or a
            modification name containing parentheses.
        """
        for mod_type, feature in (
            (ModType.ISOTOPE, "isotope labels"),
            (ModType.LABILE, "labile modifications"),
            (ModType.UNKNOWN, "unknown-position modifications"),
            (ModType.INTERVAL, "interval modifications"),
            (ModType.CHARGE, "charge and charge adducts"),
        ):
            if self.has_mods(mod_type):
                raise PeptacularError(f"IP2 format does not support {feature}.")

        annot_copy = self.copy()
        annot_copy.condense_static_mods(inplace=True)

        def group(mods: Mapping[str, int]) -> str:
            out: list[str] = []
            for name, count in mods.items():
                if count != 1:
                    raise PeptacularError("IP2 format does not support modification multipliers.")
                if "(" in name or ")" in name:
                    raise PeptacularError(f"IP2 format cannot write a modification name containing parentheses: {name!r}.")
                out.append(f"({name})")
            return "".join(out)

        sequence = annot_copy.stripped_sequence
        last = len(sequence) - 1
        internal = annot_copy._internal_mods or {}
        cterm = annot_copy._cterm_mods or {}

        # The reader rewrites each parenthesis group into square brackets and then decides
        # which group is terminal by those brackets, so a terminal name holding one of them
        # comes back as something else. An internal name is read back unchanged.
        for mods, where in ((annot_copy._nterm_mods or {}, "N-terminus"), (cterm, "C-terminus")):
            for name in mods:
                if "[" in name or "]" in name:
                    raise PeptacularError(f"IP2 format cannot write a modification name containing square brackets at the {where}: {name!r}.")

        # The same brackets defeat the reader when they sit on the last residue and a
        # C-terminal group follows it: the trailing run is then not seen as two groups and
        # the C-terminal modification comes back on the residue.
        if cterm and last >= 0:
            for name in internal.get(last, {}):
                if "[" in name or "]" in name:
                    raise PeptacularError(
                        f"IP2 format cannot write a C-terminal modification after a last-residue modification name containing square brackets: {name!r}."
                    )

        # The reader gives the first group after the last residue to that residue and every
        # later group to the C-terminus, so the writer must emit exactly that shape.
        if last >= 0:
            if cterm and last not in internal:
                raise PeptacularError("IP2 format cannot place a C-terminal modification on a last residue that has no modification of its own.")
            if last in internal and len(internal[last]) > 1:
                raise PeptacularError("IP2 format cannot place multiple modifications on the last residue, they would be read as C-terminal modifications.")

        parts: list[str] = [group(annot_copy._nterm_mods or {})]
        for index, residue in enumerate(sequence):
            parts.append(residue)
            if index in internal:
                parts.append(group(internal[index]))
        parts.append(group(cterm))
        return "".join(parts)

    @staticmethod
    def from_ip2_sequence(sequence: str) -> "ProFormaAnnotation":
        """Internal function for converting a single IP2 sequence."""
        if re.match(r"^([A-Z]|-)\..*\.([A-Z]|-)$", sequence):
            sequence = sequence[2:-2]
        sequence = re.sub(r"\(([^)]+)\)", r"[\1]", sequence)

        # A run of adjacent brackets is handled differently depending on where it sits:
        #  - leading (N-terminal) run: a single '-' after the whole run, e.g.
        #    '[mod1][mod2]PEPTIDE' -> '[mod1][mod2]-PEPTIDE' (stacked N-term mods).
        #  - trailing run at the very end of the sequence: the first bracket stays
        #    attached to the preceding residue, a single '-' follows it, and any
        #    further brackets stack directly after (ProForma allows only one '-'
        #    C-terminal separator), e.g. 'PEPTIDE[1][2][3]' -> 'PEPTIDE[1]-[2][3]'
        #    (residue mod followed by stacked C-terminal mods).
        #  - a run sandwiched between residues on both sides is left untouched, e.g.
        #    'PEP[mod1][mod2]TIDE' stays as-is (multiple mods on one residue).
        # A single prior blanket substitution (every '][' -> ']-[') got the leading
        # and trailing cases right only by coincidence and produced unparseable
        # output for the sandwiched case.
        nterm_match = re.match(r"^(\[[^\]]+\])+", sequence)
        if nterm_match and nterm_match.end() < len(sequence):
            end = nterm_match.end()
            sequence = sequence[:end] + "-" + sequence[end:]

        def _dash_before_later_brackets(m: re.Match[str]) -> str:
            brackets = re.findall(r"\[[^\]]+\]", m.group(0))
            return brackets[0] + "-" + "".join(brackets[1:])

        sequence = re.sub(r"(?:\[[^\]]+\]){2,}$", _dash_before_later_brackets, sequence)

        return ProFormaAnnotation.parse(sequence)

    def to_diann(self) -> str:
        """Convert this annotation to a DIA-NN format string.

        The sequence is enclosed in underscores, with modifications written in square brackets
        immediately following the modified residue. N-terminal modifications are written inside
        the opening underscore, and C-terminal modifications are written after the closing underscore.
        The modification string is copied through unchanged. Static modifications are written onto
        the residues they apply to.

        Two modifications at one terminus are refused rather than written. ``from_diann`` gives
        the leading bracket group a dash and rewrites a single trailing group, so a second group
        at either terminus would come back as a string it cannot parse.

        >>> import peptacular as pt
        >>> pt.parse("PEPTIDE").to_diann()
        '_PEPTIDE_'
        >>> pt.parse("PEP[Phospho]TIDE").to_diann()
        '_PEP[Phospho]TIDE_'
        >>> pt.parse("PEP[+79.9663]TIDE").to_diann()
        '_PEP[+79.9663]TIDE_'
        >>> pt.parse("[Acetyl]-PEPTIDE").to_diann()
        '_[Acetyl]PEPTIDE_'
        >>> pt.parse("PEPTIDE-[Amidated]").to_diann()
        '_PEPTIDE_[Amidated]'
        >>> pt.parse("[Acetyl]-PEP[Phospho]TIDE-[Amidated]").to_diann()
        '_[Acetyl]PEP[Phospho]TIDE_[Amidated]'

        :return: The annotation as a DIA-NN sequence.
        :rtype: str
        :raises PeptacularError: If the annotation has isotope labels, labile modifications,
            unknown-position modifications, interval modifications, charge or charge adducts,
            modification multipliers, more than one N-terminal modification, or more than one
            C-terminal modification.
        """
        for mod_type, feature in (
            (ModType.ISOTOPE, "isotope labels"),
            (ModType.LABILE, "labile modifications"),
            (ModType.UNKNOWN, "unknown-position modifications"),
            (ModType.INTERVAL, "interval modifications"),
            (ModType.CHARGE, "charge or charge adducts"),
        ):
            if self.has_mods(mod_type):
                raise PeptacularError(f"DIA-NN format does not support {feature}.")

        annot_copy = self.copy()
        annot_copy.condense_static_mods(inplace=True)

        def group(mods: Mapping[str, int]) -> str:
            out: list[str] = []
            for name, count in mods.items():
                if count != 1:
                    raise PeptacularError("DIA-NN format does not support modification multipliers.")
                out.append(f"[{name}]")
            return "".join(out)

        sequence = annot_copy.stripped_sequence
        internal = annot_copy._internal_mods or {}
        nterm = annot_copy._nterm_mods or {}
        cterm = annot_copy._cterm_mods or {}

        # The reader matches a terminal group with a bracket pattern that stops at the first
        # `]`, so a terminal name holding a square bracket is cut in the wrong place and the
        # result does not parse. An internal name is read back unchanged.
        for mods, where in ((nterm, "N-terminus"), (cterm, "C-terminus")):
            for name in mods:
                if "[" in name or "]" in name:
                    raise PeptacularError(f"DIA-NN format cannot write a modification name containing square brackets at the {where}: {name!r}.")

        # `from_diann` gives the leading bracket group a dash and rewrites a single trailing
        # group, so a second group at either terminus comes back as a string it cannot parse.
        if len(nterm) > 1:
            raise PeptacularError("DIA-NN format cannot write more than one N-terminal modification, from_diann would not read the result back.")
        if len(cterm) > 1:
            raise PeptacularError("DIA-NN format cannot write more than one C-terminal modification, from_diann would not read the result back.")

        parts: list[str] = ["_", group(nterm)]
        for index, residue in enumerate(sequence):
            parts.append(residue)
            if index in internal:
                parts.append(group(internal[index]))
        parts.append("_")
        parts.append(group(cterm))
        return "".join(parts)

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

    def to_maxquant(self) -> str:
        """Convert this annotation to a MaxQuant modified-sequence string.

        MaxQuant writes the sequence inside underscores, with each modification in round
        brackets directly after the residue it sits on. The group before the first residue is
        the N-terminal modification, and a group after the last residue is read back as a
        C-terminal modification. The modification string is copied through unchanged, so both
        the abbreviated (``ox``) and the full-name (``Oxidation (M)``) spellings are written
        as given. Static modifications are written onto the residues they apply to.

        MaxQuant writes at most one modification per position, so a second modification on one
        residue or at one terminus is refused. A modification on the last residue is refused
        too, because ``from_maxquant`` reads a trailing group as a C-terminal modification and
        the two cannot be told apart.

        >>> import peptacular as pt
        >>> pt.parse("PEPTIDE").to_maxquant()
        '_PEPTIDE_'
        >>> pt.parse("PEPM[Oxidation (M)]TIDE").to_maxquant()
        '_PEPM(Oxidation (M))TIDE_'
        >>> pt.parse("[ac]-SM[ox]K").to_maxquant()
        '_(ac)SM(ox)K_'
        >>> pt.parse("PEPTIDE-[Amidated]").to_maxquant()
        '_PEPTIDE(Amidated)_'
        >>> pt.parse("PEP[Phospho]TIDE-[Amidated]").to_maxquant()
        '_PEP(Phospho)TIDE(Amidated)_'

        :return: The annotation as a MaxQuant modified sequence.
        :rtype: str
        :raises PeptacularError: If the annotation has isotope labels, labile modifications,
            unknown-position modifications, interval modifications, charge or charge adducts,
            modification multipliers, more than one modification at a terminus or on one
            residue, a modification on the last residue, or a modification name with
            unbalanced parentheses.
        """
        for mod_type, feature in (
            (ModType.ISOTOPE, "isotope labels"),
            (ModType.LABILE, "labile modifications"),
            (ModType.UNKNOWN, "unknown-position modifications"),
            (ModType.INTERVAL, "interval modifications"),
            (ModType.CHARGE, "charge or charge adducts"),
        ):
            if self.has_mods(mod_type):
                raise PeptacularError(f"MaxQuant format does not support {feature}.")

        annot_copy = self.copy()
        annot_copy.condense_static_mods(inplace=True)

        def group(mods: Mapping[str, int], where: str) -> str:
            if len(mods) > 1:
                raise PeptacularError(f"MaxQuant format cannot write more than one modification {where}.")
            out: list[str] = []
            for name, count in mods.items():
                if count != 1:
                    raise PeptacularError("MaxQuant format does not support modification multipliers.")
                depth = 0
                for char in name:
                    depth += (char == "(") - (char == ")")
                    if depth < 0:
                        break
                if depth != 0:
                    raise PeptacularError(f"MaxQuant format cannot write a modification name with unbalanced parentheses: {name!r}.")
                out.append(f"({name})")
            return "".join(out)

        sequence = annot_copy.stripped_sequence
        last = len(sequence) - 1
        internal = annot_copy._internal_mods or {}

        # `from_maxquant` reads a group after the last residue as a C-terminal modification.
        if last in internal:
            raise PeptacularError("MaxQuant format cannot write a modification on the last residue, from_maxquant would read it as a C-terminal modification.")

        parts: list[str] = ["_", group(annot_copy._nterm_mods or {}, "at the N-terminus")]
        for index, residue in enumerate(sequence):
            parts.append(residue)
            if index in internal:
                parts.append(group(internal[index], "on one residue"))
        parts.append(group(annot_copy._cterm_mods or {}, "at the C-terminus"))
        parts.append("_")
        return "".join(parts)

    @staticmethod
    def from_maxquant(sequence: str) -> "ProFormaAnnotation":
        """Read a MaxQuant modified-sequence string into an annotation.

        MaxQuant writes the sequence inside underscores, with each modification in round
        brackets directly after the residue it sits on. A group before the first residue is
        the N-terminal modification, a group after the last residue is the C-terminal
        modification, and every other group sits on the residue before it. Nested
        parentheses are matched, so a full name such as ``Oxidation (M)`` stays in one
        piece, and the modification string is read through unchanged.

        >>> import peptacular as pt
        >>> pt.ProFormaAnnotation.from_maxquant("_PEPTIDE_").serialize()
        'PEPTIDE'
        >>> pt.ProFormaAnnotation.from_maxquant("_(ac)SM(ox)K_").serialize()
        '[ac]-SM[ox]K'
        >>> pt.ProFormaAnnotation.from_maxquant("_ALTDM(Oxidation (M))PQM(Oxidation (M))R_").serialize()
        'ALTDM[Oxidation (M)]PQM[Oxidation (M)]R'
        >>> pt.ProFormaAnnotation.from_maxquant("_P(Phospho)EPTIDE_").serialize()
        'P[Phospho]EPTIDE'
        >>> pt.ProFormaAnnotation.from_maxquant("_PEPTIDE(Amidated)_").serialize()
        'PEPTIDE-[Amidated]'

        :param sequence: A MaxQuant modified sequence, with or without the enclosing underscores.
        :type sequence: str
        :return: The sequence as an annotation.
        :rtype: ProFormaAnnotation
        :raises ProFormaFormatError: If the parentheses are unbalanced, or if the ProForma
            string the groups are rewritten into cannot be parsed.
        """
        sequence = sequence.strip("_")

        # Split into residues and round-bracket groups, matching nested parentheses so that
        # names such as "Oxidation (M)" stay in one piece.
        tokens: list[tuple[bool, str]] = []
        depth = 0
        start = 0
        for index, char in enumerate(sequence):
            if char == "(":
                if depth == 0:
                    start = index + 1
                depth += 1
            elif char == ")":
                if depth == 0:
                    raise ProFormaFormatError(f"Unbalanced parentheses in MaxQuant sequence: {sequence!r}.")
                depth -= 1
                if depth == 0:
                    tokens.append((True, sequence[start:index]))
            elif depth == 0:
                tokens.append((False, char))
        if depth != 0:
            raise ProFormaFormatError(f"Unbalanced parentheses in MaxQuant sequence: {sequence!r}.")

        # A leading group is N-terminal, a trailing group is C-terminal, the rest sit on the
        # residue before them.
        first = next((i for i, (is_mod, _) in enumerate(tokens) if not is_mod), len(tokens))
        last = next((i for i in range(len(tokens) - 1, -1, -1) if not tokens[i][0]), -1)
        out: list[str] = []
        if first:
            out.extend(f"[{name}]" for _, name in tokens[:first])
            out.append("-")
        for is_mod, text in tokens[first : last + 1]:
            out.append(f"[{text}]" if is_mod else text)
        if last >= 0 and last + 1 < len(tokens):
            out.append("-")
            out.extend(f"[{name}]" for _, name in tokens[last + 1 :])
        return ProFormaAnnotation.parse("".join(out))

    def to_skyline(self) -> str:
        """Convert this annotation to a Skyline modified-sequence string.

        Skyline writes each modification in square brackets directly after the residue it
        sits on, either as a signed mass difference (``[+57.0]``) or as a modification name
        (``[Carbamidomethyl]``). The modification string is copied through unchanged, so a
        delta mass keeps the precision it was given. Static modifications are written onto
        the residues they apply to.

        Skyline writes a terminal modification onto the first or last residue, and its
        mass-delta forms sum two modifications on one residue into a single bracket. Its name
        and UniMod forms give each modification its own bracket instead, ``[A][B]`` and
        ``(unimod:4)(unimod:35)``, which ``from_skyline`` reads. ``from_skyline`` cannot tell a
        terminal modification, or a summed mass, from an ordinary modification on that residue,
        so a terminal modification and a second modification on one residue are refused, as is
        a bare unsigned number, which ``from_skyline`` would read back as a signed delta mass.

        >>> import peptacular as pt
        >>> pt.parse("PEPTIDE").to_skyline()
        'PEPTIDE'
        >>> pt.parse("PEPT[+79.9663]IDE").to_skyline()
        'PEPT[+79.9663]IDE'
        >>> pt.parse("PEPT[-18.0106]IDE").to_skyline()
        'PEPT[-18.0106]IDE'
        >>> pt.parse("PEPC[Carbamidomethyl]IDE").to_skyline()
        'PEPC[Carbamidomethyl]IDE'

        :return: The annotation as a Skyline modified sequence.
        :rtype: str
        :raises PeptacularError: If the annotation has isotope labels, labile modifications,
            unknown-position modifications, interval modifications, charge or charge adducts,
            modification multipliers, a terminal modification, more than one modification on
            one residue, an unsigned bare number as a modification, or a modification name with
            unbalanced square brackets.
        """
        for mod_type, feature in (
            (ModType.ISOTOPE, "isotope labels"),
            (ModType.LABILE, "labile modifications"),
            (ModType.UNKNOWN, "unknown-position modifications"),
            (ModType.INTERVAL, "interval modifications"),
            (ModType.CHARGE, "charge or charge adducts"),
        ):
            if self.has_mods(mod_type):
                raise PeptacularError(f"Skyline format does not support {feature}.")

        annot_copy = self.copy()
        annot_copy.condense_static_mods(inplace=True)

        # Skyline writes a terminal modification onto the first or last residue, where
        # `from_skyline` cannot tell it from a modification on that residue.
        if annot_copy._nterm_mods or annot_copy._cterm_mods:
            raise PeptacularError(
                "Skyline format cannot write a terminal modification, Skyline writes it onto the "
                "terminal residue and from_skyline could not tell it from a modification on that residue."
            )

        sequence = annot_copy.stripped_sequence
        internal = annot_copy._internal_mods or {}

        parts: list[str] = []
        for index, residue in enumerate(sequence):
            parts.append(residue)
            mods = internal.get(index)
            if not mods:
                continue
            for count in mods.values():
                if count != 1:
                    raise PeptacularError("Skyline format does not support modification multipliers.")
            if len(mods) > 1:
                raise PeptacularError(
                    "Skyline format cannot write more than one modification on one residue: Skyline's mass-delta forms sum them into a single bracket, "
                    "which from_skyline cannot tell from one modification of that mass."
                )
            (name,) = mods
            depth = 0
            for char in name:
                depth += (char == "[") - (char == "]")
                if depth < 0:
                    break
            if depth != 0:
                raise PeptacularError(f"Skyline format cannot write a modification name with unbalanced square brackets: {name!r}.")
            if _SKYLINE_BARE_NUMBER.fullmatch(name) and name[0] not in "+-":
                raise PeptacularError(
                    f"Skyline format cannot write the unsigned number {name!r} as a modification, from_skyline would read it back as a signed delta mass."
                )
            parts.append(f"[{name}]")
        return "".join(parts)

    @staticmethod
    def from_skyline(sequence: str) -> "ProFormaAnnotation":
        """Read a Skyline modified-sequence string into an annotation.

        Skyline writes each modification directly after the residue it sits on, in square
        brackets, as a signed mass difference (``[+57]``, ``[-0.98]``) or as a modification
        name (``[Phospho (ST)]``), or as ``(unimod:N)`` in round brackets. Square brackets are
        matched as balanced pairs. A round-bracket ``(unimod:N)`` is read without regard to case
        as ``[UNIMOD:N]``, and a bracket holding a bare number gets the ``+`` ProForma requires
        when it has no sign. Skyline writes a terminal modification onto the first or last
        residue, so a Skyline string carries no terminal modification of its own. ProForma
        syntax that Skyline never writes is passed through to the parser rather than refused,
        as in the other readers, so ``from_skyline("PEPTIDE-[Amidated]")`` does return a
        C-terminal modification.

        >>> import peptacular as pt
        >>> pt.ProFormaAnnotation.from_skyline("GVDC[+57]QEVSQEK").serialize()
        'GVDC[+57]QEVSQEK'
        >>> pt.ProFormaAnnotation.from_skyline("GSMAFVPK(unimod:259)").serialize()
        'GSMAFVPK[UNIMOD:259]'
        >>> pt.ProFormaAnnotation.from_skyline("PEPT[Phospho (ST)]IDER").serialize()
        'PEPT[Phospho (ST)]IDER'
        >>> pt.ProFormaAnnotation.from_skyline("AAM[15.9949]Q[-0.98]R").serialize()
        'AAM[+15.9949]Q[-0.98]R'

        :param sequence: A Skyline modified sequence.
        :type sequence: str
        :return: The sequence as an annotation.
        :rtype: ProFormaAnnotation
        :raises ProFormaFormatError: If the brackets are unbalanced, a bracket is empty or comes
            before the first residue, a round bracket is not ``(unimod:N)``, or the ProForma
            string the brackets are rewritten into cannot be parsed.
        """
        out: list[str] = []
        index = 0
        while index < len(sequence):
            char = sequence[index]
            if char == "[":
                depth = 0
                for end in range(index, len(sequence)):
                    depth += (sequence[end] == "[") - (sequence[end] == "]")
                    if depth == 0:
                        break
                else:
                    raise ProFormaFormatError(f"Unbalanced square brackets in Skyline sequence: {sequence!r}.")
                # Stripped before the checks below: a bracket written with padding held a
                # bare number that the unstripped test read as a name, and `mass()` then
                # raised for a modification Skyline means as a delta.
                text = sequence[index + 1 : end].strip()
                if not text:
                    raise ProFormaFormatError(f"Empty modification in Skyline sequence: {sequence!r}.")
                if _SKYLINE_COMMA_NUMBER.fullmatch(text):
                    # Skyline's fallback writes the mass in the current culture
                    # (ModifiedSequence.cs:308) and its own parser accepts a comma decimal
                    # (MassModification.cs:108), so this is a mass written on a machine whose
                    # culture uses a comma. Which number it is depends on that culture, so it
                    # is refused rather than read as a name or guessed at.
                    raise ProFormaFormatError(
                        f"Comma decimal in a Skyline modification: {text!r} in {sequence!r}. Skyline writes a mass in the writer's "
                        f"culture, so the value is ambiguous here; rewrite it with a period."
                    )
                if _SKYLINE_BARE_NUMBER.fullmatch(text) and text[0] not in "+-":
                    text = f"+{text}"
                _skyline_require_residue(out, sequence)
                out.append(f"[{text}]")
                index = end + 1
            elif char == "(":
                end = sequence.find(")", index)
                match = _SKYLINE_UNIMOD.fullmatch(sequence[index + 1 : end]) if end != -1 else None
                if match is None:
                    raise ProFormaFormatError(f"Expected (unimod:N) for a round bracket in Skyline sequence: {sequence!r}.")
                _skyline_require_residue(out, sequence)
                out.append(f"[UNIMOD:{match.group(1)}]")
                index = end + 1
            elif char in "])":
                raise ProFormaFormatError(f"Unbalanced brackets in Skyline sequence: {sequence!r}.")
            else:
                out.append(char)
                index += 1
        return ProFormaAnnotation.parse("".join(out))

    def to_casanovo(self) -> str:
        """Convert this annotation to a Casanovo format string.

        Casanovo writes a modification as a signed delta mass with no brackets, directly after
        the residue it sits on, and an N-terminal modification as a signed delta mass at the
        very start of the sequence. Static modifications are condensed onto the residues they
        apply to.

        >>> import peptacular as pt
        >>> pt.parse("PEPTIDE").to_casanovo()
        'PEPTIDE'
        >>> pt.parse("PEPM[+15.995]IDE").to_casanovo()
        'PEPM+15.995IDE'
        >>> pt.parse("[+43.006]-PEPTIDE").to_casanovo()
        '+43.006PEPTIDE'
        >>> pt.parse("[+42.011]-PEPM[+15.995]IDEC[+57.021]").to_casanovo()
        '+42.011PEPM+15.995IDEC+57.021'

        :return: The annotation as a Casanovo sequence string.
        :rtype: str
        :raises PeptacularError: If the annotation has modifications that are not signed delta masses,
            multiple modifications on one residue, multiple N-terminal modifications, C-terminal modifications,
            isotope labels, labile modifications, unknown-position modifications, interval modifications,
            charge or charge adducts, or modification multipliers.
        """
        for mod_type, feature in (
            (ModType.ISOTOPE, "isotope labels"),
            (ModType.LABILE, "labile modifications"),
            (ModType.UNKNOWN, "unknown-position modifications"),
            (ModType.INTERVAL, "interval modifications"),
            (ModType.CHARGE, "charge or charge adducts"),
        ):
            if self.has_mods(mod_type):
                raise PeptacularError(f"Casanovo format does not support {feature}.")

        annot_copy = self.copy()
        annot_copy.condense_static_mods(inplace=True)

        if annot_copy._cterm_mods:
            raise PeptacularError("Casanovo format does not support C-terminal modifications.")

        nterm = annot_copy._nterm_mods or {}
        if len(nterm) > 1:
            raise PeptacularError("Casanovo format does not support multiple N-terminal modifications.")

        nterm_str = ""
        for name, count in nterm.items():
            if count != 1:
                raise PeptacularError("Casanovo format does not support modification multipliers.")
            if not _SIGNED_DELTA_RE.match(name):
                raise PeptacularError(f"Casanovo format only supports signed delta mass modifications, not {name!r}.")
            nterm_str = name

        sequence = annot_copy.stripped_sequence
        internal = annot_copy._internal_mods or {}
        parts: list[str] = [nterm_str]
        for index, residue in enumerate(sequence):
            parts.append(residue)
            if index in internal:
                mods = internal[index]
                if len(mods) > 1:
                    raise PeptacularError("Casanovo format does not support multiple modifications on one residue.")
                for name, count in mods.items():
                    if count != 1:
                        raise PeptacularError("Casanovo format does not support modification multipliers.")
                    if not _SIGNED_DELTA_RE.match(name):
                        raise PeptacularError(f"Casanovo format only supports signed delta mass modifications, not {name!r}.")
                    parts.append(name)
        return "".join(parts)

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

    def to_ms2_pip(self, *, inplace: bool = False) -> tuple[str, str]:
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
            raise PeptacularError("MS2PIP format does not support isotope, labile, unknown, interval, charge, or charge adduct modifications.")

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
                    raise PeptacularError("MS2PIP format does not support modification multipliers.")
                mod_tuples.append((0, mod_name))

        # Process C-terminal modifications
        if annot_copy._cterm_mods is not None:
            for mod_name, count in annot_copy._cterm_mods.items():
                if count != 1:
                    raise PeptacularError("MS2PIP format does not support modification multipliers.")
                mod_tuples.append((-1, mod_name))

        # Process internal modifications
        if annot_copy._internal_mods is not None:
            for index, mods_dict in annot_copy._internal_mods.items():
                if len(mods_dict) > 1:
                    raise PeptacularError("MS2PIP format does not support multiple modifications at the same site.")
                for mod_name, count in mods_dict.items():
                    if count != 1:
                        raise PeptacularError("MS2PIP format does not support modification multipliers.")
                    # MS2PIP uses 1-indexed positions
                    mod_tuples.append((index + 1, mod_name))

        unmod_sequence = annot_copy.stripped_sequence

        mod_str = "|".join(f"{loc}|{name}" for loc, name in mod_tuples)

        return unmod_sequence, mod_str

    @staticmethod
    def from_ms2_pip(sequence: str, mod_str: str, *, static_mods: Mapping[str, float | int | str] | None = None) -> "ProFormaAnnotation":
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
            raise PeptacularError(f"Invalid MS2PIP modification string format: {mod_str}")

        # Process modifications in pairs (location, name)
        for i in range(0, len(mod_parts), 2):
            loc_str = mod_parts[i]
            mod_name = mod_parts[i + 1]

            # Parse location
            try:
                loc = int(loc_str)
            except ValueError:
                raise PeptacularError(f"Invalid MS2PIP modification location {loc_str!r} in {mod_str!r}") from None

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
                    raise PeptacularError(f"Modification location {loc} is out of range for sequence of length {len(sequence)}")
                annot.append_internal_mod_at_index(internal_index, mod_name, inplace=True)

        # Add static modifications
        for static_aa, mass in (static_mods or {}).items():
            annot.add_static_mod_by_residue(static_aa, mass, inplace=True)

        return annot

    def to_msfragger(self, *, inplace: bool = False) -> tuple[str, str]:
        """Convert peptide annotation to MSFragger (Peptide, Assigned Modifications) format.

        FragPipe's ``psm.tsv`` carries the unmodified sequence in ``Peptide`` and the
        modifications assigned to it in ``Assigned Modifications`` (e.g.
        ``4M(15.9949), N-term(42.0106)``). Fixed modifications are in that column too, not only
        variable ones: Philosopher filters on ``i.Type == mod.Assigned`` alone
        (``lib/rep/rep.go:591``) and its own test expects a ``Variable: false`` entry in the list
        (``lib/id/pep_test.go:116``). Static modifications are condensed onto their residues and
        written, which is what that column holds.

        Each mass is written with four decimals and the entries are sorted as strings, as
        Philosopher writes them (``%.4f`` at ``lib/rep/rep.go:596``, ``sort.Strings`` at
        ``lib/rep/psm.go:564``), so ``10P(15.9949), 4P(15.9949)`` is the order for positions 10
        and 4, and a residue entry precedes ``N-term``.

        FragPipe's ``psm.tsv`` also contains a ``Modified Peptide`` column with residue masses
        (e.g. ``M[147]``). This method intentionally writes ``(Peptide, Assigned Modifications)``
        pairs rather than ``Modified Peptide``, because ``Modified Peptide`` records integer residue
        masses from which the original delta mass (such as ``+15.9949``) cannot be recovered or read back.

        :param inplace: Modify this object when ``True``; operate on a copy when ``False``.
        :type inplace: bool
        :return: A tuple of ``(unmodified_sequence, assigned_modifications)``.
        :rtype: tuple[str, str]
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
            raise PeptacularError("MSFragger format does not support isotope, labile, unknown, interval, charge, or charge adduct modifications.")

        if not inplace:
            annot_copy = self.copy()
            annot_copy.condense_static_mods(inplace=True)
        else:
            self.condense_static_mods(inplace=True)
            annot_copy = self

        def _format_delta_mass(name: str) -> str:
            if not _SIGNED_DELTA_RE.match(name):
                raise PeptacularError(f"MSFragger format only supports signed delta mass modifications, not {name!r}.")
            # Philosopher writes every mass with four decimals:
            # `fmt.Sprintf("%s%s(%.4f)", pos, i.AminoAcid, i.MassDiff)` (lib/rep/rep.go:596).
            # Copying the ProForma text gave "4M(16)" for `[+16]` and "4M(.5)" for `[+.5]`,
            # neither of which FragPipe writes.
            return f"{float(name):.4f}"

        mod_entries: list[str] = []

        # Process N-terminal modifications
        if annot_copy._nterm_mods is not None:
            for mod_name, count in annot_copy._nterm_mods.items():
                if count != 1:
                    raise PeptacularError("MSFragger format does not support modification multipliers.")
                mass_str = _format_delta_mass(mod_name)
                mod_entries.append(f"N-term({mass_str})")

        # Process internal modifications
        if annot_copy._internal_mods is not None:
            for index in sorted(annot_copy._internal_mods.keys()):
                mods_dict = annot_copy._internal_mods[index]
                for mod_name, count in mods_dict.items():
                    if count != 1:
                        raise PeptacularError("MSFragger format does not support modification multipliers.")
                    mass_str = _format_delta_mass(mod_name)
                    pos = index + 1
                    residue = annot_copy.sequence[index]
                    mod_entries.append(f"{pos}{residue}({mass_str})")

        # Process C-terminal modifications
        if annot_copy._cterm_mods is not None:
            for mod_name, count in annot_copy._cterm_mods.items():
                if count != 1:
                    raise PeptacularError("MSFragger format does not support modification multipliers.")
                mass_str = _format_delta_mass(mod_name)
                mod_entries.append(f"C-term({mass_str})")

        # Philosopher sorts the entry strings and joins them with ", ": `sort.Strings(assL)`
        # (lib/rep/psm.go:564) then `strings.Join(assL, ", ")` (lib/rep/psm.go:640). It is a
        # plain string sort, so "10P(...)" comes before "4P(...)" and "4M(...)" before
        # "N-term(...)".
        unmod_sequence = annot_copy.stripped_sequence
        return unmod_sequence, ", ".join(sorted(mod_entries))

    @staticmethod
    def from_msfragger(peptide: str, modifications: str) -> "ProFormaAnnotation":
        """Create ProFormaAnnotation from MSFragger peptide and assigned modifications pair.

        :param peptide: The unmodified amino acid sequence.
        :type peptide: str
        :param modifications: The comma-separated assigned modifications string
            (e.g. ``"N-term(42.0106), 4M(15.9949)"``).
        :type modifications: str
        :return: A parsed ProFormaAnnotation object.
        :rtype: ProFormaAnnotation
        """
        # `Peptide` is the stripped sequence: the FragPipe output guide calls it the
        # "peptide amino acid sequence, any modifications not included ('stripped' peptide
        # sequence)". Anything else here is not that column, and reading it anyway doubled a
        # modification already written in ProForma, or carried a charge through unchecked.
        if not isinstance(peptide, str):
            raise PeptacularError(f"MSFragger Peptide must be a string, not {type(peptide).__name__}.")
        if not _MSFRAGGER_PEPTIDE_RE.match(peptide):
            raise PeptacularError(
                f"MSFragger Peptide must be a bare upper-case residue sequence, not {peptide!r}. "
                f"That column holds the stripped sequence, so it carries no modification, charge or separator."
            )

        annot = ProFormaAnnotation(sequence=peptide)

        # An empty `Assigned Modifications` cell means no modifications: in
        # compomics/psm_utils `example_files/fragpipe.psm.tsv` the one row with an empty cell
        # (FLLEAGADQEHK) also has a blank `Modified Peptide`, which the FragPipe output guide
        # defines as "blank if peptide is unmodified". A reader handing the cell over as None,
        # or as the NaN pandas gives an empty cell, means that same empty cell.
        if modifications is None or (isinstance(modifications, float) and modifications != modifications):
            return annot
        if not isinstance(modifications, str):
            raise PeptacularError(f"MSFragger Assigned Modifications must be a string, or None or NaN for an empty cell, not {type(modifications).__name__}.")

        if modifications.strip() == "":
            return annot

        parts = [part.strip() for part in modifications.split(",")]
        for part in parts:
            if not part:
                raise PeptacularError(f"Invalid MSFragger modification entry: {part!r}")

            match = _MSFRAGGER_ENTRY_RE.match(part)
            if not match:
                raise PeptacularError(f"Invalid MSFragger modification entry: {part!r}")

            nterm, cterm, pos_str, residue, mass_str = match.groups()
            mass_val = float(mass_str)

            if nterm:
                annot.append_nterm_mod(mass_val, inplace=True)
            elif cterm:
                annot.append_cterm_mod(mass_val, inplace=True)
            else:
                pos = int(pos_str)
                if pos < 1 or pos > len(peptide):
                    raise PeptacularError(f"Modification position {pos} is out of range for peptide of length {len(peptide)}")
                expected_residue = peptide[pos - 1]
                if residue != expected_residue:
                    raise PeptacularError(
                        f"Residue mismatch at position {pos}: expected {residue!r} from modification, but found {expected_residue!r} in peptide {peptide!r}"
                    )
                annot.append_internal_mod_at_index(pos - 1, mass_val, inplace=True)

        return annot

    def isotopic_distribution(
        self,
        charge: CHARGE_TYPE | None = None,
        *,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        max_isotopes: int | None = None,
        min_abundance_threshold: float = 0.001,
    ) -> list[IsotopicData]:
        """Calculate the aggregated isotopic distribution from elemental composition.

        :param ion_type: Fragment ion type to use.
        :type ion_type: ION_TYPE
        :param charge: Charge override; uses the annotation charge when ``None``.
        :type charge: CHARGE_TYPE | None
        :param isotopes: Isotope offsets or element-count overrides.
        :type isotopes: ISOTOPE_TYPE | None
        :param deltas: Custom neutral-loss or gain formula(s).
        :type deltas: CUSTOM_LOSS_TYPE | None
        :param max_isotopes: Maximum nominal isotope window, or ``None`` for adaptive sizing.
        :type max_isotopes: int | None
        :param min_abundance_threshold: Minimum relative abundance (vs. the most abundant peak).
        :type min_abundance_threshold: float
        :return: Aggregated isotope peaks sorted by neutron offset.
        :rtype: list[IsotopicData]
        """
        frag_annot = self
        if charge is not None:  # update charge
            frag_annot = frag_annot.set_charge(charge, inplace=False)

        fragment = frag_annot.frag(ion_type=ion_type, isotopes=isotopes, deltas=deltas, calculate_with_composition=True)
        composition = fragment.composition
        assert composition is not None

        peaks = brain_isotopic_distribution(
            formula=cast(Mapping[str | ElementInfo, int | float], composition),
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
            charge=fragment.charge_state,
        )
        # The composition counts an H atom per proton; lift each to PROTON_MASS, as frag() does.
        binding = proton_binding_offset(frag_annot.charge_adducts, True)
        if binding:
            peaks = [IsotopicData(mass=peak.mass + binding, neutron_count=peak.neutron_count, abundance=peak.abundance) for peak in peaks]
        return peaks

    def estimate_isotopic_distribution(
        self,
        charge: CHARGE_TYPE | None = None,
        *,
        ion_type: ION_TYPE = IonType.PRECURSOR,
        isotopes: ISOTOPE_TYPE | None = None,
        deltas: CUSTOM_LOSS_TYPE | None = None,
        max_isotopes: int | None = None,
        min_abundance_threshold: float = 0.001,
    ) -> list[IsotopicData]:
        """Estimate an aggregated isotopic distribution based on mass."""

        mass = self.mass(ion_type=ion_type, charge=charge, isotopes=isotopes, deltas=deltas)

        return estimate_isotopic_distribution(
            neutral_mass=mass,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
        )

    @staticmethod
    def random(
        *,
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

    def left_semi_spans(self, *, min_len: int | None = None, max_len: int | None = None) -> Generator[Span, None, None]:
        """Get left semi-enzymatic sequences (N-terminus fixed)."""
        return left_semi_spans(self, min_len=min_len, max_len=max_len)

    def right_semi_spans(self, *, min_len: int | None = None, max_len: int | None = None) -> Generator[Span, None, None]:
        """Get right semi-enzymatic sequences (C-terminus fixed)."""
        return right_semi_spans(self, min_len=min_len, max_len=max_len)

    def semi_spans(self, *, min_len: int | None = None, max_len: int | None = None) -> Generator[Span, None, None]:
        """Get all semi-enzymatic sequences."""
        return semi_spans(self, min_len=min_len, max_len=max_len)

    def nonspecific_spans(self, *, min_len: int | None = None, max_len: int | None = None) -> Generator[Span, None, None]:
        """Get all non-enzymatic sequences (all possible subsequences)."""
        return nonspecific_spans(self, min_len=min_len, max_len=max_len)

    def cleavage_sites(
        self,
        enzyme: str | re.Pattern[str],
    ) -> Generator[int, None, None]:
        """Yield 0-based cleavage positions for ``enzyme``.

        :param enzyme: A protease name from tacular's ``PROTEASE_LOOKUP`` or a compiled pattern.
            A plain string is never treated as a regex.
        :type enzyme: str | re.Pattern[str]
        :raises UnknownEnzymeError: If ``enzyme`` is a string that names no known protease.
        :return: Generator of 0-based indices where cleavage occurs.
        :rtype: Generator[int, None, None]
        """
        # Call the underlying function
        return get_cleavage_sites(self, enzyme)

    def simple_cleavage_sites(
        self, cleave_on: str, *, restrict_before: str = "", restrict_after: str = "", cterminal: bool = True
    ) -> Generator[int, None, None]:
        """Get cleavage sites using simple amino acid rules."""
        pattern = generate_regex(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        )
        return self.cleavage_sites(pattern)

    def digest_spans(
        self, enzyme: str | re.Pattern[str], *, missed_cleavages: int = 0, semi: bool = False, min_len: int | None = None, max_len: int | None = None
    ) -> Generator[Span, None, None]:
        """Digest this annotation and yield the :class:`Span` of each peptide.

        :param enzyme: A protease name from tacular's ``PROTEASE_LOOKUP`` or a compiled pattern.
        :raises UnknownEnzymeError: If ``enzyme`` is a string that names no known protease.

        Use ``annotation[span]`` to get a peptide, or :func:`peptacular.digest` for
        ``(peptide, span)`` pairs.
        """
        return digest_annotation_by_regex(
            annotation=self,
            enzyme=enzyme,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )

    def simple_digest_spans(
        self,
        cleave_on: str,
        *,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Digest this annotation with amino-acid cleavage rules and yield the :class:`Span` of each peptide."""
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

    def sequential_digest_spans(
        self, enzyme_configs: list[EnzymeConfig], *, min_len: int | None = None, max_len: int | None = None
    ) -> Generator[Span, None, None]:
        """Digest with each :class:`EnzymeConfig` in turn and yield the :class:`Span` of each final peptide."""
        return sequential_digest_annotation(self, enzyme_configs, min_len=min_len, max_len=max_len)

    @property
    def prop(self) -> AnnotationProperties:
        """Get the properties of this annotation."""
        return AnnotationProperties(self.stripped_sequence)


def _adopt_mixin_qualnames() -> None:
    """Name the inherited mod accessors ``ProFormaAnnotation.<name>``, not ``_ModAccessMixin.<name>``.

    Keeps the private base class out of reprs, argument errors and tracebacks. Only
    ``ProFormaAnnotation`` inherits ``_ModAccessMixin``, so the rename is exact.
    """
    prefix = f"{_ModAccessMixin.__qualname__}."
    for member in vars(_ModAccessMixin).values():
        funcs = (member.fget, member.fset, member.fdel) if isinstance(member, property) else (member,)
        for func in funcs:
            if isinstance(func, FunctionType) and func.__qualname__.startswith(prefix):
                func.__qualname__ = f"{ProFormaAnnotation.__qualname__}.{func.__qualname__[len(prefix) :]}"


_adopt_mixin_qualnames()
del _adopt_mixin_qualnames

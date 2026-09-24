"""Modification accessors of :class:`~peptacular.annotation.ProFormaAnnotation`.

The ``validate_*``, ``has_*``, ``get_*``, ``set_*``, ``append_*``, ``extend_*``,
``remove_*``, ``pop_*`` and ``clear_*`` methods live here as :class:`_ModAccessMixin`,
which ``ProFormaAnnotation`` inherits. They are documented and used as methods of
``ProFormaAnnotation``; this module is private.
"""

from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from enum import StrEnum
from typing import TYPE_CHECKING, Any, Self

from tacular import AA_LOOKUP

from ..constants import ModType, ModTypeLiteral
from ..diagnostics import InvalidPositionError, PeptacularError
from ..proforma_components import (
    MODIFICATION_TYPE,
    FixedModification,
    GlobalChargeCarrier,
    IsotopeReplacement,
    ModificationTags,
    PositionScore,
)
from ..utils import _resolve_mod_types
from .mod import (
    Interval,
    Mod,
    Mods,
    as_mod_iterable,
    convert_moddict_input,
    convert_single_mod_input,
    is_mod_collection,
)

EMPTY_ISOTOPE_MODS = Mods[IsotopeReplacement](mod_type=ModType.ISOTOPE, _mods=None)
EMPTY_STATIC_MODS = Mods[FixedModification](mod_type=ModType.STATIC, _mods=None)
EMPTY_UNKNOWN_MODS = Mods[ModificationTags](mod_type=ModType.UNKNOWN, _mods=None)
EMPTY_LABILE_MODS = Mods[ModificationTags](mod_type=ModType.LABILE, _mods=None)
EMPTY_NTERM_MODS = Mods[ModificationTags](mod_type=ModType.NTERM, _mods=None)
EMPTY_CTERM_MODS = Mods[ModificationTags](mod_type=ModType.CTERM, _mods=None)
EMPTY_INTERNAL_MODS = Mods[ModificationTags](mod_type=ModType.INTERNAL, _mods=None)


def _concrete_position_labels(mods: "Mods | None") -> Iterable[str]:
    """Yield the ambiguous-position label (``#label``) id of every *concrete* modification
    in ``mods`` (bare ``[#label]`` references and ``PositionScore`` tags carry no concrete
    modification and are skipped). A label repeated in the output has multiple concrete
    occurrences, which ``validate_ambiguous_labels`` treats as an error."""
    if mods is None:
        return
    for mod in mods.mods:
        value = mod.value
        if not isinstance(value, ModificationTags):
            continue
        for tag in value.tags:
            position_id = getattr(tag, "position_id", None)
            if position_id is None or isinstance(tag, PositionScore):
                continue
            yield position_id


class ChargeType(StrEnum):
    INT = "int"
    ADDUCTS = "adducts"
    NONE = "none"


class _ModAccessMixin:
    """Modification accessors of :class:`~peptacular.annotation.ProFormaAnnotation` (private base class)."""

    if TYPE_CHECKING:
        # State and members that ProFormaAnnotation provides; declared for type checkers only.
        _sequence: str | None
        _compound_name: str | None
        _ion_name: str | None
        _peptide_name: str | None
        _isotope_mods: dict[str, int] | None
        _static_mods: dict[str, int] | None
        _labile_mods: dict[str, int] | None
        _unknown_mods: dict[str, int] | None
        _nterm_mods: dict[str, int] | None
        _cterm_mods: dict[str, int] | None
        _internal_mods: dict[int, dict[str, int]] | None
        _intervals: list[Interval] | None
        _charge: int | list[str] | None
        _validate: bool

        @property
        def sequence(self) -> str: ...
        @sequence.setter
        def sequence(self, value: str | None) -> None: ...
        @property
        def start_aa(self) -> str | None: ...
        @property
        def end_aa(self) -> str | None: ...
        @property
        def charge_type(self) -> ChargeType: ...
        @property
        def charge(self) -> int | Mods[GlobalChargeCarrier] | None: ...
        @charge.setter
        def charge(self, value: int | str | list[str] | Mods[GlobalChargeCarrier] | None) -> None: ...
        @property
        def charge_adducts(self) -> Mods[GlobalChargeCarrier]: ...
        @property
        def isotope_mods(self) -> Mods[IsotopeReplacement]: ...
        @isotope_mods.setter
        def isotope_mods(self, value: Any) -> None: ...
        @property
        def static_mods(self) -> Mods[FixedModification]: ...
        @static_mods.setter
        def static_mods(self, value: Any) -> None: ...
        @property
        def labile_mods(self) -> Mods[ModificationTags]: ...
        @labile_mods.setter
        def labile_mods(self, value: Any) -> None: ...
        @property
        def unknown_mods(self) -> Mods[ModificationTags]: ...
        @unknown_mods.setter
        def unknown_mods(self, value: Any) -> None: ...
        @property
        def nterm_mods(self) -> Mods[ModificationTags]: ...
        @nterm_mods.setter
        def nterm_mods(self, value: Any) -> None: ...
        @property
        def cterm_mods(self) -> Mods[ModificationTags]: ...
        @cterm_mods.setter
        def cterm_mods(self, value: Any) -> None: ...
        @property
        def internal_mods(self) -> dict[int, Mods[ModificationTags]]: ...
        @internal_mods.setter
        def internal_mods(self, value: dict[int, Any] | None) -> None: ...
        @property
        def intervals(self) -> tuple[Interval, ...]: ...
        @intervals.setter
        def intervals(self, value: list[Interval] | None) -> None: ...
        def copy(self) -> Self: ...

    """
    Validators
    """

    def validate_sequence(self) -> None:
        """Check that every residue in the sequence is a recognised amino acid.

        :raises PeptacularError: If an unrecognised amino acid code is found.
        """
        for aa in self.sequence:
            if aa not in AA_LOOKUP:
                raise PeptacularError(f"Invalid amino acid '{aa}' in sequence '{self.sequence}'")

    def validate_isotope_mods(self) -> None:
        """Check that all isotope modifications are structurally valid.

        :raises PeptacularError: If any isotope modification is invalid.
        """
        if errors := self.isotope_mods.validate():
            raise PeptacularError(f"Invalid isotope modifications: {errors}")

    def validate_static_mods(self) -> None:
        """Check that all static (fixed) modifications are structurally valid.

        :raises PeptacularError: If any static modification is invalid.
        """
        if errors := self.static_mods.validate():
            raise PeptacularError(f"Invalid static modifications: {errors}")

    def validate_labile_mods(self) -> None:
        """Check that all labile modifications are structurally valid.

        :raises PeptacularError: If any labile modification is invalid.
        """
        if errors := self.labile_mods.validate():
            raise PeptacularError(f"Invalid labile modifications: {errors}")

    def validate_unknown_mods(self) -> None:
        """Check that all unknown-localisation modifications are structurally valid.

        :raises PeptacularError: If any unknown modification is invalid.
        """
        if errors := self.unknown_mods.validate():
            raise PeptacularError(f"Invalid unknown modifications: {errors}")

    def validate_nterm_mods(self) -> None:
        """Check that all N-terminal modifications are structurally valid.

        :raises PeptacularError: If any N-terminal modification is invalid.
        """
        if errors := self.nterm_mods.validate():
            raise PeptacularError(f"Invalid N-terminal modifications: {errors}")

    def validate_cterm_mods(self) -> None:
        """Check that all C-terminal modifications are structurally valid.

        :raises PeptacularError: If any C-terminal modification is invalid.
        """
        if errors := self.cterm_mods.validate():
            raise PeptacularError(f"Invalid C-terminal modifications: {errors}")

    def validate_internal_mods(self) -> None:
        """Check that all internal (per-position) modifications are structurally valid.

        :raises PeptacularError: If any internal modification at any position is invalid.
        """
        for pos, mods in self.internal_mods.items():
            if errors := mods.validate():
                raise PeptacularError(f"Invalid internal modifications at position {pos}: {errors}")

    def validate_intervals(self) -> None:
        """Check that all intervals are valid, non-overlapping, and within sequence bounds.

        :raises PeptacularError: If any interval is invalid, intervals overlap, or an interval
            falls outside the sequence length.
        """
        intervals = self.intervals
        for interval in intervals:
            if errors := interval.validate():
                raise PeptacularError(f"Invalid interval: {errors}")

        # ensure no overlapping intervals
        sorted_intervals = sorted(intervals, key=lambda x: x.start)
        for i in range(1, len(sorted_intervals)):
            if sorted_intervals[i].start < sorted_intervals[i - 1].end:
                raise PeptacularError(f"Overlapping intervals detected: {sorted_intervals[i - 1]} and {sorted_intervals[i]}")

        # ensure that intervals dont start/end out of bounds
        seq_len = len(self.sequence) if self._sequence is not None else 0
        for interval in intervals:
            if interval.start < 0 or interval.end > seq_len:
                raise PeptacularError(f"Interval {interval} is out of bounds for sequence length {seq_len}")

    def validate_ambiguous_labels(self) -> None:
        """Check that each ambiguous-position label (``#label``) has at most one
        concrete modification among its occurrences; the rest must be bare
        references (e.g. ``[#label]``).

        :raises PeptacularError: If a label has more than one concrete occurrence.
        """
        concrete_label_counts: Counter[str] = Counter()

        def scan(mods: "Mods | None") -> None:
            for position_id in _concrete_position_labels(mods):
                concrete_label_counts[position_id] += 1

        if self.has_internal_mods:
            for mods in self.internal_mods.values():
                scan(mods)
        if self.has_nterm_mods:
            scan(self.nterm_mods)
        if self.has_cterm_mods:
            scan(self.cterm_mods)
        if self.has_intervals:
            for interval in self.intervals:
                scan(interval.mods)

        duplicated = sorted(label for label, count in concrete_label_counts.items() if count > 1)
        if duplicated:
            raise PeptacularError(
                f"Ambiguous modification label(s) {duplicated} have more than one concrete modification; "
                "exactly one occurrence of a labelled group may carry the modification text, "
                "others must be bare references (e.g. [#label])."
            )

    def validate_charge(self) -> None:
        """Check that the charge value is structurally valid.

        :raises PeptacularError: If the charge adducts are invalid or the charge type is
            unrecognised.
        """
        charge_type = self.charge_type

        match charge_type:
            case ChargeType.INT:
                pass
            case ChargeType.ADDUCTS:
                if errors := self.charge_adducts.validate():
                    raise PeptacularError(f"Invalid charge adducts: {errors}")
            case ChargeType.NONE:
                pass
            case _:
                raise PeptacularError(f"Invalid charge type: {charge_type}")

    def validate_annotation(self) -> None:
        """Run all individual validators in order; raises on the first error found.

        :raises PeptacularError: If any component of the annotation is structurally invalid.
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
        self.validate_ambiguous_labels()
        self.validate_charge()

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

    """
    Set Methods - Replace existing modifications
    """

    def set_sequence(self, sequence: str | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def set_compound_name(self, name: str | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def set_ion_name(self, name: str | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def set_peptide_name(self, name: str | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def set_isotope_mods(self, mods: dict[str, int] | Mods[IsotopeReplacement] | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def set_static_mods(self, mods: dict[str, int] | Mods[FixedModification] | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def set_labile_mods(self, mods: dict[str, int] | Mods[ModificationTags] | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def set_unknown_mods(self, mods: dict[str, int] | Mods[ModificationTags] | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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
        self, mods: dict[str, int] | Mods[ModificationTags] | None, *, inplace: bool = True, validate: bool | None = None, start_aa: str | None = None
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
        self, mods: dict[str, int] | Mods[ModificationTags] | None, *, inplace: bool = True, validate: bool | None = None, end_aa: str | None = None
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
        self, mods: dict[int, dict[str, int] | Mods[ModificationTags] | None] | None, *, inplace: bool = True, validate: bool | None = None
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
            self.validate_ambiguous_labels()
        return self

    def set_intervals(self, intervals: list[Interval] | None, *, inplace: bool = True, validate: bool | None = None) -> Self:
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
            self.validate_ambiguous_labels()
        return self

    def set_internal_mods_at_index(self, index: int, mods: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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
                raise PeptacularError(f"Invalid internal modifications at position {index}")

        if self._internal_mods is None:
            self._internal_mods = {}

        self._internal_mods[index] = mods
        # The ambiguous-label invariant is global (a label may have at most one concrete
        # modification across the whole annotation), but a full re-scan on every single-index
        # set makes residue-by-residue construction O(n^2). Mods that carry no concrete
        # position label can't introduce a new violation, so only re-validate when they do.
        new_mods = Mods[ModificationTags](mod_type=ModType.INTERNAL, _mods=mods)
        if validate and any(True for _ in _concrete_position_labels(new_mods)):
            self.validate_ambiguous_labels()
        return self

    def set_charge(
        self,
        charge: int | str | list[str] | tuple[str, ...] | Mods[GlobalChargeCarrier] | GlobalChargeCarrier | Mod[GlobalChargeCarrier] | None,
        *,
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        """Replace the charge value.

        :param charge: New charge as an integer, adduct string(s), ``Mods``, or ``None`` to clear.
        :type charge: int | str | list[str] | tuple[str, ...] | Mods[GlobalChargeCarrier] | GlobalChargeCarrier | Mod[GlobalChargeCarrier] | None
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        :raises PeptacularError: If the resolved charge value has an unsupported type.
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().set_charge(charge, inplace=True, validate=validate)

        set_value: None | int | list[str] = None
        if isinstance(charge, bool):
            # bool is an int subclass; guard it before the int branch so True/False don't
            # slip through and serialize as a garbage charge like 'PEPTIDE/True'.
            raise PeptacularError(f"Unsupported charge type: {type(charge)!r}")
        elif isinstance(charge, int):
            # A charge of 0 is a neutral peptidoform (no charge component per ProForma 2.1
            # section 11.5), so clear it to None rather than storing a literal 0.
            set_value = charge if charge != 0 else None
        elif isinstance(charge, str):
            set_value = [charge]
        elif isinstance(charge, (list, tuple)):
            if len(charge) == 0:
                set_value = None
            else:
                set_value = [str(c) for c in charge]
                if len(set_value) == 0:
                    set_value = None
        elif charge is None:
            set_value = None
        elif isinstance(charge, Mods):
            # Expand each carrier by its occurrence count so repeated adducts survive
            # the round-trip into the ``list[str]`` storage; iterating keys alone would
            # drop the count and silently reduce a multi-adduct charge to one carrier.
            set_value = [str(c) for c, n in charge._mods.items() for _ in range(n)] if charge._mods else None
        elif isinstance(charge, Mod):
            # A Mod wraps a charge carrier value; str(Mod) would emit the dataclass repr
            # (e.g. "Mod(value=GlobalChargeCarrier(...), count=1)"), which is not a valid
            # charge carrier. Serialize the wrapped carrier itself (it already encodes its
            # own occurrence, e.g. "Na:z+1^2"), repeated by the Mod's count.
            # A count of 0 is a neutral peptidoform: clear to None (matching the empty
            # list/int-zero branches) rather than storing [] and serializing "PEPTIDE/[]".
            set_value = [str(charge.value)] * charge.count if charge.count > 0 else None
        elif isinstance(charge, GlobalChargeCarrier):
            set_value = [str(charge)]
        else:
            raise PeptacularError(f"Unsupported charge type: {type(charge)!r}")

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

    def set_mods(self, mods: Mapping[ModType | ModTypeLiteral | int, Any] | None, *, inplace: bool = True) -> Self:
        """Set a modification by type, replacing any existing mods of that type"""

        if not inplace:
            return self.copy().set_mods(mods=mods, inplace=True)

        if mods is None:
            self.clear_mods(inplace=True)
            return self

        for mod_type, mod_value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise InvalidPositionError(f"Internal modification index out of range: {mod_type}")
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

        if is_mod_collection(mod):
            for item in mod:
                self._append_mod_generic(item, attr_name, validator, inplace=True, validate=validate)
            return self

        mod_str, count = convert_single_mod_input(mod)

        if validate:
            if not validator(mod_str).is_valid:
                raise PeptacularError(f"Invalid modification: {mod_str}")

        mod_dict = getattr(self, attr_name)
        if mod_dict is None:
            setattr(self, attr_name, {})
            mod_dict = getattr(self, attr_name)

        if mod_str in mod_dict:
            mod_dict[mod_str] += count
        else:
            mod_dict[mod_str] = count

        return self

    def append_isotope_mod(self, mod: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def append_static_mod(self, mod: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def append_labile_mod(self, mod: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def append_unknown_mod(self, mod: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def append_nterm_mod(self, mod: Any, *, inplace: bool = True, validate: bool | None = None, start_aa: str | None = None) -> Self:
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

    def append_cterm_mod(self, mod: Any, *, inplace: bool = True, validate: bool | None = None, end_aa: str | None = None) -> Self:
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

    def append_internal_mod_at_index(self, index: int, mod: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

        if is_mod_collection(mod):
            for item in mod:
                self.append_internal_mod_at_index(index, item, inplace=True, validate=validate)
            return self

        mod_str, count = convert_single_mod_input(mod)

        if validate:
            if not ModificationTags.from_string(mod_str).is_valid:
                raise PeptacularError(f"Invalid modification: {mod_str}")

        if self._internal_mods is None:
            self._internal_mods = {}

        if index not in self._internal_mods:
            self._internal_mods[index] = {}

        if mod_str in self._internal_mods[index]:
            self._internal_mods[index][mod_str] += count
        else:
            self._internal_mods[index][mod_str] = count

        return self

    def append_interval(self, interval: Interval | tuple[int, int, bool, Any], *, inplace: bool = True, validate: bool | None = None) -> Self:
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
                raise PeptacularError(f"Invalid interval: {interval}")

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

    def append_mods(self, mods: Mapping[ModType | ModTypeLiteral | int, Any], *, inplace: bool = True, validate: bool | None = None) -> Self:
        """Append modifications of multiple types from a mapping of mod-type to value.

        :param mods: Mapping of :class:`ModType` (or literal/index) to a modification value, or a list/tuple of values to append each of.
        :type mods: Mapping[ModType | ModTypeLiteral | int, Any]
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        :raises InvalidPositionError: If an integer key is out of range for the current sequence.
        """
        if not inplace:
            return self.copy().append_mods(mods, inplace=True, validate=validate)

        for mod_type, value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise InvalidPositionError(f"Internal modification index out of range: {mod_type}")
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
        append_method: Callable[..., Self],
        inplace: bool = True,
        validate: bool | None = None,
    ) -> Self:
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy()._extend_generic(mods, append_method, inplace=True, validate=validate)
        if mods is not None:
            for mod in as_mod_iterable(mods):
                append_method(mod, inplace=True, validate=validate)
        return self

    def extend_isotope_mods(self, mods: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def extend_static_mods(self, mods: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def extend_labile_mods(self, mods: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def extend_unknown_mods(self, mods: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def extend_nterm_mods(self, mods: Any, *, inplace: bool = True, validate: bool | None = None, start_aa: str | None = None) -> Self:
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
            for mod in as_mod_iterable(mods):
                self.append_nterm_mod(mod, inplace=True, validate=validate, start_aa=start_aa)
        return self

    def extend_cterm_mods(self, mods: Any, *, inplace: bool = True, validate: bool | None = None, end_aa: str | None = None) -> Self:
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
            for mod in as_mod_iterable(mods):
                self.append_cterm_mod(mod, inplace=True, validate=validate, end_aa=end_aa)
        return self

    def extend_internal_mods_at_index(self, index: int, mods: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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
            for mod in as_mod_iterable(mods):
                self.append_internal_mod_at_index(index, mod, inplace=True, validate=validate)
        return self

    def extend_intervals(self, intervals: Any, *, inplace: bool = True, validate: bool | None = None) -> Self:
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

    def extend_mods(self, mods: Mapping[ModType | ModTypeLiteral | int, Any], *, inplace: bool = True, validate: bool | None = None) -> Self:
        """Extend modifications of multiple types by iterating through each mapped iterable.

        :param mods: Mapping of :class:`ModType` (or literal/index) to iterable of modification values. A bare string is one modification.
        :type mods: Mapping[ModType | ModTypeLiteral | int, Any]
        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :param validate: Override the instance-level validation flag for this call only.
        :type validate: bool | None
        :return: The (possibly new) annotation.
        :rtype: Self
        :raises InvalidPositionError: If an integer key is out of range for the current sequence.
        """
        if validate is None:
            validate = self._validate
        if not inplace:
            return self.copy().extend_mods(mods, inplace=True, validate=validate)

        for mod_type, value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise InvalidPositionError(f"Internal modification index out of range: {mod_type}")
                self.extend_internal_mods_at_index(mod_type, value, inplace=True, validate=validate)
                continue
            self._extend_by_type(value, ModType(mod_type), inplace=True, validate=validate)

        return self

    """
    REMOVE Methods
    """

    def remove_mods(self, mods: Mapping[ModType | ModTypeLiteral | int, Any], *, inplace: bool = True) -> Self:
        """Remove modifications by decrementing their counts."""
        if not inplace:
            return self.copy().remove_mods(mods, inplace=True)

        for mod_type, mod_value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise InvalidPositionError(f"Internal modification index out of range: {mod_type}")
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

    def remove_isotope_mod(self, mod: Any, *, inplace: bool = True) -> Self:
        """Remove a specific isotope modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_isotope_mods", inplace)

    def remove_static_mod(self, mod: Any, *, inplace: bool = True) -> Self:
        """Remove a specific static modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_static_mods", inplace)

    def remove_labile_mod(self, mod: Any, *, inplace: bool = True) -> Self:
        """Remove a specific labile modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_labile_mods", inplace)

    def remove_unknown_mod(self, mod: Any, *, inplace: bool = True) -> Self:
        """Remove a specific unknown modification by decrementing its count."""
        return self._remove_mod_generic(mod, "_unknown_mods", inplace)

    def remove_nterm_mod(self, mod: Any, *, inplace: bool = True, start_aa: str | None = None) -> Self:
        """Remove a specific N-terminal modification by decrementing its count."""
        if start_aa is not None and self.start_aa != start_aa:
            return self if inplace else self.copy()
        return self._remove_mod_generic(mod, "_nterm_mods", inplace)

    def remove_cterm_mod(self, mod: Any, *, inplace: bool = True, end_aa: str | None = None) -> Self:
        """Remove a specific C-terminal modification by decrementing its count."""
        if end_aa is not None and self.end_aa != end_aa:
            return self if inplace else self.copy()
        return self._remove_mod_generic(mod, "_cterm_mods", inplace)

    def remove_internal_mod_at_index(self, index: int, mod: Any, *, inplace: bool = True) -> Self:
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

    def remove_interval(self, interval: Interval, *, inplace: bool = True) -> Self:
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
        mod_enums = _resolve_mod_types(mod_types)
        return any(self._has_mods_by_type(mod_enum) for mod_enum in mod_enums)

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
        mod_enums = _resolve_mod_types(mod_types)
        return {mod_enum: self._get_mods_by_type(mod_enum) for mod_enum in mod_enums if self._has_mods_by_type(mod_enum)}

    """
    Pop Methods
    """

    def pop_isotope_mods(self, *, inplace: bool = True) -> Mods[IsotopeReplacement]:
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

    def pop_static_mods(self, *, inplace: bool = True) -> Mods[FixedModification]:
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

    def pop_labile_mods(self, *, inplace: bool = True) -> Mods[ModificationTags]:
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

    def pop_unknown_mods(self, *, inplace: bool = True) -> Mods[ModificationTags]:
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

    def pop_nterm_mods(self, *, inplace: bool = True) -> Mods[ModificationTags]:
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

    def pop_cterm_mods(self, *, inplace: bool = True) -> Mods[ModificationTags]:
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

    def pop_internal_mods(self, *, inplace: bool = True) -> dict[int, Mods[ModificationTags]]:
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

    def pop_intervals(self, *, inplace: bool = True) -> list[Interval]:
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

    def pop_charge(self, *, inplace: bool = True) -> int | Mods[GlobalChargeCarrier] | None:
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

    def pop_internal_mod_at_index(self, index: int, *, inplace: bool = True) -> tuple[tuple[MODIFICATION_TYPE, int], ...]:
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
        self, mod_types: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None = None, *, inplace: bool = True
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

        mod_enums: list[ModType] = _resolve_mod_types(mod_types)

        d: dict[ModType, Any] = {}
        for mod_enum in mod_enums:
            d[mod_enum] = self._pop_mod_by_type(mod_enum)

        return d

    def filter_mods(
        self, mods: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None = None, *, inplace: bool = True, keep: bool = True
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
            return self.copy().filter_mods(mods=mods, inplace=True, keep=keep)

        if keep:
            # Keep only specified mods
            mod_types_to_keep = set(_resolve_mod_types(mods))

            all_mod_types = {mod_type for mod_type in ModType}
            mod_types_to_remove = all_mod_types - mod_types_to_keep
        else:
            # Remove only specified mods
            mod_types_to_remove = _resolve_mod_types(mods)

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

    def clear_isotope_mods(self, *, inplace: bool = True) -> Self:
        """Clear all isotope modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_isotope_mods", inplace)

    def clear_static_mods(self, *, inplace: bool = True) -> Self:
        """Clear all static modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_static_mods", inplace)

    def clear_nterm_mods(self, *, inplace: bool = True) -> Self:
        """Clear all N-terminal modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_nterm_mods", inplace)

    def clear_cterm_mods(self, *, inplace: bool = True) -> Self:
        """Clear all C-terminal modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_cterm_mods", inplace)

    def clear_labile_mods(self, *, inplace: bool = True) -> Self:
        """Clear all labile modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_labile_mods", inplace)

    def clear_unknown_mods(self, *, inplace: bool = True) -> Self:
        """Clear all unknown-localisation modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_unknown_mods", inplace)

    def clear_internal_mods(self, *, inplace: bool = True) -> Self:
        """Clear all internal (per-position) modifications.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_internal_mods", inplace)

    def clear_internal_mod_at_index(self, index: int, *, inplace: bool = True) -> Self:
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

    def clear_intervals(self, *, inplace: bool = True) -> Self:
        """Clear all ambiguous sequence intervals.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) annotation.
        :rtype: Self
        """
        return self._clear_mod_dict("_intervals", inplace)

    def clear_charge(self, *, inplace: bool = True) -> Self:
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

    def clear_mods(self, mods: ModTypeLiteral | ModType | Iterable[ModTypeLiteral] | Iterable[ModType] | None = None, *, inplace: bool = True) -> Self:
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
        mod_enums = _resolve_mod_types(mods)
        for mod_enum in mod_enums:
            self._clear_mod_by_type(mod_enum)
        return self

    def strip_mods(self, *, inplace: bool = False) -> Self:
        """Remove all modifications of every type, leaving only the bare sequence.

        :param inplace: Modify this object when ``True``; return a modified copy when ``False``.
        :type inplace: bool
        :return: The (possibly new) bare annotation.
        :rtype: Self
        """
        return self.clear_mods(None, inplace=inplace)

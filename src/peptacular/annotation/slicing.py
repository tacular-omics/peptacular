from __future__ import annotations

import random
from collections.abc import Callable, Generator, Sequence
from typing import TYPE_CHECKING, Any

from ..annotation.mod import Mods
from ..constants import Terminal
from ..proforma_components.comps import FixedModification, ModificationTags
from .parser import Interval, ProFormaParser

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def slice_annotation(
    annotation: ProFormaAnnotation,
    start: int | None,
    stop: int | None,
    inplace: bool = False,
) -> ProFormaAnnotation:
    """
    Slice the annotation sequence and return a new annotation with the sliced sequence and modifications.

    Args:
        annotation: The ProFormaAnnotation to slice
        start: Start index (inclusive)
        stop: Stop index (exclusive)
        inplace: Whether to modify the annotation in place

    Returns:
        Sliced annotation

    Raises:
        ValueError: If indices are out of bounds or invalid
    """
    if not inplace:
        return slice_annotation(
            annotation.copy(),
            start=start,
            stop=stop,
            inplace=True,
        )

    # Handle default values and negative indices
    seq_len = len(annotation)
    start = _normalize_start_index(start, seq_len)
    stop = _normalize_stop_index(stop, seq_len)
    _validate_slice_indices(start, stop, seq_len)

    new_sequence = annotation.stripped_sequence[start:stop]

    # Early return if no modifications exist
    if not annotation.has_mods():
        annotation.sequence = new_sequence
        return annotation

    # Adjust internal modifications
    new_internal_mods = None
    if annotation.has_internal_mods:
        if annotation._internal_mods is None:
            raise RuntimeError("Annotation internal_mods are None despite has_internal_mods being True.")
        new_internal_mods = _adjust_internal_mods(
            annotation._internal_mods,
            start,
            stop,
        )

    # Adjust intervals
    new_intervals = None
    if annotation.has_intervals:
        if annotation._intervals is None:
            raise RuntimeError("Annotation intervals are None despite has_intervals being True.")
        new_intervals = _adjust_intervals(annotation._intervals, start, stop)

    # Handle terminal modifications
    if start > 0:
        annotation.clear_nterm_mods(inplace=True)
    if stop < seq_len:
        annotation.clear_cterm_mods(inplace=True)

    # Adjust static modifications for removed terminals
    if annotation.has_static_mods and annotation._static_mods is not None:
        keep_nterm = start == 0
        keep_cterm = stop == seq_len
        if not (keep_nterm and keep_cterm):
            adjusted = _adjust_static_mods(
                annotation._static_mods,
                keep_nterm=keep_nterm,
                keep_cterm=keep_cterm,
            )
            annotation._static_mods = adjusted if adjusted else None

    # Update annotation
    annotation.sequence = new_sequence
    annotation.internal_mods = new_internal_mods
    annotation.intervals = new_intervals
    return annotation


def split_annotation(annotation: ProFormaAnnotation) -> list[ProFormaAnnotation]:
    """
    Split each amino acid in the sequence into a separate ProFormaAnnotation.

    Args:
        annotation: The ProFormaAnnotation to split

    Returns:
        List of single amino acid annotations

    Raises:
        ValueError: If annotation contains intervals
    """

    annot_slices: list[tuple[int, int]] = []
    for i, _ in enumerate(annotation.sequence):
        annot_slices.append((i, i + 1))

    interval_slices: list[tuple[int, int]] = []
    if annotation.has_intervals:
        for interval in annotation.intervals:
            start: int = interval.start
            end: int = interval.end
            interval_slices.append((start, end))

    for start, stop in interval_slices:
        annot_slices = [(s, e) for s, e in annot_slices if s < start or s >= stop]

    annot_slices += interval_slices

    annot_slices.sort(key=lambda x: x[0])

    result: list[ProFormaAnnotation] = []
    for s, e in annot_slices:
        sliced = slice_annotation(annotation, s, e, inplace=False)
        result.append(sliced)

    return result


def join_annotations(annotations: Sequence[ProFormaAnnotation]) -> ProFormaAnnotation:
    """
    Join multiple ProFormaAnnotations into a single annotation (most likely from a split operation).

    Args:
        annotations: List of ProFormaAnnotations to join

    Returns:
        The joined ProFormaAnnotation
    """
    if not annotations:
        raise ValueError("No annotations to join")

    if len(annotations) == 1:
        return annotations[0].copy()

    # Start with a copy of the first annotation
    result = annotations[0].copy()

    # Concatenate sequences
    full_sequence = "".join(ann.sequence for ann in annotations)
    result.sequence = full_sequence

    # Merge internal modifications, adjusting positions
    new_internal_mods: dict[int, Mods[ModificationTags]] = {}
    current_pos = 0

    for ann in annotations:
        if ann.has_internal_mods:
            for pos, mods in ann.internal_mods.items():
                new_internal_mods[current_pos + pos] = mods.copy()
        current_pos += len(ann.sequence)

    result.internal_mods = new_internal_mods if new_internal_mods else {}

    # Merge intervals, adjusting positions
    new_intervals: list[Interval] = []
    current_pos = 0

    for ann in annotations:
        if ann.has_intervals:
            for interval in ann.intervals:
                new_intervals.append(
                    interval.update(
                        start=current_pos + interval.start,
                        end=current_pos + interval.end,
                    )
                )
        current_pos += len(ann.sequence)

    result.intervals = new_intervals if new_intervals else []

    # Handle terminal modifications - only keep from first (N-term) and last (C-term)
    if len(annotations) > 1:
        last_ann = annotations[-1]
        if last_ann.has_cterm_mods:
            result.cterm_mods = last_ann.cterm_mods.copy()

        first_ann = annotations[0]
        if first_ann.has_nterm_mods:
            result.nterm_mods = first_ann.nterm_mods.copy()

    return result


def shift_annotation(
    annotation: ProFormaAnnotation,
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    inplace: bool = False,
) -> ProFormaAnnotation:
    """
    Shift the annotation by n positions in a cyclic manner.

    Args:
        annotation: The ProFormaAnnotation to shift
        n: Number of positions to shift (positive shifts right, negative shifts left)
        keep_nterm: Number of N-terminal residues to keep unchanged
        keep_cterm: Number of C-terminal residues to keep unchanged
        inplace: If True, modify the annotation in place

    Returns:
        The shifted annotation

    Raises:
        ValueError: If intervals would be cut off and slice_intervals=False
    """

    if not inplace:
        return shift_annotation(
            annotation.copy(),
            n=n,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
            inplace=True,
        )

    if not annotation.has_sequence:
        return annotation

    seq_len = len(annotation.sequence)

    # Validate keep parameters
    if keep_nterm < 0 or keep_cterm < 0:
        raise ValueError("keep_nterm and keep_cterm must be non-negative")
    if keep_nterm + keep_cterm > seq_len:
        raise ValueError("keep_nterm + keep_cterm cannot exceed sequence length")

    # If nothing to shift, return as is
    if keep_nterm + keep_cterm >= seq_len:
        return annotation

    effective_shift = n % seq_len
    if effective_shift == 0:
        return annotation

    # Extract the parts of the sequence
    nterm_part = annotation.sequence[:keep_nterm]
    middle_part = annotation.sequence[keep_nterm : seq_len - keep_cterm]
    cterm_part = annotation.sequence[seq_len - keep_cterm :]

    # Calculate effective shift for the middle part
    middle_len = len(middle_part)
    if middle_len == 0:
        return annotation

    middle_shift = effective_shift % middle_len if middle_len > 0 else 0

    # Shift only the middle part
    shifted_middle = middle_part[middle_shift:] + middle_part[:middle_shift]
    shifted_sequence = nterm_part + shifted_middle + cterm_part

    # Shift internal modifications (only in the middle section)
    new_internal_mods: dict[int, Mods[ModificationTags]] = {}
    if annotation.has_internal_mods:
        for mod_index, mods in annotation.internal_mods.items():
            if mod_index < keep_nterm:
                # Keep N-terminal mods in place
                new_internal_mods[mod_index] = mods.copy()
            elif mod_index >= seq_len - keep_cterm:
                # Keep C-terminal mods in place
                new_internal_mods[mod_index] = mods.copy()
            else:
                # Shift middle mods
                relative_pos = mod_index - keep_nterm
                shifted_relative_pos = (relative_pos - middle_shift) % middle_len
                new_internal_mods[keep_nterm + shifted_relative_pos] = mods.copy()

    # Handle intervals
    new_intervals: list[Interval] = []
    if annotation.has_intervals:
        for interval in annotation.intervals:
            # Check if interval is entirely in kept regions
            if interval.end <= keep_nterm or interval.start >= seq_len - keep_cterm:
                # Keep interval as is
                new_intervals.append(interval)
            elif interval.start >= keep_nterm and interval.end <= seq_len - keep_cterm:
                # Interval is entirely in middle - shift it
                relative_start = interval.start - keep_nterm
                relative_end = interval.end - keep_nterm
                new_start = (relative_start - middle_shift) % middle_len + keep_nterm
                new_end = (relative_end - middle_shift) % middle_len + keep_nterm

                if new_start < new_end:
                    new_intervals.append(Interval(new_start, new_end, interval.ambiguous, interval.mods))
                else:
                    raise ValueError("Shifting intervals that wrap around the sequence end is not supported.")
            else:
                # Interval spans kept and shifted regions
                raise ValueError(f"Interval [{interval.start}:{interval.end}] spans kept and shifted regions. This is not supported.")

    # Update annotation
    annotation.sequence = shifted_sequence
    annotation.internal_mods = new_internal_mods if new_internal_mods else {}
    annotation.intervals = new_intervals
    return annotation


def shuffle_annotation(
    annotation: ProFormaAnnotation,
    seed: Any = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    inplace: bool = False,
) -> ProFormaAnnotation:
    """
    Shuffle the annotation sequence randomly.

    Args:
        annotation: The ProFormaAnnotation to shuffle
        seed: Random seed for reproducible shuffling
        keep_nterm: Number of N-terminal residues to keep unchanged
        keep_cterm: Number of C-terminal residues to keep unchanged
        inplace: If True, modify the annotation in place

    Returns:
        The shuffled annotation
    """
    if not inplace:
        return shuffle_annotation(annotation.copy(), seed, keep_nterm, keep_cterm, inplace=True)

    if seed is not None:
        random.seed(seed)

    seq_len = len(annotation.stripped_sequence)

    # Validate keep parameters
    if keep_nterm < 0 or keep_cterm < 0:
        raise ValueError("keep_nterm and keep_cterm must be non-negative")
    if keep_nterm + keep_cterm > seq_len:
        raise ValueError("keep_nterm + keep_cterm cannot exceed sequence length")

    # If nothing to shuffle, return as is
    if seq_len <= 1 or keep_nterm + keep_cterm >= seq_len:
        return annotation

    # Extract parts
    nterm_part = annotation.stripped_sequence[:keep_nterm]
    middle_part = annotation.stripped_sequence[keep_nterm : seq_len - keep_cterm]
    cterm_part = annotation.stripped_sequence[seq_len - keep_cterm :]

    # Create position mapping for middle part only
    middle_positions = list(enumerate(middle_part, start=keep_nterm))

    # Shuffle only the middle positions
    middle_shuffled = [(i, aa) for i, aa in middle_positions]
    random.shuffle(middle_shuffled)

    # Build shuffled sequence
    shuffled_middle: list[str] = []
    position_mapping: dict[int, int] = {}

    # Keep N-term positions unchanged
    for i in range(keep_nterm):
        position_mapping[i] = i

    # Map shuffled middle positions
    for new_relative_pos, (original_pos, aa) in enumerate(middle_shuffled):
        new_pos = keep_nterm + new_relative_pos
        shuffled_middle.append(aa)
        position_mapping[original_pos] = new_pos

    # Keep C-term positions unchanged
    for i in range(seq_len - keep_cterm, seq_len):
        position_mapping[i] = i

    # Update sequence
    annotation.sequence = nterm_part + "".join(shuffled_middle) + cterm_part

    # Update internal modifications with new positions
    if annotation.has_internal_mods:
        new_internal_mods = {}
        for original_pos, mods in annotation.internal_mods.items():
            new_pos = position_mapping[original_pos]
            new_internal_mods[new_pos] = mods.copy()
        annotation.internal_mods = new_internal_mods if new_internal_mods else None

    # Handle intervals - check if they're in shuffled region
    if annotation.has_intervals:
        for interval in annotation.intervals:
            if interval.start < keep_nterm or interval.end > seq_len - keep_cterm:
                # Interval overlaps with kept regions - check if it's entirely in kept region
                if not (interval.end <= keep_nterm or interval.start >= seq_len - keep_cterm):
                    raise ValueError(f"Interval [{interval.start}:{interval.end}] spans kept and shuffled regions. This is not supported.")
            # If interval is entirely in middle, its positions will be shuffled via the mapping
            # Note: This could lead to non-contiguous intervals, which might not be valid
            # Consider raising an error if intervals exist in the shuffled region
            if keep_nterm < interval.start < seq_len - keep_cterm or keep_nterm < interval.end < seq_len - keep_cterm:
                raise ValueError(
                    f"Shuffling sequences with intervals in the shuffled region is not supported. "
                    f"Interval [{interval.start}:{interval.end}] would be disrupted."
                )

    return annotation


def reverse_annotation(
    annotation: ProFormaAnnotation,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    inplace: bool = False,
) -> ProFormaAnnotation:
    """
    Reverse the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to reverse
        keep_nterm: Number of N-terminal residues to keep unchanged
        keep_cterm: Number of C-terminal residues to keep unchanged
        inplace: If True, modify the annotation in place
        swap_terms: If True, swap N-term and C-term modifications

    Returns:
        The reversed annotation
    """
    if not inplace:
        return reverse_annotation(
            annotation.copy(),
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
            inplace=True,
        )

    seq_len = len(annotation)

    sequence: str = annotation.stripped_sequence

    # Validate keep parameters
    if keep_nterm < 0 or keep_cterm < 0:
        raise ValueError("keep_nterm and keep_cterm must be non-negative")
    if keep_nterm + keep_cterm > seq_len:
        raise ValueError("keep_nterm + keep_cterm cannot exceed sequence length")

    # If nothing to reverse, return as is
    if seq_len <= 1 or keep_nterm + keep_cterm >= seq_len:
        return annotation

    # Extract parts
    nterm_part = sequence[:keep_nterm]
    middle_part = sequence[keep_nterm : seq_len - keep_cterm]
    cterm_part = sequence[seq_len - keep_cterm :]

    # Reverse only the middle part
    reversed_middle = middle_part[::-1]
    annotation.sequence = nterm_part + reversed_middle + cterm_part

    # Reverse internal modifications (only in the middle section)
    if annotation.has_internal_mods:
        new_internal_mods: dict[int, Mods[ModificationTags]] = {}
        middle_len = len(middle_part)

        for pos, mods in annotation.internal_mods.items():
            if pos < keep_nterm:
                # Keep N-terminal mods in place
                new_internal_mods[pos] = mods.copy()
            elif pos >= seq_len - keep_cterm:
                # Keep C-terminal mods in place
                new_internal_mods[pos] = mods.copy()
            else:
                # Reverse middle mods
                relative_pos = pos - keep_nterm
                reversed_relative_pos = middle_len - 1 - relative_pos
                new_internal_mods[keep_nterm + reversed_relative_pos] = mods.copy()

        annotation.internal_mods = new_internal_mods if new_internal_mods else {}

    # Reverse the intervals too (only in the middle section)
    if annotation.has_intervals:
        middle_len = len(middle_part)
        new_intervals: list[Interval] = []

        for interval in annotation.intervals:
            if interval.end <= keep_nterm or interval.start >= seq_len - keep_cterm:
                # Keep intervals in kept regions as is
                new_intervals.append(interval)
            elif interval.start >= keep_nterm and interval.end <= seq_len - keep_cterm:
                # Interval is entirely in middle - reverse it
                relative_start = interval.start - keep_nterm
                relative_end = interval.end - keep_nterm
                new_start = keep_nterm + (middle_len - relative_end)
                new_end = keep_nterm + (middle_len - relative_start)
                new_intervals.append(interval.update(start=new_start, end=new_end))
            else:
                # Interval spans kept and reversed regions
                raise ValueError(f"Interval [{interval.start}:{interval.end}] spans kept and reversed regions. This is not supported.")

        annotation.intervals = new_intervals

    return annotation


def sort_annotation(
    annotation: ProFormaAnnotation,
    inplace: bool = False,
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
) -> ProFormaAnnotation:
    """
    Sort the residues in the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to sort
        inplace: If True, modify the annotation in place
        key: Function to use for sorting key (default: alphabetical)
        reverse: If True, sort in descending order

    Returns:
        The sorted annotation
    """
    if not inplace:
        return sort_annotation(annotation.copy(), inplace=True, key=key, reverse=reverse)

    seq_len = len(annotation.sequence)
    if seq_len <= 1:
        return annotation

    # Store all modification types to restore later
    n_term_mods = annotation.pop_nterm_mods()
    c_term_mods = annotation.pop_cterm_mods()
    labile_mods = annotation.pop_labile_mods()
    unknown_mods = annotation.pop_unknown_mods()
    static_mods = annotation.pop_static_mods()
    isotope_mods = annotation.pop_isotope_mods()
    charge = annotation.pop_charge()

    # Split, sort, and rejoin
    sorted_components = split_annotation(annotation)

    # Create list of (component, original_position) tuples
    sequence_elems = list(enumerate(comp.serialize() for comp in sorted_components))

    # Sort based on key function or alphabetically
    if key is None:
        sequence_elems.sort(key=lambda x: x[1], reverse=reverse)
    else:
        sequence_elems.sort(key=lambda x: key(x[1]), reverse=reverse)

    new_sequence = "".join(comp for _, comp in sequence_elems)

    first: ProFormaParser = ProFormaParser(new_sequence).parse().__next__()[0]

    annotation.sequence = "".join(first.amino_acids)
    annotation.internal_mods = first.internal_mods
    annotation.intervals = first.intervals

    # Restore terminal and global modifications
    annotation.nterm_mods = n_term_mods
    annotation.cterm_mods = c_term_mods
    annotation.labile_mods = labile_mods
    annotation.unknown_mods = unknown_mods
    annotation.static_mods = static_mods
    annotation.isotope_mods = isotope_mods
    annotation.charge = charge

    return annotation


def generate_sliding_windows(
    annotation: ProFormaAnnotation,
    window_size: int,
    reverse: bool = False,
) -> Generator[ProFormaAnnotation]:
    """
    Generate sliding windows of the annotation with a specified size.

    Args:
        annotation: The ProFormaAnnotation to generate windows from
        window_size: Size of each window
        keep_terms: Whether to keep terminal modifications in windows
        slice_intervals: Whether to slice intervals at window boundaries
        keep_labile: Whether to keep labile modifications in windows
        reverse: If True, generate windows from right to left

    Yields:
        ProFormaAnnotation instances representing each window

    Raises:
        ValueError: If window_size is invalid or annotation has no sequence
    """
    if not annotation.has_sequence:
        raise ValueError("Annotation must have a sequence to create sliding windows")

    if window_size <= 0:
        raise ValueError("Window size must be positive")

    seq_len = len(annotation.sequence)
    if window_size > seq_len:
        raise ValueError(f"Window size {window_size} cannot be greater than sequence length {seq_len}.")

    if reverse:
        # Generate windows from right to left
        for start in range(seq_len - window_size, -1, -1):
            stop = start + window_size
            yield slice_annotation(
                annotation,
                start=start,
                stop=stop,
                inplace=False,
            )
    else:
        # Generate windows from left to right
        for start in range(seq_len - window_size + 1):
            stop = start + window_size
            yield slice_annotation(
                annotation,
                start=start,
                stop=stop,
                inplace=False,
            )


def _normalize_start_index(start: int | None, seq_len: int) -> int:
    """Normalize start index handling None and negative values"""
    if start is None:
        return 0
    elif start < 0:
        start = seq_len + start
        if start < 0:
            raise ValueError("Start index is out of bounds for the sequence length.")
        return start
    elif start > seq_len:
        raise ValueError(f"Start index exceeds the sequence length. Sequence length is {seq_len}, but start is {start}.")
    return start


def _normalize_stop_index(stop: int | None, seq_len: int) -> int:
    """Normalize stop index handling None and negative values"""
    if stop is None:
        return seq_len
    elif stop < 0:
        stop = seq_len + stop
        if stop < 0:
            raise ValueError("Stop index is out of bounds for the sequence length.")
        return stop
    elif stop > seq_len:
        raise ValueError(f"Stop index exceeds the sequence length. Sequence length is {seq_len}, but stop is {stop}.")
    return stop


def _validate_slice_indices(start: int, stop: int, seq_len: int) -> None:
    """Validate that slice indices are valid"""
    if start > stop:
        raise ValueError(f"Start index {start} cannot be greater than stop index {stop}.")


def _adjust_internal_mods(internal_mods: dict[int, dict[str, int]], start: int, stop: int) -> dict[int, dict[str, int]]:
    """Adjust internal modifications for slicing"""
    new_internal_mods: dict[int, dict[str, int]] = {}
    for pos, mods in internal_mods.items():
        if start <= pos < stop:
            new_internal_mods[pos - start] = mods
    return new_internal_mods


def _adjust_static_mods(
    static_mods: dict[str, int],
    keep_nterm: bool,
    keep_cterm: bool,
) -> dict[str, int]:
    """Adjust static modifications when a terminal is removed by slicing.

    :param static_mods: Serialized FixedModification keys mapped to counts.
    :type static_mods: dict[str, int]
    :param keep_nterm: True if the N-terminus is retained in the slice.
    :type keep_nterm: bool
    :param keep_cterm: True if the C-terminus is retained in the slice.
    :type keep_cterm: bool
    :return: New dict with terminal-specific rules filtered out as appropriate.
    :rtype: dict[str, int]
    """
    new_static_mods: dict[str, int] = {}
    for mod_key, count in static_mods.items():
        # Quick path: no terminal markers in the key string means no terminal rules to filter
        has_nterm_rule = "N-term" in mod_key
        has_cterm_rule = "C-term" in mod_key
        if not has_nterm_rule and not has_cterm_rule:
            new_static_mods[mod_key] = count
            continue

        fixed_mod = FixedModification.from_string(mod_key)

        # Empty position_rules means "applies everywhere" — always keep
        if not fixed_mod.position_rules:
            new_static_mods[mod_key] = count
            continue

        filtered_rules = tuple(
            rule
            for rule in fixed_mod.position_rules
            if not (rule.terminal == Terminal.N_TERM and not keep_nterm) and not (rule.terminal == Terminal.C_TERM and not keep_cterm)
        )

        if not filtered_rules:
            continue  # all rules were for removed terminals; drop mod

        if len(filtered_rules) == len(fixed_mod.position_rules):
            new_static_mods[mod_key] = count  # no change
        else:
            new_fixed_mod = FixedModification(
                modifications=fixed_mod.modifications,
                position_rules=filtered_rules,
            )
            new_key = new_fixed_mod.serialize()
            new_static_mods[new_key] = new_static_mods.get(new_key, 0) + count

    return new_static_mods


def _adjust_intervals(intervals: list[Interval], start: int, stop: int) -> list[Interval]:
    """Adjust intervals for slicing"""
    new_intervals: list[Interval] = []
    for interval in intervals:
        # Check if interval overlaps with the slice
        if interval.start < stop and interval.end > start:
            # Check if the interval gets cut off by the slice boundaries
            interval_gets_cut = (interval.start < start) or (interval.end > stop)

            if interval_gets_cut:
                raise ValueError(f"Interval [{interval.start}:{interval.end}] would be cut by slice [{start}:{stop}]. Slicing intervals is not supported.")

            # Calculate new positions relative to slice
            new_interval_start = max(0, interval.start - start)
            new_interval_end = min(stop - start, interval.end - start)

            # Only add if the interval has meaningful content in the slice
            if new_interval_start < new_interval_end:
                new_intervals.append(
                    interval.update(
                        start=new_interval_start,
                        end=new_interval_end,
                    )
                )
    return new_intervals

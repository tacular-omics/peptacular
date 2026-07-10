from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def condense_mods_to_intervals(annotation: ProFormaAnnotation, inplace: bool = True) -> ProFormaAnnotation:
    """
    Condense modifications into intervals
    """
    if inplace is False:
        return condense_mods_to_intervals(annotation.copy(), inplace=True)

    if not annotation.has_intervals or not annotation.has_internal_mods:
        return annotation

    indexes_to_remove: list[int] = []
    for index, mods in annotation.internal_mods.items():
        for interval in annotation.intervals:
            if interval.start <= index < interval.end:
                # Merge modifications from internal_mods into interval
                for mod in mods:
                    interval.append_mod(mod, inplace=True)

                indexes_to_remove.append(index)

    for index in indexes_to_remove:
        annotation.clear_internal_mod_at_index(index, inplace=True)

    return annotation


def coverage(
    annotation: ProFormaAnnotation,
    annotations: Iterable[ProFormaAnnotation],
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> list[int]:
    cov_arr = [0] * len(annotation)

    for sub_annots in annotations:
        peptide_cov = [1] * len(sub_annots)
        if ignore_ambiguity is False:
            # Filter for ambiguous intervals only
            for interval in sub_annots.intervals:
                if interval.ambiguous:
                    for i in range(interval.start, interval.end):
                        peptide_cov[i] = 0

        for subsequence_index in find_indices(sub_annots, other=annotation, ignore_mods=ignore_mods):
            start = subsequence_index

            for i, cov in enumerate(peptide_cov):
                if accumulate:
                    cov_arr[start + i] += cov
                else:
                    cov_arr[start + i] = cov

    # apply to annotation
    if ignore_ambiguity is False:
        # Set coverage to 0 for ambiguous positions
        for interval in annotation.intervals:
            if interval.ambiguous:
                for i in range(interval.start, interval.end):
                    cov_arr[i] = 0

    return cov_arr


def percent_coverage(
    annotation: ProFormaAnnotation,
    annotations: Iterable[ProFormaAnnotation],
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> float:
    cov_arr = coverage(
        annotation,
        annotations,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )
    if len(cov_arr) == 0:
        return 0.0

    return sum(cov_arr) / len(cov_arr)


def modification_coverage(
    annotation: ProFormaAnnotation,
    annotations: Iterable[ProFormaAnnotation],
    ignore_ambiguity: bool = False,
    accumulate: bool = False,
) -> dict[int, int]:
    cov_dict = coverage(
        annotation,
        annotations,
        accumulate=accumulate,
        ignore_mods=False,
        ignore_ambiguity=ignore_ambiguity,
    )
    mod_cov: dict[int, int] = {}
    for pos in annotation.internal_mods.keys():
        mod_cov[pos] = cov_dict[pos]

    return mod_cov


def condense_static_mods(annotation: ProFormaAnnotation, inplace: bool = True) -> ProFormaAnnotation:
    """
    Condense static mods into internal mods

    Note: This function currently has limited functionality as the static mod
    parsing is not yet fully implemented in the new structure.
    """
    if inplace is False:
        return condense_static_mods(annotation.copy(), inplace=True)

    if not annotation.has_static_mods:
        return annotation

    static_mod_map = annotation.map_static_mods_to_indexes()
    for index, mods in static_mod_map.items():
        if index == -2:  # c term
            annotation.extend_cterm_mods(mods, inplace=True)
        elif index == -1:  # n term
            annotation.extend_nterm_mods(mods, inplace=True)
        else:
            annotation.extend_internal_mods_at_index(index, mods, inplace=True)

    annotation.clear_static_mods(inplace=True)

    return annotation


def condense_to_peptidoform(annotation: ProFormaAnnotation, inplace: bool = True) -> ProFormaAnnotation:
    """
    Condense all modifications into peptidoform notation
    PEP[+100]TI[+20]DE-[30]

    {+100|info:P}{+20|info:T}{-30|info:Nterm}PEPTIDE
    """

    if inplace is False:
        return condense_to_peptidoform(annotation.copy(), inplace=True)

    # pop intervals

    # Collect all modifications as strings with their counts
    if annotation.has_internal_mods:
        internal_mods = annotation.pop_internal_mods()
        for _, mods in internal_mods.items():
            annotation.extend_unknown_mods(mods)

    if annotation.has_intervals:
        intervals = annotation.pop_intervals()

        for interval in intervals:
            if not interval.has_mods:
                continue
            annotation.extend_unknown_mods(interval.mods)

    return annotation


def count_residues(annotation: ProFormaAnnotation, include_mods: bool = True) -> dict[str, int]:
    """
    Count the occurrences of each residue in the sequence.

    Args:
        annotation: The ProFormaAnnotation to analyze
        include_mods: Whether to include modifications in the count

    Returns:
        Dictionary mapping residue strings to their counts
    """

    if not include_mods:
        return dict(Counter(annotation.sequence))

    split_annotations = annotation.split()
    return dict(Counter(a.serialize() for a in split_annotations))


def percent_residues(
    annotation: ProFormaAnnotation,
    include_mods: bool = True,
) -> dict[str, float]:
    """
    Calculate the percentage of each residue in the sequence.

    Args:
        annotation: The ProFormaAnnotation to analyze
        include_mods: Whether to include modifications in the calculation
        precision: Number of decimal places to round to

    Returns:
        Dictionary mapping residue strings to their percentages
    """
    residue_counts = count_residues(annotation, include_mods=include_mods)
    total_count = sum(residue_counts.values())

    if total_count == 0:
        return {}

    percentages = {residue: (count / total_count) * 100 for residue, count in residue_counts.items()}

    return percentages


def is_subsequence(
    annotation: ProFormaAnnotation,
    other: ProFormaAnnotation,
    ignore_mods: bool = False,
    ignore_intervals: bool = True,
) -> bool:
    """
    Check if the annotation is a subsequence of another annotation.

    Args:
        annotation: The annotation to search for
        other: The annotation to search in
        ignore_mods: If True, only compare sequences (ignore modifications)
        ignore_intervals: If True, ignore interval modifications

    Returns:
        True if annotation is a subsequence of other
    """
    if annotation.sequence not in other.sequence:
        return False

    # Loop over all starting indexes where the sequence is found. Use a zero-width
    # lookahead so overlapping occurrences (e.g. 'II' in 'IIIII') are all found;
    # a plain re.finditer would advance past each match and miss overlaps.
    for match in re.finditer(f"(?={re.escape(annotation.sequence)})", other.sequence):
        start = match.start()

        # Check if all modifications are also a subsequence
        sliced_annot = other.slice(start, start + len(annotation.sequence), inplace=False)

        if ignore_mods and sliced_annot.sequence == annotation.sequence:
            return True

        if ignore_intervals:
            sliced_clean = _remove_intervals_copy(sliced_annot)
            annot_clean = _remove_intervals_copy(annotation)
            if sliced_clean == annot_clean:
                return True

        if sliced_annot == annotation:
            return True

    return False


def find_indices(
    annotation: ProFormaAnnotation,
    other: ProFormaAnnotation,
    ignore_mods: bool = False,
    ignore_intervals: bool = True,
) -> list[int]:
    """
    Find all occurrences of the annotation in another annotation.

    Args:
        annotation: The annotation to search for
        other: The annotation to search in
        ignore_mods: If True, only compare sequences
        ignore_intervals: If True, ignore interval modifications

    Returns:
        List of starting indices where annotation occurs in other

    Raises:
        TypeError: If other is not a ProFormaAnnotation
    """
    if not isinstance(other, type(annotation)):
        raise TypeError(f"other must be a {type(annotation).__name__}, got {type(other).__name__}")

    if not annotation.sequence or not other.sequence:
        return []

    indices: list[int] = []
    for match in re.finditer(f"(?={re.escape(annotation.sequence)})", other.sequence):
        start = match.start()
        sliced_other = other.slice(start, start + len(annotation.sequence), inplace=False)

        if is_subsequence(annotation, sliced_other, ignore_mods, ignore_intervals):
            indices.append(start)

    return indices


def _remove_intervals_copy(annotation: ProFormaAnnotation) -> ProFormaAnnotation:
    """Create a copy of annotation with intervals removed"""
    copy_annot = annotation.copy()
    copy_annot.intervals = None
    return copy_annot

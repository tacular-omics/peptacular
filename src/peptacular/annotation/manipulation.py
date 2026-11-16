from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation

__all__ = [
    "condense_mods_to_intervals",
    "coverage",
    "percent_coverage",
    "modification_coverage",
    "condense_static_mods",
    "condense_to_peptidoform",
    "count_residues",
    "percent_residues",
    "is_subsequence",
    "find_indices",
]


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
    subsequences: Iterable[ProFormaAnnotation],
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> list[int]:
    cov_arr = [0] * len(annotation)

    for sub_annots in subsequences:
        sub_len = len(sub_annots)
        has_ambig = False
        if ignore_ambiguity is False and sub_annots.intervals:
            for interval in sub_annots.intervals:
                if interval.ambiguous:
                    has_ambig = True
                    break

        if not has_ambig:
            for subsequence_index in find_indices(sub_annots, other=annotation, ignore_mods=ignore_mods):
                start = subsequence_index
                if accumulate:
                    for i in range(sub_len):
                        cov_arr[start + i] += 1
                else:
                    cov_arr[start : start + sub_len] = [1] * sub_len
        else:
            peptide_cov = [1] * sub_len
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
    if ignore_ambiguity is False and annotation.intervals:
        # Set coverage to 0 for ambiguous positions
        for interval in annotation.intervals:
            if interval.ambiguous:
                for i in range(interval.start, interval.end):
                    cov_arr[i] = 0

    return cov_arr


def percent_coverage(
    annotation: ProFormaAnnotation,
    subsequences: Iterable[ProFormaAnnotation],
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> float:
    cov_arr = coverage(
        annotation,
        subsequences,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )
    if len(cov_arr) == 0:
        return 0.0

    return sum(cov_arr) / len(cov_arr)


def modification_coverage(
    annotation: ProFormaAnnotation,
    subsequences: Iterable[ProFormaAnnotation],
    ignore_ambiguity: bool = False,
    accumulate: bool = False,
) -> dict[int, int]:
    cov_dict = coverage(
        annotation,
        subsequences,
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


def _residue_key(annotation: ProFormaAnnotation) -> str:
    """
    Serialize one split piece as the residue it holds and that residue's own modifications.

    The pieces that :meth:`ProFormaAnnotation.split` returns carry the parent's charge, the
    parent's labile and unlocalised modifications and the parent's names on every piece, so
    serializing a piece as-is yields keys such as ``P/2`` or ``{+79.966}P`` that are not
    residues. Those annotations belong to the peptide, not to the residue, and are dropped
    here. The internal modifications and the terminal modifications stay because they are
    written at the residue's own position, and so does a global isotope label, which applies
    to every residue.

    Args:
        annotation: One piece from :meth:`ProFormaAnnotation.split`

    Returns:
        The residue with its own modifications, for example ``M[Oxidation]``
    """

    residue = annotation.copy()
    residue.set_charge(None, inplace=True)
    residue.set_labile_mods(None, inplace=True)
    residue.set_unknown_mods(None, inplace=True)
    residue.set_compound_name(None, inplace=True)
    residue.set_ion_name(None, inplace=True)
    residue.set_peptide_name(None, inplace=True)
    return residue.serialize()


def count_residues(annotation: ProFormaAnnotation, include_mods: bool = True) -> dict[str, int]:
    """
    Count the occurrences of each residue in the sequence.

    With ``include_mods`` the keys are a residue plus the modifications localised to it,
    including any N- or C-terminal modification on that residue, so ``PEPM[Oxidation]TIDE``
    holds ``M[Oxidation]``. Annotations that describe the peptide as a whole (charge, labile
    and unlocalised modifications, and names) are not part of a key, so ``PEPTIDE/2`` counts
    the same as ``PEPTIDE``. A global isotope label is the exception and stays in the key,
    because it applies to every residue: ``<13C>PEPTIDE`` holds ``<13C>P``.

    Intervals are counted residue by residue, so the counts always sum to the sequence length:
    ``PEP(?TI)DE`` holds ``T`` and ``I`` and no ``(?TI)`` key. An interval that carries no
    modification adds nothing to the keys. An interval that carries modifications is written
    onto each of its residues as a one-residue interval with the same modifications, so
    ``PEP(?TI)[Phospho]DE`` holds ``(?T)[Phospho]`` and ``(?I)[Phospho]`` and
    ``PEP(TI)[Phospho]DE`` holds ``(T)[Phospho]`` and ``(I)[Phospho]``. That is the form
    ``serialize`` already uses for an interval, and it keeps the modification apart from a
    localised one such as ``T[Phospho]``, because the interval only says the modification is
    somewhere on its residues. ``include_mods=False`` counts plain residue letters.

    Args:
        annotation: The ProFormaAnnotation to analyze
        include_mods: Whether to include modifications in the count

    Returns:
        Dictionary mapping residue strings to their counts
    """

    if not include_mods:
        return dict(Counter(annotation.sequence))

    # split keeps an interval whole as one piece, which is right for split but not for a count,
    # so take the intervals off first, split the plain residues, then put each modified
    # interval back on the residues it spans.
    residues = annotation.copy()
    intervals = residues.pop_intervals(inplace=True)
    keys = [_residue_key(a) for a in residues.split()]
    for interval in intervals:
        if not interval.has_mods:
            continue
        prefix = "(?" if interval.ambiguous else "("
        for i in range(interval.start, interval.end):
            keys[i] = f"{prefix}{keys[i]}){interval.mod_str}"
    return dict(Counter(keys))


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

    can_skip_slicing = ignore_mods or (
        not annotation.has_mods() and not other.has_mods() and (ignore_intervals or (not annotation.intervals and not other.intervals))
    )
    if can_skip_slicing:
        return True

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

    can_skip_slicing = ignore_mods or (
        not annotation.has_mods() and not other.has_mods() and (ignore_intervals or (not annotation.intervals and not other.intervals))
    )

    if can_skip_slicing:
        indices: list[int] = []
        seq = other.sequence
        subseq = annotation.sequence
        idx = seq.find(subseq)
        while idx != -1:
            indices.append(idx)
            idx = seq.find(subseq, idx + 1)
        return indices

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

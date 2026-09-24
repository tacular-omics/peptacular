from collections.abc import Iterable

from ..annotation import (
    ProFormaAnnotation,
)
from .util import get_annotation_input

__all__ = [
    "is_subsequence",
    "find_subsequence_indices",
    "coverage",
    "percent_coverage",
    "modification_coverage",
]


def is_subsequence(subsequence: str | ProFormaAnnotation, sequence: str | ProFormaAnnotation, *, order: bool = True, ignore_mods: bool = False) -> bool:
    """
    Checks if the input subsequence is a subsequence of the input sequence. If order is True, the subsequence must be in
    the same order as in the sequence. If order is False, the subsequence can be in any order.

    .. code-block:: python

        >>> is_subsequence('PEP', 'PEPTIDE')
        True

        >>> is_subsequence('PET', 'PEPTIDE', order=False)
        True

        >>> is_subsequence('PET', 'PEPTIDE', order=True)
        False

        >>> is_subsequence('<13C>PEP', '<13C>PEPTIDE', order=True)
        True

        >>> is_subsequence('<13C>PEP[1.0]', '<13C>PEP[1.0]TIDE', order=True)
        True

        >>> is_subsequence('<13C>PEP', '<13C>PEP[1.0]TIDE', order=True)
        False

    """
    from .basic import count_residues

    if order is True:
        return len(find_subsequence_indices(sequence, subsequence, ignore_mods=ignore_mods)) != 0

    subsequence_counts = count_residues(subsequence, include_mods=not ignore_mods)
    sequence_counts = count_residues(sequence, include_mods=not ignore_mods)
    return all(count <= sequence_counts.get(aa, 0) for aa, count in subsequence_counts.items())


def find_subsequence_indices(sequence: str | ProFormaAnnotation, subsequence: str | ProFormaAnnotation, *, ignore_mods: bool = False) -> list[int]:
    """
    Retrieves all starting indexes of a given subsequence within a sequence.

    .. code-block:: python

        # Find the starting indexes of a subsequence
        >>> find_subsequence_indices("PEPTIDE", "PEP")
        [0]
        >>> find_subsequence_indices("PEPTIDE", "EPT")
        [1]
        >>> find_subsequence_indices("PEPTIDE", "E")
        [1, 6]

        # By default the function will not ignore modifications
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "PEP")
        []
        >>> find_subsequence_indices("<13C>PEPTIDE", "PEP")
        []
        >>> find_subsequence_indices("<13C>PEP[1][Phospho]TIDE", "<13C>PEP[Phospho][1]")
        [0]
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "[Acetyl]-PEP")
        [0]
        >>> find_subsequence_indices("PEPTIDE", "[Acetyl]-PEP")
        []
        >>> find_subsequence_indices("PEPTIDE[1.0]", "IDE")
        []
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE-[Amide]", "[Acetyl]-PEPTIDE-[Amide]")
        [0]

        # If ignore_mods is set to True, the function will ignore modifications
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "PEP", ignore_mods=True)
        [0]
        >>> find_subsequence_indices("PEPTIDE", "[Acetyl]-PEP", ignore_mods=True)
        [0]

    """
    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequence_annot = get_annotation_input(subsequence, copy=False)
    return subsequence_annot.find_indices(other=sequence_annot, ignore_mods=ignore_mods)


def coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: Iterable[str | ProFormaAnnotation],
    *,
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> list[int]:
    """
    Calculate the sequence coverage given a list of subsequecnes.

    The coverage is represented as a binary list where each position in the protein sequence is marked as 1 if it
    is covered by at least one peptide and 0 otherwise.

    .. code-block:: python

        >>> coverage("PEPTIDE", ["PEP"])
        [1, 1, 1, 0, 0, 0, 0]

        >>> coverage("PEPTIDE", ["PEP", "EPT"])
        [1, 1, 1, 1, 0, 0, 0]

        # If accumulate is True, overlapping indecies will be accumulated
        >>> coverage("PEPTIDE", ["PEP", "EPT"], accumulate=True)
        [1, 2, 2, 1, 0, 0, 0]

        # By default ambiguity does not add to coverage
        >>> coverage("PEPTIDE", ["P(?EP)"])
        [1, 0, 0, 0, 0, 0, 0]

    """

    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequences_annot = [get_annotation_input(subseq, copy=False) for subseq in subsequences]

    return sequence_annot.coverage(
        subsequences_annot,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )


def percent_coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: Iterable[str | ProFormaAnnotation],
    *,
    ignore_mods: bool = False,
    accumulate: bool = False,
    ignore_ambiguity: bool = False,
) -> float:
    """
    Calculates the coverage given a list of subsequences.

    if accumulate is True, overlapping indecies will be accumulated.
    if ignore_mods is True, modifications will be ignored when calculating coverage (only the amino acid sequence will be considered).
    if ignore_ambiguity is True, ambiguous regions will not be counted towards coverage.

    .. code-block:: python

        >>> round(percent_coverage("PEPTIDE", ["PEP"]), 3)
        0.429

        # ambiguity does not add to coverage by default
        >>> round(percent_coverage("PEPTIDE", ["P(?EP)"]), 3)
        0.143

        >>> round(percent_coverage("PEPTIDE", ["PEP", "EPT"]), 3)
        0.571

        >>> round(percent_coverage("PEPTIDE", ["PEPTIDE", "PEPTIDE"], accumulate=True), 3)
        2.0

    """

    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequences_annot = [get_annotation_input(subseq, copy=False) for subseq in subsequences]

    return sequence_annot.percent_coverage(
        subsequences_annot,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )


def modification_coverage(sequence: str | ProFormaAnnotation, subsequences: list[str | ProFormaAnnotation], *, accumulate: bool = False) -> dict[int, int]:
    """
    Calculate the modification coverage given a list of subsequences.

    This function identifies which modifications in the main sequence are covered by
    subsequences. It returns a dictionary where each key is a modification position/type
    (matching the format from get_mods) and each value is the number of subsequences
    that cover that modification.

    .. code-block:: python

        >>> modification_coverage("PEPTIDE[Phospho]", ["TIDE"])
        {6: 0}

        >>> modification_coverage("PEPTIDE[Phospho]", ["TIDE[Phospho]"])
        {6: 1}

        >>> modification_coverage("PEP[Phospho]TIDE[Methyl]", ["PEP[Phospho]", "TIDE[Methyl]"])
        {2: 1, 6: 1}

        >>> modification_coverage("PEP[Phospho]TIDE[Methyl]", ["PEP[Phospho]", "TIDE[Methyl]", "PEP[Phospho]"], accumulate=True)
        {2: 2, 6: 1}

    """
    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequence_annots = [get_annotation_input(subseq, copy=False) for subseq in subsequences]

    return sequence_annot.modification_coverage(annotations=subsequence_annots, accumulate=accumulate)

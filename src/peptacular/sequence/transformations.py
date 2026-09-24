from collections.abc import Callable, Sequence
from typing import Any, overload

from ..annotation import (
    ProFormaAnnotation,
)
from ..constants import ParallelMethod, ParallelMethodLiteral
from ..spans import Span
from .parallel import parallel_apply_internal
from .util import get_annotation_input

__all__ = [
    "reverse",
    "shuffle",
    "shift",
    "span_to_sequence",
    "split",
    "sort",
    "join",
]


def _reverse_single(
    sequence: str | ProFormaAnnotation,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
) -> str:
    return get_annotation_input(sequence=sequence, copy=True).reverse(inplace=True, keep_nterm=keep_nterm, keep_cterm=keep_cterm).serialize()


@overload
def reverse(
    sequence: str | ProFormaAnnotation,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def reverse(
    sequence: Sequence[str | ProFormaAnnotation],
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def reverse(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Reverses the sequence, while preserving the position of any modifications.

    swap_terms: If True, the N- and C-terminal modifications will be swapped.
    keep_nterm: Number of N-terminal residues to keep unchanged. Default is 0.
    keep_cterm: Number of C-terminal residues to keep unchanged. Default is 0.

    .. code-block:: python

        # Single sequence
        >>> reverse('PEPTIDE')
        'EDITPEP'

        # Keep first 2 residues unchanged
        >>> reverse('PEPTIDE', keep_nterm=2)
        'PEEDITP'

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _reverse_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
        )
    else:
        return _reverse_single(sequence, keep_nterm, keep_cterm)


def _shuffle_single(
    sequence: str | ProFormaAnnotation,
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
) -> str:
    """Internal function for shuffling a single sequence."""
    return get_annotation_input(sequence=sequence, copy=True).shuffle(seed=seed, keep_nterm=keep_nterm, keep_cterm=keep_cterm, inplace=True).serialize()


@overload
def shuffle(
    sequence: str | ProFormaAnnotation,
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def shuffle(
    sequence: Sequence[str | ProFormaAnnotation],
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def shuffle(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Shuffles the sequence, while preserving the position of any modifications.

    keep_nterm: Number of N-terminal residues to keep unchanged. Default is 0.
    keep_cterm: Number of C-terminal residues to keep unchanged. Default is 0.

    .. code-block:: python

        # Single sequence
        >>> shuffle('PEPTIDE', seed=0)
        'IPEPDTE'

        # Keep first 2 residues unchanged
        >>> shuffle('PEPTIDE', seed=0, keep_nterm=2)
        'PEITPED'

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _shuffle_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            seed=seed,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
        )
    else:
        return _shuffle_single(sequence, seed, keep_nterm, keep_cterm)


def _shift_single(
    sequence: str | ProFormaAnnotation,
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
) -> str:
    """Internal function for shifting a single sequence."""
    return get_annotation_input(sequence=sequence, copy=True).shift(n=n, keep_nterm=keep_nterm, keep_cterm=keep_cterm, inplace=True).serialize()


@overload
def shift(
    sequence: str | ProFormaAnnotation,
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def shift(
    sequence: Sequence[str | ProFormaAnnotation],
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def shift(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    keep_nterm: Number of N-terminal residues to keep unchanged. Default is 0.
    keep_cterm: Number of C-terminal residues to keep unchanged. Default is 0.

    .. code-block:: python

        # Single sequence
        >>> shift('PEPTIDE', 2)
        'PTIDEPE'

        # Keep first 2 residues unchanged
        >>> shift('PEPTIDE', 2, keep_nterm=2)
        'PEIDEPT'

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _shift_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            n=n,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
        )
    else:
        return _shift_single(sequence, n, keep_nterm, keep_cterm)


def _span_to_sequence_single(
    sequence: str | ProFormaAnnotation,
    span: tuple[int, int, int],
) -> str:
    return get_annotation_input(sequence=sequence, copy=True).slice(span[0], span[1], inplace=True).serialize()


@overload
def span_to_sequence(
    sequence: str | ProFormaAnnotation,
    span: tuple[int, int, int],
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def span_to_sequence(
    sequence: Sequence[str | ProFormaAnnotation],
    span: tuple[int, int, int],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def span_to_sequence(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    span: tuple[int, int, int] | Span,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Extracts a subsequence from the input sequence based on the provided span.

    The span is defined as a tuple of three integers: (start, end, step).

    .. code-block:: python

        # Single sequence
        >>> span_to_sequence('PEPTIDE', (0, 4, 0))
        'PEPT'

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _span_to_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            span=span,
        )
    else:
        return _span_to_sequence_single(sequence, span)


def _split_single(
    sequence: str | ProFormaAnnotation,
) -> list[str]:
    return [a.serialize() for a in get_annotation_input(sequence=sequence, copy=True).split()]


@overload
def split(
    sequence: str | ProFormaAnnotation,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def split(
    sequence: Sequence[str | ProFormaAnnotation],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def split(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Splits sequence into a list of amino acids, preserving modifications.

    .. code-block:: python

        # Single sequence
        >>> split('PEPTIDE')
        ['P', 'E', 'P', 'T', 'I', 'D', 'E']

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _split_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _split_single(sequence)


def _sort_single(
    sequence: str | ProFormaAnnotation,
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
) -> str:
    return get_annotation_input(sequence=sequence, copy=True).sort(inplace=True, key=key, reverse=reverse).serialize()


@overload
def sort(
    sequence: str | ProFormaAnnotation,
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def sort(
    sequence: Sequence[str | ProFormaAnnotation],
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def sort(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Sorts the input sequence using the provided sort function. Terminal sequences are kept in place.

    key: A function that serves as a key for the sort comparison. Default is None.
    reverse: If True, the sorted sequence is reversed (descending order). Default is False.

    .. code-block:: python

        # Single sequence
        >>> sort('PEPTIDE')
        'DEEIPPT'

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _sort_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            key=key,
            reverse=reverse,
        )
    else:
        return _sort_single(sequence, key, reverse)


def _join_single(
    annotations: Sequence[ProFormaAnnotation | str],
) -> str:
    annotations = [get_annotation_input(a, copy=False) for a in annotations]
    return ProFormaAnnotation.join(annotations).serialize()


@overload
def join(
    annotations: Sequence[ProFormaAnnotation | str],
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def join(
    annotations: Sequence[Sequence[ProFormaAnnotation | str]],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def join(
    annotations: Sequence[ProFormaAnnotation | str] | Sequence[Sequence[ProFormaAnnotation | str]],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Joins a list of annotations into a single annotation.

    .. code-block:: python

        # Single list of annotations
        >>> join(['PEPTIDE', 'MODIFIED'])
        'PEPTIDEMODIFIED'

    """
    if isinstance(annotations, Sequence) and isinstance(annotations[0], Sequence) and not isinstance(annotations[0], (str, ProFormaAnnotation)):
        return parallel_apply_internal(
            _join_single,
            annotations,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _join_single(annotations)  # type: ignore

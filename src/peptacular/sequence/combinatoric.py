from collections.abc import Sequence
from typing import overload

from ..annotation import ProFormaAnnotation
from ..constants import ParallelMethod, ParallelMethodLiteral
from .parallel import parallel_apply_internal
from .util import get_annotation_input

__all__ = [
    "permutations",
    "product",
    "combinations",
    "combinations_with_replacement",
]


def _permutations_single(
    sequence: str | ProFormaAnnotation,
    size: int | None = None,
) -> list[str]:
    annotation = get_annotation_input(sequence, copy=False)
    return [a.serialize() for a in annotation.permutations(size=size)]


@overload
def permutations(
    sequence: str | ProFormaAnnotation,
    size: int | None = None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def permutations(
    sequence: Sequence[str | ProFormaAnnotation],
    size: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def permutations(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    size: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Generates all permutations of the input sequence. Terminal sequence are kept in place.

    The size of the permutations. If None, uses the length of the sequence.

    .. code-block:: python

        >>> permutations('PET')
        ['PET', 'PTE', 'EPT', 'ETP', 'TPE', 'TEP']

        >>> permutations('[3]-PET-[1]')
        ['[3]-PET-[1]', '[3]-PTE-[1]', '[3]-EPT-[1]', '[3]-ETP-[1]', '[3]-TPE-[1]', '[3]-TEP-[1]']

        >>> permutations('PE[3.14]T')
        ['PE[3.14]T', 'PTE[3.14]', 'E[3.14]PT', 'E[3.14]TP', 'TPE[3.14]', 'TE[3.14]P']

        >>> permutations('<13C>PET')
        ['<13C>PET', '<13C>PTE', '<13C>EPT', '<13C>ETP', '<13C>TPE', '<13C>TEP']

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _permutations_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            size=size,
        )
    else:
        return _permutations_single(
            sequence=sequence,
            size=size,
        )


def _product_single(
    sequence: str | ProFormaAnnotation,
    repeat: int | None,
) -> list[str]:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return [a.serialize() for a in annotation.product(repeat=repeat)]


@overload
def product(
    sequence: str | ProFormaAnnotation,
    repeat: int | None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def product(
    sequence: Sequence[str | ProFormaAnnotation],
    repeat: int | None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def product(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    repeat: int | None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Generates all cartesian products of the input sequence of a given size. Terminal sequence are kept in place.

    The size of the combinations to be generated. If None, uses the length of the sequence.

    .. code-block:: python

        >>> product('PET', 2)
        ['PP', 'PE', 'PT', 'EP', 'EE', 'ET', 'TP', 'TE', 'TT']

        >>> product('[3]-PET-[1]', 2)[:5]
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EP-[1]', '[3]-EE-[1]']

        >>> product('PE[3.14]T', 2)
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]P', 'E[3.14]E[3.14]', 'E[3.14]T', 'TP', 'TE[3.14]', 'TT']

        >>> product('<13C>PET', 2)
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EP', '<13C>EE', '<13C>ET', '<13C>TP', '<13C>TE', '<13C>TT']

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _product_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            repeat=repeat,
        )
    else:
        return _product_single(
            sequence=sequence,
            repeat=repeat,
        )


def _combinations_single(
    sequence: str | ProFormaAnnotation,
    size: int | None,
) -> list[str]:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return [a.serialize() for a in annotation.combinations(r=size)]


@overload
def combinations(
    sequence: str | ProFormaAnnotation,
    size: int | None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def combinations(
    sequence: Sequence[str | ProFormaAnnotation],
    size: int | None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def combinations(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    size: int | None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Generates all combinations of the input sequence of a given size. Terminal sequence are kept in place.

    size: The size of the combinations to be generated. If None, uses the length of the sequence.

    .. code-block:: python

        >>> combinations('PET', 2)
        ['PE', 'PT', 'ET']

        >>> combinations('[3]-PET-[1]', 2)
        ['[3]-PE-[1]', '[3]-PT-[1]', '[3]-ET-[1]']

        >>> combinations('PE[3.14]T', 2)
        ['PE[3.14]', 'PT', 'E[3.14]T']

        >>> combinations('<13C>PET', 2)
        ['<13C>PE', '<13C>PT', '<13C>ET']

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _combinations_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            size=size,
        )
    else:
        return _combinations_single(
            sequence=sequence,
            size=size,
        )


def _combinations_with_replacement_single(
    sequence: str | ProFormaAnnotation,
    size: int | None,
) -> list[str]:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return [a.serialize() for a in annotation.combinations_with_replacement(r=size)]


@overload
def combinations_with_replacement(
    sequence: str | ProFormaAnnotation,
    size: int | None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def combinations_with_replacement(
    sequence: Sequence[str | ProFormaAnnotation],
    size: int | None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def combinations_with_replacement(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    size: int | None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Generates all combinations with replacement of the input sequence of a given size. Terminal sequence are kept
    in place.

    size: The size of the combinations to be generated. If None, uses the length of the sequence.

    .. code-block:: python

        >>> combinations_with_replacement('PET', 2)
        ['PP', 'PE', 'PT', 'EE', 'ET', 'TT']

        >>> combinations_with_replacement('[3]-PET-[1]', 2)
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EE-[1]', '[3]-ET-[1]', '[3]-TT-[1]']

        >>> combinations_with_replacement('PE[3.14]T', 2)
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]E[3.14]', 'E[3.14]T', 'TT']

        >>> combinations_with_replacement('<13C>PET', 2)
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EE', '<13C>ET', '<13C>TT']

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _combinations_with_replacement_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            size=size,
        )
    else:
        return _combinations_with_replacement_single(
            sequence=sequence,
            size=size,
        )

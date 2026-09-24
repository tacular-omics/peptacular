from collections.abc import Sequence
from typing import overload

from ..annotation import ProFormaAnnotation
from ..constants import ParallelMethod, ParallelMethodLiteral
from ..digestion.core import generate_regex
from ..spans import Span
from .parallel import parallel_apply_internal
from .util import get_annotation_input

__all__ = [
    "left_semi_digest",
    "right_semi_digest",
    "semi_digest",
    "nonspecific_digest",
    "cleavage_sites",
    "simple_cleavage_sites",
    "digest",
    "simple_digest",
]


def _left_semi_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[tuple[str, Span]]:
    annot = get_annotation_input(sequence, copy=False)
    return [
        (annot[span].serialize(), span)
        for span in annot.left_semi_spans(
            min_len=min_len,
            max_len=max_len,
        )
    ]


@overload
def left_semi_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]]: ...


@overload
def left_semi_digest(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[tuple[str, Span]]]: ...


def left_semi_digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]] | list[list[tuple[str, Span]]]:
    """Semi-enzymatic sequences that keep the N-terminus of ``sequence`` (every prefix shorter than the full sequence, within the length limits).

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param min_len: Minimum length.
    :param max_len: Maximum length.
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: ``(sequence, Span)`` tuples, or a list of such lists for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _left_semi_digest,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _left_semi_digest(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _right_semi_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[tuple[str, Span]]:
    annot = get_annotation_input(sequence, copy=False)
    return [
        (annot[span].serialize(), span)
        for span in annot.right_semi_spans(
            min_len=min_len,
            max_len=max_len,
        )
    ]


@overload
def right_semi_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]]: ...


@overload
def right_semi_digest(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[tuple[str, Span]]]: ...


def right_semi_digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]] | list[list[tuple[str, Span]]]:
    """Semi-enzymatic sequences that keep the C-terminus of ``sequence`` (every suffix shorter than the full sequence, within the length limits).

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param min_len: Minimum length.
    :param max_len: Maximum length.
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: ``(sequence, Span)`` tuples, or a list of such lists for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _right_semi_digest,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _right_semi_digest(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _semi_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[tuple[str, Span]]:
    annot = get_annotation_input(sequence, copy=False)
    return [
        (annot[span].serialize(), span)
        for span in annot.semi_spans(
            min_len=min_len,
            max_len=max_len,
        )
    ]


@overload
def semi_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]]: ...


@overload
def semi_digest(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[tuple[str, Span]]]: ...


def semi_digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]] | list[list[tuple[str, Span]]]:
    """
    Builds all semi-enzymatic sequences from the given input `sequence`.
    Equivalent to combining left and right semi-enzymatic sequences.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _semi_digest,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _semi_digest(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _nonspecific_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[tuple[str, Span]]:
    annot = get_annotation_input(sequence, copy=False)
    return [
        (annot[span].serialize(), span)
        for span in annot.nonspecific_spans(
            min_len=min_len,
            max_len=max_len,
        )
    ]


@overload
def nonspecific_digest(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]]: ...


@overload
def nonspecific_digest(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[tuple[str, Span]]]: ...


def nonspecific_digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]] | list[list[tuple[str, Span]]]:
    """
    Builds all non-enzymatic sequences from the given input `sequence`.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _nonspecific_digest,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _nonspecific_digest(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _cleavage_sites(sequence: str | ProFormaAnnotation, enzyme_regex: str) -> list[int]:
    return list(get_annotation_input(sequence, copy=False).cleavage_sites(enzyme=enzyme_regex))


@overload
def cleavage_sites(
    sequence: str | ProFormaAnnotation,
    enzyme_regex: str,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[int]: ...


@overload
def cleavage_sites(
    sequence: Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[int]]: ...


def cleavage_sites(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[int] | list[list[int]]:
    """
    Return positions where cleavage occurs in input `sequence` based on the provided enzyme regex.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _cleavage_sites,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            enzyme_regex=enzyme_regex,
        )
    else:
        return _cleavage_sites(
            sequence=sequence,
            enzyme_regex=enzyme_regex,
        )


def _simple_cleavage_sites(
    sequence: str | ProFormaAnnotation,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
) -> list[int]:
    enzyme_regex = generate_regex(
        cleave_on=cleave_on,
        restrict_before=restrict_before,
        restrict_after=restrict_after,
        cterminal=cterminal,
    )
    return list(get_annotation_input(sequence, copy=False).cleavage_sites(enzyme=enzyme_regex))


@overload
def simple_cleavage_sites(
    sequence: str | ProFormaAnnotation,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[int]: ...


@overload
def simple_cleavage_sites(
    sequence: Sequence[str | ProFormaAnnotation],
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[int]]: ...


def simple_cleavage_sites(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[int] | list[list[int]]:
    """
    Get cleavage sites using simple amino acid rules.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_cleavage_sites,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        )
    else:
        return _simple_cleavage_sites(
            sequence=sequence,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        )


def _digest(
    sequence: str | ProFormaAnnotation,
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[tuple[str, Span]]:
    annot = get_annotation_input(sequence, copy=False)
    return [
        (annot[span].serialize(), span)
        for span in annot.digest(
            enzyme=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )
    ]


@overload
def digest(
    sequence: str | ProFormaAnnotation,
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]]: ...


@overload
def digest(
    sequence: Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[tuple[str, Span]]]: ...


def digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]] | list[list[tuple[str, Span]]]:
    """
    Returns digested sequences using a regular expression to define cleavage sites.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _digest,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _digest(
            sequence=sequence,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )


def _digest_single(
    sequence: str | ProFormaAnnotation,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[tuple[str, Span]]:
    annot = get_annotation_input(sequence, copy=False)
    return [
        (annot[span].serialize(), span)
        for span in annot.simple_digest(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )
    ]


@overload
def simple_digest(
    sequence: str | ProFormaAnnotation,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]]: ...


@overload
def simple_digest(
    sequence: Sequence[str | ProFormaAnnotation],
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[tuple[str, Span]]]: ...


def simple_digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[tuple[str, Span]] | list[list[tuple[str, Span]]]:
    """
    Returns digested sequences using amino acid specifications with optional restrictions.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _digest_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _digest_single(
            sequence=sequence,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )

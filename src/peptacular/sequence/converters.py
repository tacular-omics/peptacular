from typing import overload

from ..annotation import ProFormaAnnotation
from ..constants import ParallelMethod, ParallelMethodLiteral
from .parallel import parallel_apply_internal

__all__ = [
    "convert_ip2_sequence",
    "convert_diann_sequence",
    "convert_casanovo_sequence",
]


def _convert_ip2_sequence_single(sequence: str) -> str:
    return ProFormaAnnotation.from_ip2_sequence(sequence).serialize()


@overload
def convert_ip2_sequence(
    sequence: str,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def convert_ip2_sequence(
    sequence: list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def convert_ip2_sequence(
    sequence: str | list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Converts a IP2-Like sequence to a proforma2.0 compatible sequence.

    .. code-block:: python

        >>> convert_ip2_sequence('K.PEP(phospho)TIDE.K')
        'PEP[phospho]TIDE'

    """
    if isinstance(sequence, list):
        return parallel_apply_internal(
            _convert_ip2_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _convert_ip2_sequence_single(sequence)


def _convert_diann_sequence_single(sequence: str) -> str:
    return ProFormaAnnotation.from_diann(sequence).serialize()


@overload
def convert_diann_sequence(
    sequence: str,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def convert_diann_sequence(
    sequence: list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def convert_diann_sequence(
    sequence: str | list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Converts a DIANN-Like sequence to a proforma2.0 compatible sequence.

    .. code-block:: python

        >>> convert_diann_sequence('_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]_[Methyl]')
        '[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]-[Methyl]'

    """
    if isinstance(sequence, list):
        return parallel_apply_internal(
            _convert_diann_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _convert_diann_sequence_single(sequence)


def _convert_casanovo_sequence_single(sequence: str) -> str:
    return ProFormaAnnotation.from_casanovo(sequence).serialize()


@overload
def convert_casanovo_sequence(
    sequence: str,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str: ...


@overload
def convert_casanovo_sequence(
    sequence: list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[str]: ...


def convert_casanovo_sequence(
    sequence: str | list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Converts a Casanovo sequence with modifications to a proforma2.0 compatible sequence.

    :return: Proforma2.0 compatable sequence or list of sequences.
    :rtype: str | list[str]

    .. code-block:: python

        >>> convert_casanovo_sequence('+43.006P+100EPTIDE')
        '[+43.006]-P[+100]EPTIDE'

    """
    if isinstance(sequence, list):
        return parallel_apply_internal(
            _convert_casanovo_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _convert_casanovo_sequence_single(sequence)

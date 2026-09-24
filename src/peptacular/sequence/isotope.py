from collections.abc import Sequence
from typing import overload

from tacular import IonType

from ..annotation import ProFormaAnnotation
from ..annotation.annotation import (
    CHARGE_TYPE,
    CUSTOM_LOSS_TYPE,
    ION_TYPE,
    ISOTOPE_TYPE,
)
from ..constants import ParallelMethod, ParallelMethodLiteral
from ..isotope import IsotopicData
from .parallel import parallel_apply_internal
from .util import HasSequence, get_annotation_input

__all__ = [
    "isotopic_distribution",
    "estimate_isotopic_distribution",
]


def _isotopic_distribution_single(
    annotation: str | ProFormaAnnotation | HasSequence,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,  # based on the most abundant peak
) -> list[IsotopicData]:
    return get_annotation_input(annotation).isotopic_distribution(
        ion_type=ion_type,
        charge=charge,
        isotopes=isotopes,
        deltas=deltas,
        max_isotopes=max_isotopes,
        min_abundance_threshold=min_abundance_threshold,
    )


@overload
def isotopic_distribution(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[IsotopicData]]: ...


@overload
def isotopic_distribution(
    sequence: str | ProFormaAnnotation | HasSequence,
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[IsotopicData]: ...


def isotopic_distribution(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[IsotopicData] | list[list[IsotopicData]]:
    """Exact isotopic distribution of a peptide ion from its elemental composition.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param ion_type: Ion type whose composition is used; defaults to the precursor.
    :param charge: Charge state or charge carriers; ``None`` uses the sequence's own charge.
    :param isotopes: Isotope labels to apply.
    :param deltas: Extra mass or formula deltas to apply.
    :param max_isotopes: Keep at most this many peaks.
    :param min_abundance_threshold: Drop peaks below this relative abundance.
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A list of ``IsotopicData`` peaks, or a list of such lists for list input.
    :raises CompositionError: The ion's composition is not available (e.g. a mass-only modification).
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, (str, ProFormaAnnotation)):
        return parallel_apply_internal(
            _isotopic_distribution_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
        )
    else:
        return _isotopic_distribution_single(
            sequence,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
        )


def _estimate_isotopic_distribution_single(
    annotation: str | ProFormaAnnotation | HasSequence,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
) -> list[IsotopicData]:
    return get_annotation_input(annotation).estimate_isotopic_distribution(
        ion_type=ion_type,
        charge=charge,
        isotopes=isotopes,
        deltas=deltas,
        max_isotopes=max_isotopes,
        min_abundance_threshold=min_abundance_threshold,
    )


@overload
def estimate_isotopic_distribution(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[IsotopicData]]: ...


@overload
def estimate_isotopic_distribution(
    sequence: str | ProFormaAnnotation | HasSequence,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[IsotopicData]: ...


def estimate_isotopic_distribution(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[IsotopicData] | list[list[IsotopicData]]:
    if isinstance(sequence, Sequence) and not isinstance(sequence, (str, ProFormaAnnotation)):
        return parallel_apply_internal(
            _estimate_isotopic_distribution_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
        )
    else:
        return _estimate_isotopic_distribution_single(
            sequence,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
        )

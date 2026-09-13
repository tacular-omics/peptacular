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
from ..constants import parallelMethod, parallelMethodLiteral
from ..isotope import IsotopicData
from .parallel import parallel_apply_internal
from .util import get_annotation_input


def _isotopic_distribution_single(
    annotation: str | ProFormaAnnotation,
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
    annotations: Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[list[IsotopicData]]: ...


@overload
def isotopic_distribution(
    annotations: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[IsotopicData]: ...


def isotopic_distribution(
    annotations: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[IsotopicData] | list[list[IsotopicData]]:
    if isinstance(annotations, Sequence) and not isinstance(annotations, (str, ProFormaAnnotation)):
        return parallel_apply_internal(
            _isotopic_distribution_single,
            annotations,
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
            annotations,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
        )


def _estimate_isotopic_distribution_single(
    annotation: str | ProFormaAnnotation,
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
    annotations: Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[list[IsotopicData]]: ...


@overload
def estimate_isotopic_distribution(
    annotations: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[IsotopicData]: ...


def estimate_isotopic_distribution(
    annotations: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = 0.001,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[IsotopicData] | list[list[IsotopicData]]:
    if isinstance(annotations, Sequence) and not isinstance(annotations, (str, ProFormaAnnotation)):
        return parallel_apply_internal(
            _estimate_isotopic_distribution_single,
            annotations,
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
            annotations,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
            max_isotopes=max_isotopes,
            min_abundance_threshold=min_abundance_threshold,
        )

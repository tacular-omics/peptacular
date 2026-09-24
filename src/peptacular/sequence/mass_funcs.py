from collections import Counter
from collections.abc import Sequence
from typing import overload

from tacular import ElementInfo, IonType

from ..annotation import ProFormaAnnotation
from ..annotation.annotation import (
    CHARGE_TYPE,
    CUSTOM_LOSS_TYPE,
    ION_TYPE,
    ISOTOPE_TYPE,
)
from ..constants import (
    ParallelMethod,
    ParallelMethodLiteral,
)
from .parallel import parallel_apply_internal
from .util import get_annotation_input

__all__ = [
    "mass",
    "mz",
    "comp",
]


def _mass_single(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE,
    charge: CHARGE_TYPE | None,
    monoisotopic: bool,
    isotopes: ISOTOPE_TYPE | None,
    deltas: CUSTOM_LOSS_TYPE | None,
    calculate_with_composition: bool,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.mass(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        calculate_with_composition=calculate_with_composition,
    )


@overload
def mass(
    sequence: str | ProFormaAnnotation,
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def mass(
    sequence: Sequence[str | ProFormaAnnotation],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def mass(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """
    Calculate the mass of an amino acid 'sequence'.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _mass_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_with_composition=calculate_with_composition,
        )
    else:
        return _mass_single(
            sequence=sequence,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_with_composition=calculate_with_composition,
        )


def _mz_single(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE,
    charge: CHARGE_TYPE | None,
    monoisotopic: bool,
    isotopes: ISOTOPE_TYPE | None,
    deltas: CUSTOM_LOSS_TYPE | None,
    calculate_with_composition: bool,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.mz(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        calculate_with_composition=calculate_with_composition,
    )


@overload
def mz(
    sequence: str | ProFormaAnnotation,
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def mz(
    sequence: Sequence[str | ProFormaAnnotation],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def mz(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """
    Calculate the m/z (mass-to-charge ratio) of an amino acid 'sequence'.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _mz_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_with_composition=calculate_with_composition,
        )
    else:
        return _mz_single(
            sequence=sequence,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_with_composition=calculate_with_composition,
        )


def _comp_single(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
) -> Counter[ElementInfo]:
    annotation = get_annotation_input(sequence=sequence, copy=True)

    return annotation.comp(
        ion_type=ion_type,
        charge=charge,
        isotopes=isotopes,
        deltas=deltas,
    )


@overload
def comp(
    sequence: str | ProFormaAnnotation,
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> Counter[ElementInfo]: ...


@overload
def comp(
    sequence: Sequence[str | ProFormaAnnotation],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[Counter[ElementInfo]]: ...


def comp(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    charge: CHARGE_TYPE | None = None,
    *,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> Counter[ElementInfo] | list[Counter[ElementInfo]]:
    """
    Calculates the elemental composition of a peptide sequence, including modifications.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _comp_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
        )
    else:
        return _comp_single(
            sequence=sequence,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            deltas=deltas,
        )

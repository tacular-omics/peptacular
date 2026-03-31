from collections.abc import Sequence
from typing import overload

from tacular import IonType

from ..annotation import ProFormaAnnotation
from ..annotation.annotation import (
    CHARGE_TYPE,
    CUSTOM_LOSS_TYPE,
    ION_TYPE,
    ISOTOPE_TYPE,
    LOSS_TYPE,
)
from ..annotation.utils import Fragment
from ..constants import parallelMethod, parallelMethodLiteral
from .parallel import parallel_apply_internal
from .util import get_annotation_input

FRAGMENT_MASSES_RETURN = dict[tuple[IonType, int], list[float]]


def _fragment_single(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE] = (),
    calculate_composition: bool = False,
    max_ndeltas: int = 1,
) -> list[Fragment]:
    annotation = get_annotation_input(sequence=sequence, copy=False)

    return annotation.fragment(
        ion_types=ion_types,
        charges=charges,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        neutral_deltas=neutral_deltas,
        calculate_composition=calculate_composition,
        max_ndeltas=max_ndeltas,
    )


@overload
def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] | ION_TYPE = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (None,),
    calculate_composition: bool = False,
    max_ndeltas: int = 1,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[Fragment]: ...


@overload
def fragment(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] | ION_TYPE = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (None,),
    calculate_composition: bool = False,
    max_ndeltas: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[list[Fragment]]: ...


def fragment(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE] = (),
    calculate_composition: bool = False,
    max_ndeltas: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[Fragment] | list[list[Fragment]]:
    """
    Builds fragment ions from a given input sequence or list of sequences.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _fragment_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            max_ndeltas=max_ndeltas,
            calculate_composition=calculate_composition,
        )
    else:
        return _fragment_single(
            sequence=sequence,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            max_ndeltas=max_ndeltas,
            calculate_composition=calculate_composition,
        )


def _frag_single(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_composition: bool = False,
    position: int | tuple[int, int] | None = None,
) -> Fragment:
    annotation = get_annotation_input(sequence=sequence, copy=False)

    return annotation.frag(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        calculate_composition=calculate_composition,
        position=position,
    )


@overload
def frag(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_composition: bool = False,
    position: int | tuple[int, int] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> Fragment: ...


@overload
def frag(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_composition: bool = False,
    position: int | tuple[int, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[Fragment]: ...


def frag(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_composition: bool = False,
    position: int | tuple[int, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> Fragment | list[Fragment]:
    """
    Calculate a single fragment from a sequence or multiple sequences.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _frag_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_type=ion_type,
            charge=charge,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_composition=calculate_composition,
            position=position,
        )
    else:
        return _frag_single(
            sequence=sequence,
            ion_type=ion_type,
            charge=charge,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            calculate_composition=calculate_composition,
            position=position,
        )


def _fast_fragment_single(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[int] | None = None,
    monoisotopic: bool = True,
) -> FRAGMENT_MASSES_RETURN:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.fast_fragment(
        ion_types=ion_types,
        charges=charges,
        monoisotopic=monoisotopic,
    )


@overload
def fast_fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[int] | None = None,
    monoisotopic: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> FRAGMENT_MASSES_RETURN: ...


@overload
def fast_fragment(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[int] | None = None,
    monoisotopic: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[FRAGMENT_MASSES_RETURN]: ...


def fast_fragment(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[int] | None = None,
    monoisotopic: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> FRAGMENT_MASSES_RETURN | list[FRAGMENT_MASSES_RETURN]:
    """Compute fragment ion m/z values for a sequence or list of sequences.

    Uses a fast prefix/suffix-sum approach. Returns a dict mapping
    ``(IonType, charge)`` to a list of m/z values of length equal to the
    sequence length, ordered from fragment position 1 to N. Neutral losses,
    isotope shifts, and custom deltas are not supported.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _fast_fragment_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
        )
    else:
        return _fast_fragment_single(
            sequence=sequence,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
        )

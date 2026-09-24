import os
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, overload

from tacular import IonType

from ..annotation import ProFormaAnnotation
from ..annotation.annotation import (
    CHARGE_TYPE,
    CUSTOM_LOSS_TYPE,
    ION_TYPE,
    ISOTOPE_TYPE,
    LOSS_TYPE,
)
from ..annotation.frag_arrays import FRAGMENT_ARRAY_KEYS
from ..annotation.frag_arrays import fragment_arrays as _fragment_arrays
from ..annotation.utils import Fragment
from ..constants import ParallelMethod, ParallelMethodLiteral
from .parallel import parallel_apply_internal
from .util import HasSequence, get_annotation_input

if TYPE_CHECKING:
    import numpy as np

__all__ = [
    "FRAGMENT_MASSES_RETURN",
    "fragment",
    "frag",
    "fast_fragment",
    "fragment_arrays",
    "FRAGMENT_ARRAY_KEYS",
]

FRAGMENT_MASSES_RETURN = dict[tuple[IonType, int], list[float]]


def _fragment_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (),
    calculate_with_composition: bool = False,
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
        calculate_with_composition=calculate_with_composition,
        max_ndeltas=max_ndeltas,
    )


@overload
def fragment(
    sequence: str | ProFormaAnnotation | HasSequence,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (None,),
    calculate_with_composition: bool = False,
    max_ndeltas: int = 1,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[Fragment]: ...


@overload
def fragment(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (None,),
    calculate_with_composition: bool = False,
    max_ndeltas: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[Fragment]]: ...


def fragment(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (),
    calculate_with_composition: bool = False,
    max_ndeltas: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
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
            calculate_with_composition=calculate_with_composition,
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
            calculate_with_composition=calculate_with_composition,
        )


def _frag_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    position: int | tuple[int, int] | None = None,
) -> Fragment:
    annotation = get_annotation_input(sequence=sequence, copy=False)

    return annotation.frag(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        calculate_with_composition=calculate_with_composition,
        position=position,
    )


@overload
def frag(
    sequence: str | ProFormaAnnotation | HasSequence,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    position: int | tuple[int, int] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> Fragment: ...


@overload
def frag(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    position: int | tuple[int, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[Fragment]: ...


def frag(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    position: int | tuple[int, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
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
            calculate_with_composition=calculate_with_composition,
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
            calculate_with_composition=calculate_with_composition,
            position=position,
        )


def _fast_fragment_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: int | Sequence[int] | None = None,
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
    sequence: str | ProFormaAnnotation | HasSequence,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: int | Sequence[int] | None = None,
    *,
    monoisotopic: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> FRAGMENT_MASSES_RETURN: ...


@overload
def fast_fragment(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: int | Sequence[int] | None = None,
    *,
    monoisotopic: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[FRAGMENT_MASSES_RETURN]: ...


def fast_fragment(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: int | Sequence[int] | None = None,
    *,
    monoisotopic: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
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


def _fragment_arrays_chunk(chunk: list[ProFormaAnnotation], **kwargs: Any) -> "dict[str, np.ndarray]":
    return _fragment_arrays(chunk, **kwargs)


def fragment_arrays(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (),
    calculate_with_composition: bool = False,
    max_ndeltas: int = 1,
    min_length: int | None = None,
    max_length: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> "dict[str, np.ndarray]":
    """The ions of :func:`fragment` as numpy columns: a dict of equal-length arrays, one row per ion.

    Needs numpy (``pip install "peptacular[numpy]"``); without it this raises
    :class:`~peptacular.interop.MissingOptionalDependencyError`. Takes the arguments of
    :func:`fragment` and returns the same ions, in the same order, with the same values.
    ``pl.DataFrame(result)``, ``pa.table(result)`` and ``pd.DataFrame(result)`` accept the
    dict as is. The keys are :data:`FRAGMENT_ARRAY_KEYS`:

    - ``peptide_index`` (int64): index of the input peptide (always 0 for a single sequence).
    - ``ion_type`` (str): the ion letter, e.g. ``"b"``, ``"y"``, ``"p"``.
    - ``position`` (int64): the ion number (3 for b3), the start of an internal ion, or 0 for
      an ion with no position (precursor, neutral).
    - ``end_position`` (int64): the end of an internal ion, else 0.
    - ``charge_state`` (int64), ``mz`` (float64; the mass when the charge is 0) and ``mass``
      (float64; the charged mass, as :attr:`Fragment.mass`).
    - ``isotope`` (int64): number of 13C swapped in (the ``isotopes=`` offset), 0 for none.
    - ``isotope_label`` (str): all isotope swaps as text (``"13C^2"``, ``"15N"``), ``""`` for none.
    - ``delta_label`` (str): neutral losses and gains as text, in the form of
      :func:`fragment_records` (water loss is ``"H-2O-1"``, a mass keeps its sign:
      ``"-17.0^2"``), ``""`` for none.
    - ``delta_mass`` (float64): the total mass those deltas add (negative for a loss).

    String columns are numpy ``object`` arrays of ``str``, which numpy 1.26+, polars,
    pyarrow and pandas all accept. Plain a/b/c/x/y/z series (no isotope swap, formula delta
    or neutral loss) are computed with numpy prefix sums in the same arithmetic order as
    :func:`fragment`; other ions are built by :func:`fragment` and copied in.

    A list runs in this process unless ``n_workers`` or ``method`` is given; then it is split
    into contiguous chunks (``chunksize`` peptides each, default one chunk per worker) that
    run in parallel and are joined in input order.

    :raises MissingOptionalDependencyError: If numpy is not installed.

    >>> import peptacular as pt
    >>> cols = pt.fragment_arrays(["PEPTIDE", "PEM[Oxidation]K"], ion_types="b", charges=1)
    >>> cols["peptide_index"].tolist(), cols["position"].tolist()
    ([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1], [1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4])
    >>> cols["mz"][:2].round(4).tolist()
    [98.06, 227.1026]
    """
    kwargs: dict[str, Any] = {
        "ion_types": ion_types,
        "charges": charges,
        "monoisotopic": monoisotopic,
        "isotopes": isotopes,
        "deltas": deltas,
        "neutral_deltas": neutral_deltas,
        "calculate_with_composition": calculate_with_composition,
        "max_ndeltas": max_ndeltas,
        "min_length": min_length,
        "max_length": max_length,
    }
    if not isinstance(sequence, Sequence) or isinstance(sequence, str):
        return _fragment_arrays([get_annotation_input(sequence, copy=False)], **kwargs)
    annotations = [get_annotation_input(item, copy=False) for item in sequence]
    if (method is None and n_workers is None) or len(annotations) < 2:
        return _fragment_arrays(annotations, **kwargs)

    workers = n_workers or getattr(os, "process_cpu_count", os.cpu_count)() or 1
    size = chunksize or max(1, -(-len(annotations) // workers))
    chunks = [annotations[start : start + size] for start in range(0, len(annotations), size)]
    parts = parallel_apply_internal(_fragment_arrays_chunk, chunks, n_workers=n_workers, chunksize=1, method=method, **kwargs)
    import numpy

    offset = 0
    for chunk, part in zip(chunks, parts, strict=True):
        part["peptide_index"] += offset
        offset += len(chunk)
    return {key: numpy.concatenate([part[key] for part in parts]) for key in FRAGMENT_ARRAY_KEYS}

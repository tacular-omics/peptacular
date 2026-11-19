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
from .util import HasSequence, get_annotation_input

__all__ = [
    "mass",
    "mz",
    "comp",
]


def _mass_single(
    sequence: str | ProFormaAnnotation | HasSequence,
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
    sequence: str | ProFormaAnnotation | HasSequence,
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
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
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
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
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

    A positive charge adds protons and a negative charge removes them, so ``PEPTIDE/-2`` is
    the deprotonated ion :math:`[M-2H]^{2-}`.

    ProForma 2.1 words this differently. Section 11.5 reads: "Any positive charges are
    defined as addition of protons H+, any negative charges are defined as addition of
    electrons e-." peptacular follows the specification for positive charges and departs
    from it for negative ones, because the ion negative-mode electrospray produces is the
    deprotonated one. The two readings differ by a neutral hydrogen atom per charge,
    ``PROTON_MASS + ELECTRON_MASS`` or about 1.007825 Da, which is 2.01565 Da at charge 2:

        >>> import peptacular as pt
        >>> round(pt.mass("PEPTIDE/-2"), 6)                              # [M-2H]2-, what this returns
        797.345411
        >>> round(pt.parse("PEPTIDE").neutral_mass(), 6)                 # the neutral peptide
        799.359964
        >>> round(pt.parse("PEPTIDE").neutral_mass() + 2 * pt.ELECTRON_MASS, 6)  # [M+2e]2-, the specification's reading
        799.361061
        >>> round(pt.PROTON_MASS + pt.ELECTRON_MASS, 6)                  # the gap, per charge
        1.007825

    There is no charge-carrier notation for the electron-addition reading, because a charge
    carrier is a charged formula of elements and an electron is not an element: ``/[e:z-1]``
    and ``/[e-1:z-1]`` both raise ``ProFormaFormatError``. A carrier can only spell out the
    deprotonation this already does, and gives the same number:

        >>> round(pt.mass("PEPTIDE/[H-1:z-1^2]"), 6)
        797.345411

    To obtain the specification's reading, add the electrons to the neutral mass yourself,
    as the third example above does.
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
    sequence: str | ProFormaAnnotation | HasSequence,
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
    sequence: str | ProFormaAnnotation | HasSequence,
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
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
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
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
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

    The denominator is the absolute charge, so a negative-mode ion still gets a positive
    *m/z*. The numerator is the deprotonated mass, as in :func:`mass`: ProForma 2.1 section
    11.5 words a negative charge as "addition of electrons e-", and peptacular models
    deprotonation instead, because that is the ion negative-mode electrospray produces. The
    difference is about 1.007825 Da per charge, one neutral hydrogen atom:

        >>> import peptacular as pt
        >>> round(pt.mz("PEPTIDE/-2"), 6)
        398.672706
        >>> round((pt.parse("PEPTIDE").neutral_mass() + 2 * pt.ELECTRON_MASS) / 2, 6)
        399.680531
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
    sequence: str | ProFormaAnnotation | HasSequence,
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
    sequence: str | ProFormaAnnotation | HasSequence,
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
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
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
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
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

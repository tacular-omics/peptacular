from collections.abc import Sequence
from typing import overload

from ..annotation import ProFormaAnnotation
from ..constants import ParallelMethod, ParallelMethodLiteral
from ..property.data import (
    HPLCScale,
    HydrophobicityScale,
    PhysicalPropertyScale,
    PolarityScale,
    SecondaryStructureMethod,
    SecondaryStructureType,
    SurfaceAccessibilityScale,
)
from ..property.types import (
    AggregationMethod,
    AggregationMethodLiteral,
    MissingAAHandling,
    MissingAAHandlingLiteral,
    WeightingMethods,
    WeightingMethodsLiteral,
)
from .parallel import parallel_apply_internal
from .util import HasSequence, get_annotation_input

__all__ = [
    "calc_property",
    "hydrophobicity",
    "flexibility",
    "hydrophilicity",
    "surface_accessibility",
    "polarity",
    "mutability",
    "codons",
    "bulkiness",
    "recognition_factors",
    "transmembrane_tendency",
    "average_buried_area",
    "hplc",
    "refractivity",
    "calc_window_property",
    "charge_at_ph",
    "pi",
    "aa_property_percentage",
    "DEFAULT_AROMATIC_RESIDUES",
    "aromaticity",
    "secondary_structure",
    "alpha_helix_percent",
    "beta_sheet_percent",
    "beta_turn_percent",
    "coil_percent",
    "property_partitions",
]


def _calc_property_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    scale: str | dict[str, float],
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.ERROR,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> float:
    """Calculate property for a single sequence"""
    return get_annotation_input(sequence=sequence, copy=True).prop.calc_property(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )


@overload
def calc_property(
    sequence: str | ProFormaAnnotation | HasSequence,
    scale: str | dict[str, float],
    *,
    missing_aa_handling: MissingAAHandlingLiteral | MissingAAHandling = MissingAAHandling.ERROR,
    aggregation_method: AggregationMethodLiteral | AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WeightingMethodsLiteral | WeightingMethods = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def calc_property(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    scale: str | dict[str, float],
    *,
    missing_aa_handling: MissingAAHandlingLiteral | MissingAAHandling = MissingAAHandling.ERROR,
    aggregation_method: AggregationMethodLiteral | AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WeightingMethodsLiteral | WeightingMethods = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def calc_property(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    scale: str | dict[str, float],
    *,
    missing_aa_handling: MissingAAHandlingLiteral | MissingAAHandling = MissingAAHandling.ERROR,
    aggregation_method: AggregationMethodLiteral | AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WeightingMethodsLiteral | WeightingMethods = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """
    Calculate a physicochemical property for a sequence or list of sequences.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _calc_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )
    else:
        return _calc_property_single(
            sequence=sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )


# Helper function for simple property calculations
def _simple_property_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    scale: str | dict[str, float],
) -> float:
    return get_annotation_input(sequence=sequence, copy=True).prop.calc_property(
        scale=scale,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )


@overload
def hydrophobicity(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def hydrophobicity(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def hydrophobicity(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average hydrophobicity (Kyte-Doolittle).

    The mean of the ``HydrophobicityScale.KYTE_DOOLITTLE`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=HydrophobicityScale.KYTE_DOOLITTLE,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=HydrophobicityScale.KYTE_DOOLITTLE,
        )


@overload
def flexibility(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def flexibility(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def flexibility(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average backbone flexibility (Vihinen).

    The mean of the ``PhysicalPropertyScale.FLEXIBILITY_VIHINEN`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.FLEXIBILITY_VIHINEN,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.FLEXIBILITY_VIHINEN,
        )


@overload
def hydrophilicity(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def hydrophilicity(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def hydrophilicity(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average hydrophilicity (Hopp-Woods).

    The mean of the ``PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD,
        )


@overload
def surface_accessibility(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def surface_accessibility(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def surface_accessibility(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average surface accessibility (Vergoten).

    The mean of the ``SurfaceAccessibilityScale.VERGOTEN`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=SurfaceAccessibilityScale.VERGOTEN,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=SurfaceAccessibilityScale.VERGOTEN,
        )


@overload
def polarity(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def polarity(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def polarity(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average polarity (Grantham).

    The mean of the ``PolarityScale.GRANTHAM`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PolarityScale.GRANTHAM,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PolarityScale.GRANTHAM,
        )


@overload
def mutability(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def mutability(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def mutability(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average relative mutability.

    The mean of the ``PhysicalPropertyScale.MUTABILITY`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.MUTABILITY,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.MUTABILITY,
        )


@overload
def codons(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def codons(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def codons(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average number of codons per residue.

    The mean of the ``PhysicalPropertyScale.CODONS`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.CODONS,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.CODONS,
        )


@overload
def bulkiness(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def bulkiness(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def bulkiness(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average side-chain bulkiness.

    The mean of the ``PhysicalPropertyScale.BULKINESS`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.BULKINESS,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.BULKINESS,
        )


@overload
def recognition_factors(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def recognition_factors(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def recognition_factors(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average recognition factor.

    The mean of the ``PhysicalPropertyScale.RECOGNITION_FACTORS`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.RECOGNITION_FACTORS,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.RECOGNITION_FACTORS,
        )


@overload
def transmembrane_tendency(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def transmembrane_tendency(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def transmembrane_tendency(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average transmembrane tendency.

    The mean of the ``PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY,
        )


@overload
def average_buried_area(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def average_buried_area(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def average_buried_area(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average buried surface area.

    The mean of the ``SurfaceAccessibilityScale.AVERAGE_BURIED_AREA`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=SurfaceAccessibilityScale.AVERAGE_BURIED_AREA,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=SurfaceAccessibilityScale.AVERAGE_BURIED_AREA,
        )


@overload
def hplc(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def hplc(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def hplc(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average HPLC retention coefficient (Meek, pH 2.1).

    The mean of the ``HPLCScale.MEEK_2_1`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=HPLCScale.MEEK_2_1,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=HPLCScale.MEEK_2_1,
        )


@overload
def refractivity(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def refractivity(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def refractivity(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Average refractivity.

    The mean of the ``PhysicalPropertyScale.REFRACTIVITY`` scale over the residues, normalised to 0-1 across the scale. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    :raises PeptacularError: A residue has no value in the scale.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.REFRACTIVITY,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.REFRACTIVITY,
        )


def calc_window_property(
    sequence: str | ProFormaAnnotation | HasSequence,
    scale: str | dict[str, float],
    *,
    window_size: int = 9,
    missing_aa_handling: MissingAAHandlingLiteral | MissingAAHandling = MissingAAHandling.ERROR,
    aggregation_method: AggregationMethodLiteral | AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WeightingMethodsLiteral | WeightingMethods = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> list[float]:
    """Property values over a sliding window along one sequence.

    :param sequence: A ProForma string or annotation.
    :param scale: A scale name from ``PROPERTY_SCALES`` (e.g. ``HydrophobicityScale.KYTE_DOOLITTLE``) or a residue-to-value dict.
    :param window_size: Residues per window.
    :param missing_aa_handling: What to do with residues missing from the scale.
    :param aggregation_method: How the values in a window are combined.
    :param normalize: Scale values to 0-1 across the scale first.
    :param weighting_scheme: Per-position weights inside a window.
    :param min_weight: Smallest weight for non-uniform schemes.
    :param max_weight: Largest weight for non-uniform schemes.
    :return: One value per window, ``len(sequence) - window_size + 1`` values.
    """
    return get_annotation_input(sequence=sequence, copy=True).prop.property_windows(
        scale=scale,
        window_size=window_size,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )


def _charge_at_ph_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    pH: float = 7.0,
) -> float:
    return get_annotation_input(sequence=sequence, copy=False).prop.charge_at_ph(
        pH=pH,
    )


@overload
def charge_at_ph(
    sequence: str | ProFormaAnnotation | HasSequence,
    pH: float = 7.0,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def charge_at_ph(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    pH: float = 7.0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def charge_at_ph(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    pH: float = 7.0,
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Net charge at a given pH, from Henderson-Hasselbalch pKa values.

    Sums the N-terminal, C-terminal and ionisable side-chain contributions. Modifications are ignored.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param pH: The pH.
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: The net charge, or a list of charges for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _charge_at_ph_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            pH=pH,
        )
    else:
        return _charge_at_ph_single(
            sequence=sequence,
            pH=pH,
        )


def _pi_single(
    sequence: str | ProFormaAnnotation | HasSequence,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.prop.pi


@overload
def pi(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def pi(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def pi(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Isoelectric point: the pH at which :func:`charge_at_ph` is zero, found by bisection.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: The pI, or a list of pI values for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _pi_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _pi_single(
            sequence=sequence,
        )


def _aa_property_percentage_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    residues: list[str],
) -> float:
    return get_annotation_input(sequence=sequence, copy=False).prop.aa_property_percentage(
        residues=residues,
    )


@overload
def aa_property_percentage(
    sequence: str | ProFormaAnnotation | HasSequence,
    residues: list[str],
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def aa_property_percentage(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    residues: list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def aa_property_percentage(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    residues: list[str],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Fraction (0-1) of residues that are in ``residues``.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param residues: One-letter amino-acid codes to count.
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _aa_property_percentage_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            residues=residues,
        )
    else:
        return _aa_property_percentage_single(
            sequence=sequence,
            residues=residues,
        )


DEFAULT_AROMATIC_RESIDUES = ["Y", "W", "F"]


@overload
def aromaticity(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    aromatic_residues: list[str] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def aromaticity(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    aromatic_residues: list[str] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def aromaticity(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    aromatic_residues: list[str] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Fraction (0-1) of aromatic residues (Y, W, F by default).

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param aromatic_residues: Residues counted as aromatic; defaults to ``['Y', 'W', 'F']``.
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    """
    if aromatic_residues is None:
        aromatic_residues = DEFAULT_AROMATIC_RESIDUES
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _aa_property_percentage_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            residues=aromatic_residues,
        )
    else:
        return _aa_property_percentage_single(
            sequence=sequence,
            residues=aromatic_residues,
        )


def _secondary_structure_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
) -> dict[str, float]:
    return get_annotation_input(sequence=sequence, copy=True).prop.secondary_structure(
        scale=scale,
    )


@overload
def secondary_structure(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> dict[str, float]: ...


@overload
def secondary_structure(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[dict[str, float]]: ...


def secondary_structure(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> dict[str, float] | list[dict[str, float]]:
    """Predicted secondary-structure fractions.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param scale: A ``SecondaryStructureMethod``; defaults to Deleage-Roux.
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A dict of ``SecondaryStructureType`` to fraction (alpha helix, beta sheet, beta turn, coil), or a list of dicts for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _secondary_structure_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=scale,
        )
    else:
        return _secondary_structure_single(
            sequence=sequence,
            scale=scale,
        )


def _alpha_helix_percent_single(
    sequence: str | ProFormaAnnotation | HasSequence,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.ALPHA_HELIX]


@overload
def alpha_helix_percent(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def alpha_helix_percent(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def alpha_helix_percent(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Predicted alpha helix fraction (0-1) by the Deleage-Roux method.

    The ``alpha_helix`` entry of :func:`secondary_structure` with the default scale.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _alpha_helix_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _alpha_helix_percent_single(
            sequence=sequence,
        )


def _beta_sheet_percent_single(
    sequence: str | ProFormaAnnotation | HasSequence,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.BETA_SHEET]


@overload
def beta_sheet_percent(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def beta_sheet_percent(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def beta_sheet_percent(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Predicted beta sheet fraction (0-1) by the Deleage-Roux method.

    The ``beta_sheet`` entry of :func:`secondary_structure` with the default scale.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _beta_sheet_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _beta_sheet_percent_single(
            sequence=sequence,
        )


def _beta_turn_percent_single(
    sequence: str | ProFormaAnnotation | HasSequence,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.BETA_TURN]


@overload
def beta_turn_percent(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def beta_turn_percent(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def beta_turn_percent(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Predicted beta turn fraction (0-1) by the Deleage-Roux method.

    The ``beta_turn`` entry of :func:`secondary_structure` with the default scale.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _beta_turn_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _beta_turn_percent_single(
            sequence=sequence,
        )


def _coil_percent_single(
    sequence: str | ProFormaAnnotation | HasSequence,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.COIL]


@overload
def coil_percent(
    sequence: str | ProFormaAnnotation | HasSequence,
    *,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float: ...


@overload
def coil_percent(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


def coil_percent(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    *,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> float | list[float]:
    """Predicted coil fraction (0-1) by the Deleage-Roux method.

    The ``coil`` entry of :func:`secondary_structure` with the default scale.

    :param sequence: A ProForma string or annotation, or a list of them (lists run in parallel above ``AUTO_PARALLEL_MIN_ITEMS``).
    :param n_workers: Worker count for list input.
    :param chunksize: Items per worker task for list input.
    :param method: Parallel backend for list input (``process``, ``thread``, ``sequential``); ``None`` chooses automatically.
    :return: A float, or a list of floats for list input.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _coil_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _coil_percent_single(
            sequence=sequence,
        )


def _property_partitions_single(
    sequence: str | ProFormaAnnotation | HasSequence,
    scale: str | dict[str, float],
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.AVG,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> list[float]:
    return get_annotation_input(sequence=sequence, copy=True).prop.property_partitions(
        scale=scale,
        num_windows=num_windows,
        aa_overlap=aa_overlap,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )


@overload
def property_partitions(
    sequence: str | ProFormaAnnotation | HasSequence,
    scale: str | dict[str, float],
    *,
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: MissingAAHandlingLiteral | MissingAAHandling = MissingAAHandling.AVG,
    aggregation_method: AggregationMethodLiteral | AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WeightingMethodsLiteral | WeightingMethods = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: None = None,
    chunksize: None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float]: ...


@overload
def property_partitions(
    sequence: Sequence[str | ProFormaAnnotation | HasSequence],
    scale: str | dict[str, float],
    *,
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: MissingAAHandlingLiteral | MissingAAHandling = MissingAAHandling.AVG,
    aggregation_method: AggregationMethodLiteral | AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WeightingMethodsLiteral | WeightingMethods = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[list[float]]: ...


def property_partitions(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    scale: str | dict[str, float],
    *,
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: MissingAAHandlingLiteral | MissingAAHandling = MissingAAHandling.AVG,
    aggregation_method: AggregationMethodLiteral | AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WeightingMethodsLiteral | WeightingMethods = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParallelMethod | ParallelMethodLiteral | None = None,
) -> list[float] | list[list[float]]:
    """Generate property values for N number of sliding windows across the sequence.

    Divides the sequence into N overlapping windows and calculates property values
    for each window. Useful for analyzing local variations in peptide properties.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _property_partitions_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=scale,
            num_windows=num_windows,
            aa_overlap=aa_overlap,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )
    else:
        return _property_partitions_single(
            sequence=sequence,
            scale=scale,
            num_windows=num_windows,
            aa_overlap=aa_overlap,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )

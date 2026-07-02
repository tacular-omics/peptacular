from collections.abc import Sequence
from typing import overload

from ..annotation import ProFormaAnnotation
from ..constants import parallelMethod, parallelMethodLiteral
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
from .util import get_annotation_input


def _calc_property_single(
    sequence: str | ProFormaAnnotation,
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
    sequence: str | ProFormaAnnotation,
    scale: str | dict[str, float],
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.ERROR,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def calc_property(
    sequence: Sequence[str | ProFormaAnnotation],
    scale: str | dict[str, float],
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.ERROR,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def calc_property(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    scale: str | dict[str, float],
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.ERROR,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
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
    sequence: str | ProFormaAnnotation,
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def hydrophobicity(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def hydrophobicity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def flexibility(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def flexibility(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def hydrophilicity(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def hydrophilicity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def surface_accessibility(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def surface_accessibility(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def polarity(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def polarity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def mutability(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def mutability(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def codons(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def codons(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def bulkiness(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def bulkiness(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def recognition_factors(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def recognition_factors(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def transmembrane_tendency(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def transmembrane_tendency(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def average_buried_area(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def average_buried_area(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def hplc(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def hplc(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def refractivity(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def refractivity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    scale: str | dict[str, float],
    window_size: int = 9,
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.ERROR,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> list[float]:
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
    sequence: str | ProFormaAnnotation,
    pH: float = 7.0,
) -> float:
    return get_annotation_input(sequence=sequence, copy=False).prop.charge_at_ph(
        pH=pH,
    )


@overload
def charge_at_ph(
    sequence: str | ProFormaAnnotation,
    pH: float = 7.0,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def charge_at_ph(
    sequence: Sequence[str | ProFormaAnnotation],
    pH: float = 7.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def charge_at_ph(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    pH: float = 7.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.prop.pi


@overload
def pi(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def pi(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def pi(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    residues: list[str],
) -> float:
    return get_annotation_input(sequence=sequence, copy=False).prop.aa_property_percentage(
        residues=residues,
    )


@overload
def aa_property_percentage(
    sequence: str | ProFormaAnnotation,
    residues: list[str],
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def aa_property_percentage(
    sequence: Sequence[str | ProFormaAnnotation],
    residues: list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def aa_property_percentage(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    residues: list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    aromatic_residues: list[str] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def aromaticity(
    sequence: Sequence[str | ProFormaAnnotation],
    aromatic_residues: list[str] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def aromaticity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    aromatic_residues: list[str] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
) -> dict[str, float]:
    return get_annotation_input(sequence=sequence, copy=True).prop.secondary_structure(
        scale=scale,
    )


@overload
def secondary_structure(
    sequence: str | ProFormaAnnotation,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> dict[str, float]: ...


@overload
def secondary_structure(
    sequence: Sequence[str | ProFormaAnnotation],
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[dict[str, float]]: ...


def secondary_structure(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> dict[str, float] | list[dict[str, float]]:
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
    sequence: str | ProFormaAnnotation,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.ALPHA_HELIX]


@overload
def alpha_helix_percent(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def alpha_helix_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def alpha_helix_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.BETA_SHEET]


@overload
def beta_sheet_percent(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def beta_sheet_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def beta_sheet_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.BETA_TURN]


@overload
def beta_turn_percent(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def beta_turn_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def beta_turn_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
) -> float:
    d = _secondary_structure_single(sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX)
    return d[SecondaryStructureType.COIL]


@overload
def coil_percent(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float: ...


@overload
def coil_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


def coil_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> float | list[float]:
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
    sequence: str | ProFormaAnnotation,
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
    sequence: str | ProFormaAnnotation,
    scale: str | dict[str, float],
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.AVG,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[float]: ...


@overload
def property_partitions(
    sequence: Sequence[str | ProFormaAnnotation],
    scale: str | dict[str, float],
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.AVG,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[list[float]]: ...


def property_partitions(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    scale: str | dict[str, float],
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.AVG,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
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

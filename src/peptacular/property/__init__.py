"""Protein property calculation module."""

from typing import Any

from .data import (
    COMPOSITION_SCALES,
    FLEXIBILITY_SCALES,
    HPLC_SCALES,
    HYDROPHILICITY_SCALES,
    HYDROPHOBICITY_SCALES,
    PHYSICAL_PROPERTY_SCALES,
    POLARITY_SCALES,
    PROPERTY_SCALES,
    SURFACE_ACCESSIBILITY_SCALES,
    BetaStrandScale,
    ChargeScale,
    CompositionScale,
    HPLCScale,
    HydrophobicityScale,
    PhysicalPropertyScale,
    PolarityScale,
    PropertyScale,
    SecondaryStructureMethod,
    SecondaryStructureScale,
    SecondaryStructureType,
    SurfaceAccessibilityScale,
)
from .prop import AnnotationProperties
from .types import (
    AggregationMethod,
    AggregationMethodLiteral,
    MissingAAHandling,
    MissingAAHandlingLiteral,
    WeightingMethods,
    WeightingMethodsLiteral,
)

__all__ = [
    # main class
    "AnnotationProperties",
    # data
    "SecondaryStructureMethod",
    "SecondaryStructureType",
    "PropertyScale",
    "HydrophobicityScale",
    "SecondaryStructureScale",
    "SurfaceAccessibilityScale",
    "ChargeScale",
    "PolarityScale",
    "HPLCScale",
    "BetaStrandScale",
    "PhysicalPropertyScale",
    "CompositionScale",
    "PROPERTY_SCALES",
    "HYDROPHOBICITY_SCALES",
    "SURFACE_ACCESSIBILITY_SCALES",
    "HPLC_SCALES",
    "HYDROPHILICITY_SCALES",
    "FLEXIBILITY_SCALES",
    "POLARITY_SCALES",
    "COMPOSITION_SCALES",
    "PHYSICAL_PROPERTY_SCALES",
    # types
    "MissingAAHandling",
    "MissingAAHandlingLiteral",
    "AggregationMethod",
    "AggregationMethodLiteral",
    "WeightingMethods",
    "WeightingMethodsLiteral",
]


def __getattr__(name: str) -> Any:
    """Forward deprecated aliases (such as ``FLIXIBILITY_SCALES``) to :mod:`peptacular.property.data`."""
    from . import data

    if name in data._DEPRECATED_ALIASES:
        return data._resolve_deprecated_alias(name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

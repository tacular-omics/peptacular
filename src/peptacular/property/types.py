from __future__ import annotations

from enum import StrEnum
from typing import Literal

from ..diagnostics import PeptacularError

__all__ = [
    "MissingAAHandling",
    "MissingAAHandlingLiteral",
    "AggregationMethod",
    "AggregationMethodLiteral",
    "WeightingMethods",
    "WeightingMethodsLiteral",
]


class MissingAAHandling(StrEnum):
    """Strategy for handling missing amino acid values"""

    ZERO = "zero"
    AVG = "avg"
    MIN = "min"
    MAX = "max"
    MEDIAN = "median"
    ERROR = "error"
    SKIP = "skip"

    @staticmethod
    def from_str(label: str) -> MissingAAHandling:
        """Convert a string to a MissingAAHandling enum member."""
        label = label.lower()
        for method in MissingAAHandling:
            if method.value == label:
                return method
        raise PeptacularError(f"Unknown MissingAAHandling: {label}")


MissingAAHandlingLiteral = Literal["zero", "avg", "min", "max", "median", "error", "skip"]


class AggregationMethod(StrEnum):
    """Strategy for aggregating amino acid values"""

    SUM = "sum"
    AVG = "avg"

    @staticmethod
    def from_str(label: str) -> AggregationMethod:
        """Convert a string to an AggregationMethod enum member."""
        label = label.lower()
        for method in AggregationMethod:
            if method.value == label:
                return method
        raise PeptacularError(f"Unknown AggregationMethod: {label}")


AggregationMethodLiteral = Literal["sum", "avg"]


class WeightingMethods(StrEnum):
    """Enumeration of possible weighting methods for amino acid properties"""

    UNIFORM = "uniform"
    LINEAR = "linear"
    EXPONENTIAL = "exponential"
    GAUSSIAN = "gaussian"
    SIGMOID = "sigmoid"
    COSINE = "cosine"
    SINUSOIDAL = "sinusoidal"

    @staticmethod
    def from_str(label: str) -> WeightingMethods:
        """Convert a string to a WeightingMethods enum member."""
        label = label.lower()
        for method in WeightingMethods:
            if method.value == label:
                return method
        raise PeptacularError(f"Unknown WeightingMethods: {label}")


WeightingMethodsLiteral = Literal["uniform", "linear", "exponential", "gaussian", "sigmoid", "cosine", "sinusoidal"]

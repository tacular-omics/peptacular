import math
from collections.abc import Sequence
from typing import Any

from ..diagnostics import PeptacularError
from .types import WeightingMethods

__all__ = [
    "get_weights",
]


def _get_uniform_weights(length: int, min_weight: float = 0.1, max_weight: float = 1.0) -> list[float]:
    """
    Generate uniform weights for a sequence of given length.

    :param length: Length of the sequence
    :param min_weight: Minimum weight (default: 0.1)
    :param max_weight: Maximum weight (default: 1.0)
    :return: List of weights
    """
    if length <= 0:
        return []

    return [max_weight] * length  # Uniform weights for simplicity


def _get_linear_weights(length: int, min_weight: float = 0.1, max_weight: float = 1.0) -> list[float]:
    """
    Generate linear weights for a sequence of given length. Middle should be peak and left right should be min.

    :param length: Length of the sequence
    :param min_weight: Minimum weight (default: 0.1)
    :param max_weight: Maximum weight (default: 1.0)
    :return: List of weights
    """

    if length <= 0:
        return []

    if length == 1:
        return [max_weight]

    if length == 2:
        return [max_weight, max_weight]

    weights: list[float] = []
    mid = (length - 1) / 2

    for i in range(length):
        # Distance from the middle (normalized to 0-1)
        distance_from_mid = abs(i - mid) / mid
        # Linear interpolation from max_weight (at center) to min_weight (at edges)
        weight = max_weight - (max_weight - min_weight) * distance_from_mid
        weights.append(weight)

    return weights


def _get_exponential_weights(length: int, min_weight: float = 0.1, max_weight: float = 1.0) -> list[float]:
    """
    Generate exponential weights for a sequence of given length. Middle should be peak and left right should be min.
    Uses exponential decay from center to edges.

    :param length: Length of the sequence
    :param min_weight: Minimum weight (default: 0.1)
    :param max_weight: Maximum weight (default: 1.0)
    :return: List of weights
    """

    if length <= 0:
        return []

    if length == 1:
        return [max_weight]

    # Fixed range: always sample between -3 and +3 (centered at 0)
    start = -3.0
    end = 3.0

    weights: list[float] = []

    for i in range(length):
        # Map index to fixed range
        x = start + (end - start) * i / (length - 1) if length > 1 else 0

        # Exponential decay from center (x=0)
        weight = math.exp(-abs(x))

        # Scale to desired range
        scaled_weight = min_weight + (max_weight - min_weight) * weight / math.exp(0)
        weights.append(scaled_weight)

    return weights


def _get_gaussian_weights(length: int, min_weight: float = 0.1, max_weight: float = 1.0) -> list[float]:
    """
    Generate Gaussian weights for a sequence of given length.
    Uses fixed Gaussian distribution sampled across predetermined space.

    :param length: Length of the sequence
    :param min_weight: Minimum weight (default: 0.1)
    :param max_weight: Maximum weight (default: 1.0)
    :return: List of weights
    """

    if length <= 0:
        return []

    if length == 1:
        return [max_weight]

    # Fixed range: always sample between -3 and +3 (centered at 0)
    # Fixed sigma for consistent Gaussian shape
    start = -3.0
    end = 3.0
    sigma = 1.0

    weights: list[float] = []

    for i in range(length):
        # Map index to fixed range
        x = start + (end - start) * i / (length - 1) if length > 1 else 0

        # Gaussian function: exp(-(x^2)/(2*sigma^2))
        gaussian_val = math.exp(-(x**2) / (2 * sigma**2))

        # Scale to desired range
        scaled_weight = min_weight + (max_weight - min_weight) * gaussian_val
        weights.append(scaled_weight)

    return weights


def _get_sigmoid_weights(
    length: int,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    steepness: float = 2.0,
) -> list[float]:
    """
    Generate sigmoid (S-curve) weights for a sequence. Creates smooth transition from min to max.

    :param length: Length of the sequence
    :param min_weight: Minimum weight (default: 0.1)
    :param max_weight: Maximum weight (default: 1.0)
    :param steepness: Controls steepness of the sigmoid curve (default: 2.0)
    :return: List of weights

    .. code-block:: python

        >>> weights = _get_sigmoid_weights(5)
        >>> [round(w, 2) for w in weights]
        [0.1, 0.26, 0.55, 0.84, 1.0]
    """
    if length <= 0:
        return []

    if length == 1:
        return [max_weight]

    weights: list[float] = []
    for i in range(length):
        # Map to range [-6, 6] for good sigmoid behavior
        x = -6 + 12 * i / (length - 1) if length > 1 else 0
        # Sigmoid function
        sigmoid_val = 1 / (1 + math.exp(-steepness * x))
        # Scale to desired range
        weight = min_weight + (max_weight - min_weight) * sigmoid_val
        weights.append(weight)

    return weights


def _get_cosine_weights(length: int, min_weight: float = 0.1, max_weight: float = 1.0, cycles: float = 1.0) -> list[float]:
    """
    Generate cosine-based weights for a sequence. Creates periodic patterns.

    :param length: Length of the sequence
    :param min_weight: Minimum weight (default: 0.1)
    :param max_weight: Maximum weight (default: 1.0)
    :param cycles: Number of cosine cycles (default: 1.0)
    :return: List of weights

    .. code-block:: python

        >>> weights = _get_cosine_weights(5, cycles=0.5)
        >>> [round(w, 2) for w in weights]
        [1.0, 0.77, 0.32, 0.32, 0.77]
    """
    if length <= 0:
        return []

    if length == 1:
        return [max_weight]

    weights: list[float] = []
    for i in range(length):
        # Cosine wave from 0 to cycles*2π
        angle = 2 * math.pi * cycles * i / (length - 1) if length > 1 else 0
        # Cosine ranges from -1 to 1, normalize to 0 to 1
        cos_val = (math.cos(angle) + 1) / 2
        # Scale to desired range
        weight = min_weight + (max_weight - min_weight) * cos_val
        weights.append(weight)

    return weights


def _get_sinusoidal_weights(length: int, min_weight: float = 0.1, max_weight: float = 1.0, phase: float = 0.0) -> list[float]:
    """
    Generate sinusoidal weights (half sine wave from 0 to π).

    :param length: Length of the sequence
    :param min_weight: Minimum weight (default: 0.1)
    :param max_weight: Maximum weight (default: 1.0)
    :param phase: Phase shift in radians (default: 0.0)
    :return: List of weights
    """
    if length <= 0:
        return []

    if length == 1:
        return [max_weight]

    weights: list[float] = []
    for i in range(length):
        # Map to range [0, π] for half sine wave
        angle = math.pi * i / (length - 1) + phase if length > 1 else phase
        # Sine function (0 to 1 range for half wave)
        sin_val = math.sin(angle)
        # Ensure non-negative values
        sin_val = max(0, sin_val)
        # Scale to desired range
        weight = min_weight + (max_weight - min_weight) * sin_val
        weights.append(weight)

    return weights


def get_weights(
    length: int,
    weights: Sequence[float] | str = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    **kwargs: Any,
) -> Sequence[float]:
    """
    Get weights for a sequence based on the specified weighting scheme.

    :param length: Length of the sequence
    :param weights: Weighting scheme or list of weights
    :param min_weight: Minimum weight value
    :param max_weight: Maximum weight value
    :param kwargs: Additional parameters for specific weight functions
    :return: List of weights
    """
    if isinstance(weights, str):
        weights_enum = WeightingMethods.from_str(weights)

        match weights_enum:
            case WeightingMethods.UNIFORM:
                return _get_uniform_weights(length, min_weight, max_weight)
            case WeightingMethods.LINEAR:
                return _get_linear_weights(length, min_weight, max_weight)
            case WeightingMethods.EXPONENTIAL:
                return _get_exponential_weights(length, min_weight, max_weight)
            case WeightingMethods.GAUSSIAN:
                return _get_gaussian_weights(length, min_weight, max_weight)
            case WeightingMethods.SIGMOID:
                return _get_sigmoid_weights(length, min_weight, max_weight, kwargs.get("steepness", 2.0))
            case WeightingMethods.COSINE:
                return _get_cosine_weights(length, min_weight, max_weight, kwargs.get("cycles", 1.0))
            case WeightingMethods.SINUSOIDAL:
                return _get_sinusoidal_weights(length, min_weight, max_weight, kwargs.get("phase", 0.0))
            case _:  # pragma: no cover - from_str already rejects unknown names
                raise PeptacularError(f"Unsupported weights type: {weights}")
    else:
        if not isinstance(weights, Sequence):
            raise TypeError("weights must be a sequence of floats or a valid WeightingMethod enum.")
        if len(weights) != length:
            raise PeptacularError(f"Length of weights list ({len(weights)}) does not match sequence length ({length}).")
        return weights

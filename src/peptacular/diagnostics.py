"""Structured, serializable diagnostics for sequence calculations."""

from dataclasses import dataclass
from typing import Literal

__all__ = ["Diagnostic", "CompositionError", "InvalidAdjustmentError", "UnsupportedOperationError", "UnknownModificationError"]


class UnknownModificationError(ValueError):
    """A modification could not be resolved in the reference vocabularies."""


class CompositionError(ValueError):
    """The requested elemental composition is not available."""


class InvalidAdjustmentError(ValueError):
    """An adjustment has invalid counts or produces an impossible composition."""


class UnsupportedOperationError(ValueError):
    """The requested operation does not support this input."""


@dataclass(frozen=True)
class Diagnostic:
    """An input failure, with a stable code and the original exception message.

    :param code: Machine-readable category, independent of message wording.
    :param stage: Whether parsing, validation, or calculation failed.
    :param message: Human-readable explanation from the original exception.
    :param exception_type: Original exception class name.
    """

    code: str
    stage: Literal["parse", "validate", "calculate"]
    message: str
    exception_type: str


def diagnostic_from_exception(exc: Exception, stage: Literal["parse", "validate", "calculate"]) -> Diagnostic:
    if isinstance(exc, UnknownModificationError):
        code = "unresolved_modification"
    elif isinstance(exc, CompositionError):
        code = "unavailable_composition"
    elif isinstance(exc, InvalidAdjustmentError):
        code = "invalid_adjustment"
    elif isinstance(exc, UnsupportedOperationError):
        code = "unsupported_operation"
    elif stage == "parse":
        code = "invalid_notation"
    elif stage == "validate":
        code = "invalid_annotation"
    elif isinstance(exc, KeyError):
        code = "unresolved_reference"
    else:
        code = "calculation_error"
    return Diagnostic(code, stage, str(exc), type(exc).__name__)

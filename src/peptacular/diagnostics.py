"""Structured, serializable diagnostics for sequence calculations."""

from dataclasses import dataclass
from typing import Literal

__all__ = [
    "Diagnostic",
    "PeptacularError",
    "ProFormaFormatError",
    "CompositionError",
    "InvalidAdjustmentError",
    "InvalidPositionError",
    "FastaFormatError",
    "UnsupportedOperationError",
    "UnknownModificationError",
]


class PeptacularError(ValueError):
    """Base class for peptacular's typed errors.

    It subclasses ``ValueError``, so code that catches ``ValueError`` keeps working.
    Catch ``PeptacularError`` to handle any peptacular input error at once.
    """


class ProFormaFormatError(PeptacularError):
    """The input is not valid ProForma notation.

    Raised by ``parse`` and ``parse_chimeric``, and by calculations that parse a
    modification, glycan, isotope label or adduct lazily (e.g. ``mass("<113C>PEP")``).
    """


class UnknownModificationError(PeptacularError):
    """A modification could not be resolved in the reference vocabularies."""


class CompositionError(PeptacularError):
    """The requested elemental composition or mass is not available (e.g. an empty sequence)."""


class InvalidAdjustmentError(PeptacularError):
    """An adjustment has invalid counts or produces an impossible composition."""


class InvalidPositionError(PeptacularError):
    """A slice index or fragment position is outside the sequence."""


class FastaFormatError(PeptacularError):
    """The input is not valid FASTA text."""


class UnsupportedOperationError(PeptacularError):
    """The requested operation does not support this input (e.g. an unknown ion type)."""


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
    elif isinstance(exc, ProFormaFormatError):
        code = "invalid_notation"
    elif stage == "parse":
        code = "invalid_notation"
    elif stage == "validate":
        code = "invalid_annotation"
    elif isinstance(exc, KeyError):
        code = "unresolved_reference"
    else:
        code = "calculation_error"
    return Diagnostic(code, stage, str(exc), type(exc).__name__)

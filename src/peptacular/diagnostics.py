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
    "UnknownEnzymeError",
    "UnknownElementError",
    "PeptacularKeyError",
]


class PeptacularError(ValueError):
    """Base class for every error peptacular raises for bad input.

    It subclasses ``ValueError``, so code that catches ``ValueError`` keeps working.
    Catch ``PeptacularError`` to handle any peptacular input error at once. A
    ``TypeError`` (an argument of the wrong Python type) is not wrapped.
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


class InvalidPositionError(PeptacularError, IndexError):
    """A slice index, residue index or fragment position is outside the sequence.

    Also an ``IndexError``, so ``except IndexError`` keeps working.
    """


class PeptacularKeyError(PeptacularError, KeyError):
    """A name was not found in a lookup table (a protease, an element symbol, ...).

    Also a ``KeyError``, so ``except KeyError`` keeps working.
    """

    def __str__(self) -> str:
        # KeyError.__str__ would repr() the message; keep it readable.
        return Exception.__str__(self)


class UnknownEnzymeError(PeptacularKeyError):
    """An ``enzyme`` string names no protease in tacular's ``PROTEASE_LOOKUP``.

    Also a ``KeyError``. To digest with a custom cleavage rule, pass a compiled
    pattern (``re.compile(...)``) instead of a string.
    """


class UnknownElementError(PeptacularKeyError):
    """An element or isotope symbol is not in tacular's ``ELEMENT_LOOKUP``. Also a ``KeyError``."""


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
    elif isinstance(exc, UnknownEnzymeError):
        code = "unknown_enzyme"
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


def lookup_element(symbol: str):
    """Return tacular's ``ElementInfo`` for ``symbol``, raising :class:`UnknownElementError` if unknown."""
    from tacular import ELEMENT_LOOKUP

    try:
        return ELEMENT_LOOKUP[symbol]
    except KeyError:
        raise UnknownElementError(f"{symbol!r} is not a known element or isotope symbol") from None

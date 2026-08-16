"""Exceptions and warnings shared by optional interoperability adapters."""


class InteropError(Exception):
    """Base class for interoperability failures."""


class MissingOptionalDependencyError(InteropError, ImportError):
    """Raised when an adapter's third-party package is not installed."""


class InteropConversionError(InteropError, ValueError):
    """Raised when a value cannot be represented by the target package."""


class LossyConversionWarning(UserWarning):
    """Warn that unsupported annotation information was discarded."""

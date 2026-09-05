"""Information-loss policies for interoperability adapters."""

from enum import StrEnum


class LossPolicy(StrEnum):
    """How an adapter should handle information its target cannot represent."""

    ERROR = "error"
    WARN = "warn"
    DROP = "drop"

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Protocol, Self

__all__ = [
    "EnzymeConfig",
    "DigestProtocol",
]


@dataclass(frozen=True, slots=True)
class EnzymeConfig:
    """One enzyme step of a sequential digestion.

    :param enzyme: A protease name from tacular's ``PROTEASE_LOOKUP`` or a compiled pattern.
    :param missed_cleavages: Maximum missed cleavages for this step.
    :param semi_enzymatic: Also produce semi-enzymatic peptides in this step.
    :param complete_digestion: If False, the undigested input is kept as well.
    """

    enzyme: str | re.Pattern[str]
    missed_cleavages: int = 0
    semi_enzymatic: bool = False
    complete_digestion: bool = True


class DigestProtocol(Protocol):
    """Protocol defining the interface for objects that can be digested."""

    @property
    def stripped_sequence(self) -> str:
        """The sequence without modifications (for calculations)."""
        ...

    def slice(
        self,
        start: int | None,
        stop: int | None,
        inplace: bool = False,
    ) -> Self: ...

    def serialize(self) -> str:
        """Serialize the annotation to a ProForma string."""
        ...

    def __len__(self) -> int:
        """Return the length of the sequence."""
        ...

    def has_mods(self) -> bool:
        """Return True if the sequence has modifications."""
        ...

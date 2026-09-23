import sys
from collections.abc import Iterable
from typing import Protocol, runtime_checkable

from .constants import ModType, ModTypeLiteral

__all__ = [
    "SupportsStr",
    "handle_number_and_intern_mod",
    "get_mod_type",
]


@runtime_checkable
class SupportsStr(Protocol):
    """Protocol for any object that can be converted to string"""

    def __str__(self) -> str: ...


def handle_number_and_intern_mod(mod: SupportsStr | float | int) -> str:
    """Validate and intern a modification string"""
    if isinstance(mod, (float, int)):
        mod_str = f"{mod:+}"
    else:
        mod_str = str(mod).strip()
        if not mod_str:
            raise ValueError("Empty modification string is not allowed")
    return sys.intern(mod_str)


def get_mod_type(mod: ModTypeLiteral | ModType | str) -> ModType:
    """
    Convert a modification type name to its :class:`ModType` member.

    >>> from peptacular.utils import get_mod_type
    >>> get_mod_type("internal")
    <ModType.INTERNAL: 'internal'>
    >>> get_mod_type("nterm")
    <ModType.NTERM: 'nterm'>

    :param mod: A ModType, or its string value (e.g. ``"nterm"``, ``"internal"``).
    :type mod: ModTypeLiteral | ModType | str
    :return: The matching ModType member.
    :rtype: ModType
    :raises ValueError: If ``mod`` is not a string or ModType, or names no ModType.
    """
    if isinstance(mod, ModType):
        return mod

    if isinstance(mod, str):
        for mod_type in ModType:
            if mod_type.value == mod:
                return mod_type
    else:
        raise ValueError(f"mod must be a string or ModType, got {type(mod)}")

    raise ValueError(f"Unknown mod type: {mod}")


def get_mods(
    mods: Iterable[ModTypeLiteral] | Iterable[ModType] | ModType | ModTypeLiteral | None,
) -> list[ModType]:
    """
    Get the list of modification types from the input.

    :param mods: Modification types as a ModType, iterable of ModTypes, or None.
    :type mods: None | ModType | Iterable[ModType]
    :return: List of ModType Enum values.
    :rtype: list[ModType]
    :raises ValueError: If mods is not None, ModType, or iterable of ModTypes.
    """

    if mods is None:
        return [mod_type for mod_type in ModType]
    elif isinstance(mods, (str, ModType)):
        # Single modification type
        return [get_mod_type(mods)]
    elif isinstance(mods, Iterable):
        # List of modification types
        return [get_mod_type(mod) for mod in mods]

    raise ValueError(f"mods parameter must be str, list of str, or None, got {type(mods)}")

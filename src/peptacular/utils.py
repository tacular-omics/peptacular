"""Helpers for :class:`~peptacular.constants.ModType` arguments."""

from collections.abc import Iterable

from .constants import ModType, ModTypeLiteral
from .diagnostics import PeptacularError

__all__ = [
    "get_mod_type",
]


def get_mod_type(mod: ModTypeLiteral | ModType | str) -> ModType:
    """Convert a modification type name to its :class:`ModType` member.

    ``ModType`` members pass through unchanged; strings are matched against the member
    values (``"nterm"``, ``"cterm"``, ``"internal"``, ``"interval"``, ``"isotope"``,
    ``"static"``, ``"labile"``, ``"unknown"``, ``"charge"``). Matching is exact and
    case-sensitive.

    >>> from peptacular import get_mod_type
    >>> get_mod_type("internal")
    <ModType.INTERNAL: 'internal'>
    >>> get_mod_type("nterm")
    <ModType.NTERM: 'nterm'>

    :param mod: A ModType, or its string value.
    :type mod: ModTypeLiteral | ModType | str
    :return: The matching ModType member.
    :rtype: ModType
    :raises TypeError: If ``mod`` is not a string or ModType.
    :raises PeptacularError: If ``mod`` is a string that names no ModType.
    """
    if isinstance(mod, ModType):
        return mod
    if not isinstance(mod, str):
        raise TypeError(f"mod must be a string or ModType, got {type(mod).__name__}")
    try:
        return ModType(mod)
    except ValueError:
        valid = ", ".join(repr(m.value) for m in ModType)
        raise PeptacularError(f"Unknown mod type: {mod!r} (expected one of {valid})") from None


def _resolve_mod_types(
    mods: Iterable[ModTypeLiteral] | Iterable[ModType] | ModType | ModTypeLiteral | None,
) -> list[ModType]:
    """Normalise a ``mods``/``mod_types`` argument to a list of ModType members (``None`` means all)."""
    if mods is None:
        return list(ModType)
    if isinstance(mods, (str, ModType)):
        return [get_mod_type(mods)]
    if isinstance(mods, Iterable):
        return [get_mod_type(mod) for mod in mods]
    raise TypeError(f"mods must be a ModType, a string, an iterable of them, or None, got {type(mods).__name__}")

"""Guard against descriptors leaking into dataclass field defaults.

:mod:`dataclasses` resolves a field's default with ``getattr(cls, name, MISSING)``, which
walks the MRO. An attribute inherited from a base class -- a property on a mixin, say --
is therefore picked up as that field's default. Slot descriptors are the one kind
:mod:`dataclasses` filters out; a property, a ``classmethod`` or a ``staticmethod`` is
stored as the default and the field silently becomes optional.

This walks every module in the package so a dataclass added later is covered the day it
lands.

:func:`_is_descriptor` is deliberately broad, so it also matches a callable: a plain
function or a lambda, because ``FunctionType`` carries ``__get__``, and from Python 3.13
a :func:`functools.partial` too (on 3.12 it is not a descriptor). No field default in
the package is a callable. If one ever legitimately is, loosen this test rather than
changing the code.
"""

import dataclasses
import importlib
import pkgutil
import types

import pytest

import peptacular as pt

# dataclasses replaces a slot descriptor default with MISSING, so it can never leak.
SAFE_DESCRIPTOR_TYPES = (types.MemberDescriptorType,)


def _module_names() -> list[str]:
    names = {"peptacular"}
    names.update(info.name for info in pkgutil.walk_packages(pt.__path__, prefix="peptacular."))
    return sorted(names)


def _import_or_none(modname: str) -> types.ModuleType | None:
    """Import a module, or return None if an optional extra is not installed."""
    try:
        return importlib.import_module(modname)
    except ImportError:
        return None


def _collect_dataclasses() -> dict[str, type]:
    found: dict[str, type] = {}
    for modname in _module_names():
        module = _import_or_none(modname)
        if module is None:
            continue
        for obj in vars(module).values():
            if not isinstance(obj, type) or not dataclasses.is_dataclass(obj):
                continue
            if not obj.__module__.startswith("peptacular"):
                continue
            found[f"{obj.__module__}.{obj.__qualname__}"] = obj
    return found


DATACLASSES = _collect_dataclasses()

# A floor, not the exact count: a collection bug that imports almost nothing would
# otherwise leave this file passing while checking next to nothing.
MINIMUM_DATACLASSES = 30


def _is_descriptor(value: object) -> bool:
    cls = type(value)
    return any(hasattr(cls, attr) for attr in ("__get__", "__set__", "__delete__"))


def test_collection_finds_the_package_dataclasses() -> None:
    assert len(DATACLASSES) >= MINIMUM_DATACLASSES, sorted(DATACLASSES)


def _assert_no_descriptor_defaults(name: str, cls: type) -> None:
    for field in dataclasses.fields(cls):
        default = field.default
        if default is dataclasses.MISSING or isinstance(default, SAFE_DESCRIPTOR_TYPES):
            continue
        assert not isinstance(default, property), (
            f"{name}.{field.name} defaults to a property object. A base class most likely "
            f"declares {field.name} as a property, which dataclasses picked up as the default."
        )
        assert not isinstance(default, classmethod | staticmethod), f"{name}.{field.name} defaults to a {type(default).__name__}"
        assert not _is_descriptor(default), f"{name}.{field.name} defaults to a {type(default).__name__} descriptor"


@pytest.mark.parametrize("name", sorted(DATACLASSES))
def test_no_field_default_is_a_descriptor(name: str) -> None:
    _assert_no_descriptor_defaults(name, DATACLASSES[name])


def test_guard_detects_a_property_default() -> None:
    """The hazard this file exists for, rebuilt inline, must still be caught.

    Without this, a green run only shows that the walk found nothing. It does not show
    that the walk would notice the shape that put ``PositionScore`` at risk.
    """

    class _ScoreMixin:
        @property
        def position_id(self) -> str | None:
            return None

    @dataclasses.dataclass(frozen=True, slots=True)
    class _Scored(_ScoreMixin):
        position_id: str | None
        _: dataclasses.KW_ONLY
        score: float = 0.0

    # The hazard is silent: position_id is declared without a default, yet this
    # constructs, because dataclasses took the inherited property as the default.
    assert isinstance(_Scored().position_id, property)

    with pytest.raises(AssertionError, match="defaults to a property object"):
        _assert_no_descriptor_defaults("_Scored", _Scored)

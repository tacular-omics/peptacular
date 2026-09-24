"""Coverage for peptacular.utils: ModType helpers."""

import pytest

import peptacular as pt
from peptacular.constants import ModType
from peptacular.utils import _resolve_mod_types, get_mod_type


class TestGetModType:
    def test_modtype_passthrough(self):
        assert get_mod_type(ModType.INTERNAL) is ModType.INTERNAL

    def test_string_lookup(self):
        assert get_mod_type("internal") is ModType.INTERNAL

    def test_non_string_non_modtype_raises_type_error(self):
        with pytest.raises(TypeError, match="mod must be a string or ModType"):
            get_mod_type(5)  # type: ignore[arg-type]

    def test_unknown_string_raises_peptacular_error(self):
        with pytest.raises(pt.PeptacularError, match="Unknown mod type: 'bogus'"):
            get_mod_type("bogus")

    def test_exported_on_pt(self):
        assert pt.get_mod_type is get_mod_type


class TestResolveModTypes:
    def test_none_returns_all(self):
        assert _resolve_mod_types(None) == list(ModType)

    def test_single_modtype(self):
        assert _resolve_mod_types(ModType.INTERNAL) == [ModType.INTERNAL]

    def test_single_string(self):
        assert _resolve_mod_types("internal") == [ModType.INTERNAL]

    def test_iterable_of_mixed(self):
        assert _resolve_mod_types([ModType.INTERNAL, "static"]) == [ModType.INTERNAL, ModType.STATIC]

    def test_invalid_type_raises(self):
        with pytest.raises(TypeError, match="mods must be"):
            _resolve_mod_types(5)  # type: ignore[arg-type]

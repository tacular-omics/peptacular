"""Coverage for peptacular.utils: modification-string interning and ModType helpers."""

import pytest

from peptacular.constants import ModType
from peptacular.utils import get_mod_type, get_mods, handle_number_and_intern_mod


class TestHandleNumberAndInternMod:
    def test_int(self):
        assert handle_number_and_intern_mod(5) == "+5"

    def test_negative_float(self):
        assert handle_number_and_intern_mod(-5.5) == "-5.5"

    def test_string(self):
        assert handle_number_and_intern_mod("Oxidation") == "Oxidation"

    def test_string_is_stripped(self):
        assert handle_number_and_intern_mod("  Oxidation  ") == "Oxidation"

    def test_empty_string_raises(self):
        with pytest.raises(ValueError, match="Empty modification string"):
            handle_number_and_intern_mod("   ")


class TestGetModType:
    def test_modtype_passthrough(self):
        assert get_mod_type(ModType.INTERNAL) is ModType.INTERNAL

    def test_string_lookup(self):
        assert get_mod_type("internal") is ModType.INTERNAL

    def test_non_string_non_modtype_raises(self):
        with pytest.raises(ValueError, match="mod must be a string or ModType"):
            get_mod_type(5)  # type: ignore[arg-type]

    def test_unknown_string_raises(self):
        with pytest.raises(ValueError, match="Unknown mod type"):
            get_mod_type("bogus")


class TestGetMods:
    def test_none_returns_all(self):
        assert get_mods(None) == list(ModType)

    def test_single_modtype(self):
        assert get_mods(ModType.INTERNAL) == [ModType.INTERNAL]

    def test_single_string(self):
        assert get_mods("internal") == [ModType.INTERNAL]

    def test_iterable_of_mixed(self):
        assert get_mods([ModType.INTERNAL, "static"]) == [ModType.INTERNAL, ModType.STATIC]

    def test_invalid_type_raises(self):
        with pytest.raises(ValueError, match="mods parameter must be str"):
            get_mods(5)  # type: ignore[arg-type]

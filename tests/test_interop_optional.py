import pytest

from peptacular.interop import MissingOptionalDependencyError, _optional


def test_missing_optional_dependency_names_its_install_extra(monkeypatch):
    def missing(name):
        raise ModuleNotFoundError(name=name.partition(".")[0])

    monkeypatch.setattr(_optional, "import_module", missing)
    with pytest.raises(MissingOptionalDependencyError, match=r"peptacular\[alphabase\]"):
        _optional.require_dependency("alphabase.peptide.precursor", "alphabase")


def test_broken_transitive_dependency_is_not_reported_as_missing_extra(monkeypatch):
    def missing(name):
        raise ModuleNotFoundError(name="broken_transitive_dependency")

    monkeypatch.setattr(_optional, "import_module", missing)
    with pytest.raises(ModuleNotFoundError) as error:
        _optional.require_dependency("alphabase", "alphabase")
    assert error.value.name == "broken_transitive_dependency"

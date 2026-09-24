"""The top-level ``pt`` namespace is exactly ``pt.__all__`` (plus subpackages)."""

import ast
import importlib
import pkgutil
import types
from pathlib import Path

import pytest

import peptacular as pt
from peptacular.property import data as property_data
from peptacular.sequence import mod_builder


def _public_modules() -> list[str]:
    names = []
    for info in pkgutil.walk_packages(pt.__path__, prefix="peptacular."):
        parts = info.name.split(".")
        if any(part.startswith("_") for part in parts[1:]):
            continue
        if parts[1] in {"mcp", "interop"} and len(parts) > 2 and parts[2] not in {"__init__"}:
            # optional-extra modules import third-party packages lazily; checked where installed
            continue
        names.append(info.name)
    return names


def test_all_has_no_duplicates():
    assert len(pt.__all__) == len(set(pt.__all__))


def test_every_all_name_resolves():
    for name in pt.__all__:
        assert hasattr(pt, name), name


def test_pt_namespace_is_all_plus_subpackages():
    public = {n for n in dir(pt) if not n.startswith("_")}
    extra = {n for n in public - set(pt.__all__) if not isinstance(getattr(pt, n), types.ModuleType)}
    assert extra == set()


@pytest.mark.parametrize(
    "name",
    ["Any", "Counter", "Literal", "groupby", "UNIMOD_LOOKUP", "PROTEASE_LOOKUP", "ELEMENT_LOOKUP", "AminoAcid", "ElementInfo", "FRAGMENT_ION_LOOKUP"],
)
def test_stdlib_and_tacular_internals_do_not_leak(name):
    assert not hasattr(pt, name)


@pytest.mark.parametrize("name", ["IonType", "Protease", "NeutralDelta"])
def test_tacular_argument_enums_are_reexported(name):
    import tacular

    assert getattr(pt, name) is getattr(tacular, name)


@pytest.mark.parametrize(
    "name",
    [
        "parallelMethod",
        "parallelMethodLiteral",
        "FLIXIBILITY_SCALES",
        "get_regex_match_indices",
        "get_regex_match_range",
        "handle_number_and_intern_mod",
        "SupportsStr",
        "CV_TO_NAME_PREFIX",
        "ReadableProtocol",
        "MassPropertyMixin",
        "SEQUENCE_TYPE",
    ],
)
def test_removed_names_are_gone(name):
    assert not hasattr(pt, name)
    assert name not in pt.__all__


def test_flixibility_alias_removed_everywhere():
    for module in (pt, pt.property, property_data):
        with pytest.raises(AttributeError):
            module.FLIXIBILITY_SCALES  # noqa: B018


@pytest.mark.parametrize("modname", _public_modules())
def test_every_public_module_defines_all(modname):
    module = importlib.import_module(modname)
    assert isinstance(getattr(module, "__all__", None), list), modname
    for name in module.__all__:
        assert hasattr(module, name), f"{modname}.{name}"


def test_pt_get_mods_is_functional_api():
    assert pt.get_mods is mod_builder.get_mods
    assert pt.get_mods("PEM[Oxidation]TIDE") == pt.parse("PEM[Oxidation]TIDE").get_mods()


def test_flexibility_scales_spelling():
    assert pt.FLEXIBILITY_SCALES is property_data.FLEXIBILITY_SCALES
    assert pt.PhysicalPropertyScale.FLEXIBILITY_VIHINEN in pt.FLEXIBILITY_SCALES


def test_property_data_has_no_duplicate_module_constants():
    tree = ast.parse(Path(property_data.__file__).read_text(encoding="utf-8"))
    names = [
        target.id
        for node in tree.body
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        for target in (node.targets if isinstance(node, ast.Assign) else [node.target])
        if isinstance(target, ast.Name)
    ]
    duplicates = {name for name in names if names.count(name) > 1}
    assert not duplicates


def test_polarity_scales_are_enum_keyed():
    assert set(pt.POLARITY_SCALES) == set(pt.PolarityScale)
    assert pt.POLARITY_SCALES[pt.PolarityScale.GRANTHAM]["R"] == pytest.approx(10.5)


def test_parallel_method_enum():
    assert [m.value for m in pt.ParallelMethod] == ["process", "thread", "sequential"]
    assert pt.mass(["PEPTIDE", "PEP"], method=pt.ParallelMethod.SEQUENTIAL) == pt.mass(["PEPTIDE", "PEP"])

"""Regression tests for the top-level ``pt`` namespace and property-scale tables."""

import ast
import warnings
from pathlib import Path

import pytest

import peptacular as pt
from peptacular import chem, constants, regex_utils, spans, utils
from peptacular.property import data as property_data
from peptacular.sequence import mod_builder


def test_pt_get_mods_is_functional_api():
    # utils.get_mods (a ModType helper) used to shadow the documented sequence-level get_mods
    assert pt.get_mods is mod_builder.get_mods
    assert pt.get_mods("PEM[Oxidation]TIDE") == pt.parse("PEM[Oxidation]TIDE").get_mods()


@pytest.mark.parametrize("module", [chem, constants, regex_utils, spans, utils])
def test_star_import_modules_define_all(module):
    assert hasattr(module, "__all__")
    for name in module.__all__:
        assert hasattr(module, name), name


LEAKED_STDLIB_NAMES = ["Counter", "Mapping", "Sequence", "overload", "re", "sys", "warnings", "groupby", "NamedTuple", "Protocol", "StrEnum", "Final"]


@pytest.mark.parametrize("name", LEAKED_STDLIB_NAMES)
def test_stdlib_names_do_not_leak_into_pt(name):
    assert not hasattr(pt, name)


def test_flexibility_scales_spelling():
    assert pt.FLEXIBILITY_SCALES is property_data.FLEXIBILITY_SCALES
    assert pt.PhysicalPropertyScale.FLEXIBILITY_VIHINEN in pt.FLEXIBILITY_SCALES


@pytest.mark.parametrize("module", [pt, pt.property, property_data])
def test_flixibility_scales_is_deprecated_alias(module):
    with pytest.warns(DeprecationWarning, match="FLEXIBILITY_SCALES"):
        old = module.FLIXIBILITY_SCALES
    assert old is property_data.FLEXIBILITY_SCALES


def test_flixibility_scales_from_import_still_works():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        from peptacular.property import FLIXIBILITY_SCALES
    assert FLIXIBILITY_SCALES is property_data.FLEXIBILITY_SCALES


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

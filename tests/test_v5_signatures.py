"""5.0 signature convention: only the natural leading arguments are positional.

Rule (docs/migration.rst, "Keyword-only options"): the input (sequence, formula, ...) and any
required argument stay positional. Of the optional parameters, only the ones listed below stay
positional, because they read naturally in second place (``pt.mz(seq, 2)``). Every other
option, including all parallel options and every ``inplace``/``validate`` flag, is keyword-only.
"""

import inspect

import pytest

import peptacular as pt

PT_POSITIONAL_OPTIONS = {
    "mass": {"charge"},
    "mz": {"charge"},
    "comp": {"charge"},
    "isotopic_distribution": {"charge"},
    "frag": {"ion_type", "charge"},
    "fragment": {"ion_types", "charges"},
    "fast_fragment": {"ion_types", "charges"},
    "permutations": {"size"},
    "charge_at_ph": {"pH"},
    "set_start_method": {"method"},
    "diagnose": {"operation"},
    "filter_mods": {"mods"},
    "get_mods": {"mods"},
    "pop_mods": {"mods"},
    "remove_mods": {"mods"},
    "strip_mods": {"mods"},
    "generate_random": {"count"},
}

ANNOTATION_POSITIONAL_OPTIONS = {
    "mass": {"charge"},
    "mz": {"charge"},
    "comp": {"charge"},
    "isotopic_distribution": {"charge"},
    "estimate_isotopic_distribution": {"charge"},
    "frag": {"ion_type", "charge"},
    "fragment": {"ion_types", "charges"},
    "fast_fragment": {"ion_types", "charges"},
    "permutations": {"size"},
    "combinations": {"r"},
    "combinations_with_replacement": {"r"},
    "product": {"repeat"},
    "clear_mods": {"mods"},
    "filter_mods": {"mods"},
    "get_mods": {"mod_types"},
    "has_mods": {"mod_types"},
    "pop_mods": {"mod_types"},
}

PARALLEL_OPTIONS = {"n_workers", "chunksize", "method", "reuse_pool"}

PT_FUNCTIONS = sorted(name for name in pt.__all__ if inspect.isfunction(getattr(pt, name)))
ANNOTATION_METHODS = sorted(
    name for name, obj in vars(pt.ProFormaAnnotation).items() if not name.startswith("_") and inspect.isfunction(getattr(obj, "__func__", obj))
)
PARALLEL_FUNCTIONS = [name for name in PT_FUNCTIONS if "n_workers" in inspect.signature(getattr(pt, name)).parameters]


def _positional_options(func) -> set[str]:
    params = inspect.signature(func).parameters.values()
    return {p.name for p in params if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD) and p.default is not p.empty}


@pytest.mark.parametrize("name", PT_FUNCTIONS)
def test_pt_function_options_are_keyword_only(name):
    assert _positional_options(getattr(pt, name)) <= PT_POSITIONAL_OPTIONS.get(name, set())


@pytest.mark.parametrize("name", ANNOTATION_METHODS)
def test_annotation_method_options_are_keyword_only(name):
    obj = vars(pt.ProFormaAnnotation)[name]
    assert _positional_options(getattr(obj, "__func__", obj)) <= ANNOTATION_POSITIONAL_OPTIONS.get(name, set())


@pytest.mark.parametrize("name", ANNOTATION_METHODS)
def test_inplace_and_validate_are_keyword_only(name):
    obj = vars(pt.ProFormaAnnotation)[name]
    params = inspect.signature(getattr(obj, "__func__", obj)).parameters
    for flag in ("inplace", "validate"):
        if flag in params:
            assert params[flag].kind is inspect.Parameter.KEYWORD_ONLY, flag


def test_every_parallel_function_is_covered():
    assert {"mass", "mz", "comp", "digest", "fragment", "parse", "chem_comp", "chem_mass", "chem_formula", "parse_formula"} <= set(PARALLEL_FUNCTIONS)


@pytest.mark.parametrize("name", PARALLEL_FUNCTIONS)
def test_parallel_options_are_keyword_only(name):
    params = inspect.signature(getattr(pt, name)).parameters
    for option in PARALLEL_OPTIONS & params.keys():
        assert params[option].kind is inspect.Parameter.KEYWORD_ONLY, option


def test_positional_options_raise_type_error():
    with pytest.raises(TypeError):
        pt.mass(["PEPTIDE"], 0, "p")  # ty: ignore[too-many-positional-arguments]
    with pytest.raises(TypeError):
        pt.chem_mass("C2H6O", True)  # ty: ignore[too-many-positional-arguments]
    with pytest.raises(TypeError):
        pt.digest("PEPTIDEKAAR", "trypsin", 1)  # ty: ignore[too-many-positional-arguments]
    with pytest.raises(TypeError):
        pt.parse("PEPTIDE").reverse(True)  # ty: ignore[too-many-positional-arguments]


def test_charge_is_the_natural_second_argument():
    assert pt.mz("PEPTIDE", 2) == pt.mz("PEPTIDE", charge=2)
    assert pt.mass("PEPTIDE", 2) == pt.mass("PEPTIDE", charge=2)
    assert pt.parse("PEPTIDE").mz(2) == pt.mz("PEPTIDE", 2)
    assert pt.mass(["PEPTIDE"], n_workers=1) == [pt.mass("PEPTIDE")]

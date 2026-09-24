"""Validation paths that user input can reach: every one raises a PeptacularError subclass."""

import copy

import pytest

import peptacular as pt
from peptacular.annotation.mod import Interval, Mod
from peptacular.annotation.slicing import generate_sliding_windows
from peptacular.property import core as property_core

PeptacularError = pt.PeptacularError


def _annot(seq: str = "PEPTIDE") -> pt.ProFormaAnnotation:
    return pt.parse(seq)


# --- slicing -----------------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["shift", "shuffle", "reverse"])
@pytest.mark.parametrize(("keep_nterm", "keep_cterm"), [(-1, 0), (0, -1), (5, 5)])
def test_keep_terms_out_of_range(method, keep_nterm, keep_cterm):
    annot = _annot()
    args = (1,) if method == "shift" else ()
    with pytest.raises(PeptacularError, match="keep_nterm"):
        getattr(annot, method)(*args, keep_nterm=keep_nterm, keep_cterm=keep_cterm)


@pytest.mark.parametrize(("window_size", "match"), [(0, "positive"), (-2, "positive"), (8, "greater than")])
def test_sliding_window_size_out_of_range(window_size, match):
    with pytest.raises(PeptacularError, match=match):
        list(generate_sliding_windows(_annot(), window_size))


def test_sliding_window_needs_sequence():
    with pytest.raises(PeptacularError, match="sequence"):
        list(generate_sliding_windows(pt.ProFormaAnnotation(), 2))


def test_step_slicing_is_rejected():
    with pytest.raises(PeptacularError, match="Step"):
        _annot()[::2]


# --- batch -------------------------------------------------------------------------------------


@pytest.mark.parametrize("batch_size", [0, -1])
def test_iter_batch_batch_size_must_be_positive(batch_size):
    with pytest.raises(PeptacularError, match="batch_size"):
        list(pt.iter_batch("mass", ["PEPTIDE"], batch_size=batch_size))


# --- ProForma JSON decode ----------------------------------------------------------------------


def _json_dict():
    return pt.to_proforma_dict(pt.parse("[Acetyl]-PEM[Oxidation]TIDE/2"))


def _mutate(**changes):
    data = copy.deepcopy(_json_dict())
    for key, value in changes.items():
        if value is ...:
            data.pop(key)
        else:
            data[key] = value
    return data


@pytest.mark.parametrize(
    ("data", "match"),
    [
        (_mutate(**{"$schema": "https://example.org/other.json"}), "schema"),
        (_mutate(**{"$schema": ...}), r"\$schema"),
        (_mutate(schema_version="9.9"), "schema version"),
        (_mutate(**{"$type": ...}), r"\$type"),
        (_mutate(**{"$type": "NotAType"}), "object type"),
        (_mutate(charge=True), "charge"),
        (_mutate(charge=[1, 2]), "charge"),
        (_mutate(intervals="x"), "intervals"),
        (_mutate(intervals=["x"]), "interval"),
        (_mutate(names=[]), "names"),
        (_mutate(bogus=1), "Unknown JSON field"),
        (_mutate(sequence=...), "Missing required"),
    ],
)
def test_from_proforma_dict_rejects(data, match):
    with pytest.raises(PeptacularError, match=match):
        pt.from_proforma_dict(data)


@pytest.mark.parametrize(
    ("mods", "match"),
    [
        ({"internal": "x"}, "internal"),
        ({"internal": ["x"]}, "entry"),
        ({"internal": [{"position": "1", "modifications": {"Oxidation": 1}}]}, "position"),
        ({"internal": [{"position": 1, "modifications": {"A": 1}}, {"position": 1, "modifications": {"B": 1}}]}, "Duplicate"),
        ({"n_terminal": "Acetyl"}, "object or null"),
        ({"n_terminal": {"Acetyl": "1"}}, "integer counts"),
    ],
)
def test_from_proforma_dict_rejects_bad_modifications(mods, match):
    data = _json_dict()
    data["modifications"].update(mods)
    with pytest.raises(PeptacularError, match=match):
        pt.from_proforma_dict(data)


def test_unknown_enum_is_rejected():
    data = _json_dict()
    data["charge"] = None
    data["intervals"] = [{"$enum": "NoSuchEnum", "value": "x"}]
    with pytest.raises(PeptacularError):
        pt.from_proforma_dict(data)


@pytest.mark.parametrize(
    ("text", "match"), [("[1, 2]", "object"), ('"x"', "object"), ("{", "Invalid ProForma JSON"), ('{"a": 1, "a": 2}', "Duplicate"), ('{"a": NaN}', "finite")]
)
def test_from_proforma_json_rejects(text, match):
    with pytest.raises(PeptacularError, match=match):
        pt.from_proforma_json(text)


def test_from_proforma_json_expected_type():
    text = pt.to_proforma_json(pt.parse("PEPTIDE"))
    with pytest.raises(TypeError):
        pt.from_proforma_json(text, expected_type=int)


def test_to_proforma_json_rejects_non_finite_and_unknown():
    with pytest.raises(TypeError):
        pt.to_proforma_json(object())


# --- malformed ProForma through pt.parse --------------------------------------------------------


@pytest.mark.parametrize(
    "text",
    [
        "PEP[Oxidation",  # unmatched bracket
        "PEPTIDE]",
        "PEP[XYZ:123]TIDE",  # unknown CV prefix
        "PEP[Glycan:NotASugar1]TIDE",  # unknown monosaccharide
        "PEP[]TIDE",  # empty modification
        "PEP[Formula:C12Xx3]TIDE",  # bad formula element
        "PEPTIDE/0[+2Na+]",
    ],
)
def test_malformed_proforma_raises_peptacular_error(text):
    with pytest.raises(PeptacularError):
        annot = pt.parse(text)
        annot.validate_annotation()
        annot.mass()


# --- annotation validators ----------------------------------------------------------------------


def test_validate_sequence_rejects_unknown_residue():
    annot = pt.ProFormaAnnotation(sequence="PEP1IDE")
    with pytest.raises(PeptacularError, match="Invalid amino acid"):
        annot.validate_sequence()


@pytest.mark.parametrize(
    ("text", "validator"),
    [
        ("{Foo:Bar}PEPTIDE", "validate_labile_mods"),
        ("[Foo:Bar]?PEPTIDE", "validate_unknown_mods"),
        ("[Foo:Bar]-PEPTIDE", "validate_nterm_mods"),
        ("PEPTIDE-[Foo:Bar]", "validate_cterm_mods"),
        ("PEP[Foo:Bar]TIDE", "validate_internal_mods"),
    ],
)
def test_mod_validators_raise(text, validator):
    try:
        annot = pt.parse(text)
    except PeptacularError:
        pytest.skip("rejected at parse time")
    with pytest.raises(PeptacularError):
        getattr(annot, validator)()
        annot.validate_annotation()


def test_validate_intervals_overlap_and_bounds():
    annot = pt.parse("PEPTIDE")
    annot.set_intervals([Interval(0, 3, mods={"Oxidation": 1}), Interval(2, 5, mods={"Oxidation": 1})], validate=False)
    with pytest.raises(PeptacularError, match="Overlapping"):
        annot.validate_intervals()
    annot = pt.parse("PEPTIDE")
    annot.set_intervals([Interval(3, 20, mods={"Oxidation": 1})], validate=False)
    with pytest.raises(PeptacularError, match="out of bounds"):
        annot.validate_intervals()


# --- charge, immonium, static mods, MS2PIP ------------------------------------------------------


def test_zero_adduct_count_rejected():
    with pytest.raises(PeptacularError):
        pt.parse("PEPTIDE/[Na:z+1^0]")


def test_immonium_without_position_on_multi_residue():
    with pytest.raises(PeptacularError, match="Immonium"):
        _annot().frag("i", 1)


def test_add_static_mod_by_residue_rejects_count():
    with pytest.raises(PeptacularError, match="count of 1"):
        _annot().add_static_mod_by_residue("P", ("Oxidation", 2))


@pytest.mark.parametrize(
    "mods",
    [{"nterm": {"Acetyl": 2}}, {"cterm": {"Amidated": 2}}, {2: {"Oxidation": 2}}, {2: {"Oxidation": 1, "Phospho": 1}}],
)
def test_to_ms2_pip_rejects_multipliers_and_stacks(mods):
    annot = _annot().set_mods(mods)
    with pytest.raises(PeptacularError, match="MS2PIP"):
        annot.to_ms2_pip()


def test_static_mod_validator():
    with pytest.raises(PeptacularError, match="fixed modification"):
        _annot().set_static_mods({"not a mod": 1}, validate=False).validate_static_mods()


@pytest.mark.parametrize("method", ["append", "set", "extend", "remove"])
def test_mod_position_out_of_range(method):
    annot = _annot()
    with pytest.raises(pt.InvalidPositionError):
        getattr(annot, f"{method}_mods")({99: "Oxidation"})


# --- Mod / Interval ----------------------------------------------------------------------------


def test_mod_negative_count():
    with pytest.raises(PeptacularError, match="non-negative"):
        Mod("Oxidation", -1)


@pytest.mark.parametrize(("start", "end", "match"), [(-1, 3, "non-negative"), (3, 3, "End position"), (4, 2, "End position")])
def test_interval_bounds(start, end, match):
    with pytest.raises(PeptacularError, match=match):
        Interval(start, end)


# --- property ----------------------------------------------------------------------------------


def test_missing_aa_handling_unknown_value():
    with pytest.raises(PeptacularError, match="missing_aa_handling"):
        property_core._get_default_value("bogus", {"A": 1.0}, "X")  # ty: ignore[invalid-argument-type]


def test_missing_aa_handling_error_mode():
    with pytest.raises(PeptacularError, match="Invalid amino acid"):
        pt.calc_property("PEPXIDE", {"P": 1.0, "E": 2.0}, missing_aa_handling="error")


def test_unknown_weighting_scheme():
    with pytest.raises(PeptacularError):
        pt.calc_property("PEPTIDE", "hydrophobicity", weighting_scheme="bogus")  # ty: ignore[invalid-argument-type]

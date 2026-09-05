import json

import pytest
from pydantic import ValidationError

import peptacular as pt
from peptacular.mcp import contracts as c
from peptacular.mcp.operations import find_modifications, reference_rows, run_operation
from peptacular.mcp.outputs import ROWS


def calculate(name, annotations, **settings):
    inputs = [{"annotation": a} for a in annotations]
    request = c.SCIENTIFIC[name].model_validate({"inputs": inputs, **settings})
    records = [{"annotation": a, "id": "duplicate", "source_index": i, "source_key": f"input:{i}"} for i, a in enumerate(annotations)]
    result = run_operation(name, request.model_dump(), records)
    for row in result["records"]:
        ROWS[name].model_validate(row)
    json.dumps(result, allow_nan=False)
    return result


def test_inspection_does_not_resolve_unknown_chemistry():
    result = calculate("inspect_peptides", ["M[unknown modification]PEPTIDE"], detail="full")
    assert result["records"][0]["status"] == "complete"
    assert result["records"][0]["length"] == 8
    assert result["records"][0]["annotation"]["schema_version"] == "1.0"


def test_inspection_preserves_explicit_charge_carriers():
    row = calculate("inspect_peptides", ["PEPTIDE/[Na:z+1]"])["records"][0]
    assert row["encoded_charge"] == ["Na:z+1"]


def test_measurements_preserve_partial_success():
    result = calculate("analyze_peptides", ["PEP[+15.5]TIDE/2", "PEPTIDE", "INVALID!"], measurements=["neutral_mass_da", "mz", "composition"])
    good, uncharged, invalid = result["records"]
    assert good["neutral_mass_da"] == pytest.approx(pt.parse("PEPTIDE").neutral_mass() + 15.5)
    assert good["mz"] > 0
    assert good["composition"] is None
    assert good["diagnostics"][0]["code"] == "unavailable_composition"
    assert uncharged["diagnostics"][0]["code"] == "missing_charge"
    assert invalid["status"] == "error"


@pytest.mark.parametrize("charge", [-3, -1, 1, 2, 4])
def test_precursor_matches_core(charge):
    a = pt.parse("PEPTIDE").set_charge(charge)
    row = calculate("analyze_peptides", [a.serialize()], measurements=["mz", "ion_mass_da", "neutral_mass_da", "composition"])["records"][0]
    assert row["mz"] == a.mz()
    assert row["ion_mass_da"] == a.mass()
    assert row["neutral_mass_da"] == a.neutral_mass()
    assert row["external_charge"] + row["intrinsic_charge"] == row["charge"] == charge


def test_charge_conflicts_and_override():
    row = calculate("analyze_peptides", ["PEPTIDE/2"], charges=[3])["records"][0]
    assert row["diagnostics"][0]["code"] == "charge_conflict"
    rows = calculate("analyze_peptides", ["PEPTIDE/2"], charges=[2, 3], charge_policy="override")["records"]
    assert [r["charge"] for r in rows] == [2, 3]


@pytest.mark.parametrize("ion", ["a", "b", "c", "x", "y", "z", "p"])
@pytest.mark.parametrize("charge", [-1, 1, 2])
def test_fragment_values_and_spans(ion, charge):
    rows = calculate("fragment_peptides", ["PEPTIDE"], ion_series=[ion], charges=[charge], include=["composition", "sequence", "label"])["records"]
    expected = pt.parse("PEPTIDE").fragment(ion_types=[ion], charges=[charge])
    assert len(rows) == len(expected)
    for row, frag in zip(rows, expected, strict=True):
        assert row["status"] == "complete", row
        assert row["mz"] == frag.mz
        assert row["ion_mass_da"] == frag.mass
        assert row["end"] - row["start"] == (7 if ion == "p" else row["ordinal"])
        assert row["ordinal"] is None if ion == "p" else row["ordinal"] > 0


def test_numeric_fragment_delta_keeps_numbers():
    rows = calculate("fragment_peptides", ["PEPTIDE/2"], deltas=[{"kind": "mass", "value": -18.0}], include=["label", "composition"])["records"]
    assert all("mz" in row for row in rows)
    assert any(row["diagnostics"] for row in rows)


@pytest.mark.parametrize("axis", ["neutral_mass_da", "ion_mass_da", "mz", "neutron_offset"])
def test_isotope_axis(axis):
    a = pt.parse("PEPTIDE/2")
    rows = calculate("isotope_envelopes", [a.serialize()], axis=axis)["records"]
    effective = a.set_charge(0, inplace=False) if axis == "neutral_mass_da" else a
    core = effective.isotopic_distribution(use_neutron_count=axis == "neutron_offset")
    assert rows[0]["position"] == pytest.approx(core[0].mass / (2 if axis == "mz" else 1))
    assert max(row["relative_abundance"] for row in rows) == 1
    assert all(row["retained_probability"] is None for row in rows)


def test_properties():
    row = calculate("analyze_peptides", ["PEPTIDE"], measurements=["property"])["records"][0]
    assert row["property"] == pt.calc_property("PEPTIDE", "hphob_kyte_doolittle", method="sequential")


def test_comparison_delta():
    row = calculate(
        "compare_peptides", ["M[Oxidation]PEPTIDE"], reference={"annotation": "MPEPTIDE"}, measurements=["annotation", "neutral_mass_da", "composition"]
    )["records"][0]
    assert row["same_sequence"]
    assert row["composition_delta"]["O"] == 1
    assert row["neutral_mass_da"]["delta_input_minus_reference"] == pytest.approx(15.99491461957)


@pytest.mark.parametrize("specificity", ["full", "semi", "nonspecific"])
def test_digest_bounds(specificity):
    result = calculate("digest_proteins", ["AKPEPTIDERAAK"], specificity=specificity, min_length=2, max_length=8)
    for row in result["records"]:
        assert row["sequence"] == "AKPEPTIDERAAK"[row["start"] : row["end"]]
        assert 2 <= row["length"] <= 8


def test_unknown_enzyme_and_bounded_digest():
    row = calculate("digest_proteins", ["PEPTIDE"], enzyme=".*")["records"][0]
    assert row["diagnostics"][0]["code"] == "unknown_enzyme"
    result = calculate("digest_proteins", ["PEPTIDE"], specificity="nonspecific", max_rows=3)
    assert len(result["records"]) == 3
    assert result["computation"] == {"complete": False, "stop_reason": "row_limit"}


@pytest.mark.parametrize(
    "edit,expected",
    [
        ({"action": "add", "location": "internal", "index": 0, "modification": "Oxidation"}, "M[Oxidation]PEPTIDE"),
        ({"action": "add", "location": "nterm", "modification": "Acetyl"}, "[Acetyl]-MPEPTIDE"),
        ({"action": "charge", "charge": -2}, "MPEPTIDE/-2"),
        ({"action": "slice", "start": 1, "end": 4}, "PEP"),
    ],
)
def test_edits(edit, expected):
    row = calculate("edit_peptides", ["MPEPTIDE"], edits=[edit])["records"][0]
    assert row["proforma"] == expected
    assert row["original_proforma"] == "MPEPTIDE"


def test_edits_atomic_and_expand():
    row = calculate("edit_peptides", ["PEPTIDE"], edits=[{"action": "charge", "charge": 2}, {"action": "slice", "start": 1, "end": 99}])["records"][0]
    assert row["status"] == "error"
    assert "proforma" not in row
    row = calculate("edit_peptides", ["<[Carbamidomethyl]@C>AC"], edits=[{"action": "expand_static"}])["records"][0]
    assert row["proforma"] == "AC[Carbamidomethyl]"


def test_candidate_limit():
    result = calculate("enumerate_modifications", ["MMMM"], rules=[{"residues": "M", "modification": "Oxidation"}], max_candidates=2)
    assert len(result["records"]) == 2
    assert result["computation"]["complete"] is False
    assert result["computation"]["stop_reason"] == "candidate_limit"


def test_overlapping_mapping():
    inputs = [{"annotation": "AAA"}]
    request = c.Map(inputs=inputs, proteins=inputs)
    result = run_operation(
        "map_peptides",
        request.model_dump(),
        [{"annotation": "AAA", "source_key": "p:0", "source_index": 0}],
        [{"annotation": "AAAAA", "source_key": "protein:0", "source_index": 0}],
    )
    assert [r["start"] for r in result["records"]] == [0, 1, 2]
    assert all(r["match_count"] == 3 and r["ambiguous"] for r in result["records"])


@pytest.mark.parametrize("target", ["proforma", "stable_json"])
def test_conversion(target):
    row = calculate("convert_annotations", ["M[Oxidation]PEPTIDE/2"], target=target)["records"][0]
    assert row["status"] == "complete"
    if target == "stable_json":
        assert pt.ProFormaAnnotation.from_dict(row["value"]).serialize() == row["proforma"]


def test_reference_and_modification_search():
    rows = find_modifications(c.FindModifications(query_type="mass", query=15.9949, tolerance=0.001))
    assert any(row["name"] == "Oxidation" for row in rows)
    assert all(abs(row["mass_error_da"]) <= 0.001 for row in rows)
    assert reference_rows(c.GetReference(topic="enzymes"), {})
    assert reference_rows(c.GetReference(topic="scales"), {})


@pytest.mark.parametrize(
    "settings",
    [
        {"max_rows": True},
        {"max_rows": 0},
        {"charges": [True]},
        {"charges": [21]},
        {"extra": 1},
        {"monoisotopic": "false"},
        {"execution": {"timeout_seconds": 0}},
        {"inputs": []},
    ],
)
def test_reject_invalid_scientific_contract(settings):
    with pytest.raises(ValidationError):
        c.Analyze.model_validate({"inputs": [{"annotation": "PEPTIDE"}], **settings})


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, "15.99"])
def test_reject_invalid_mass_query(value):
    with pytest.raises(ValidationError):
        c.FindModifications(query_type="mass", query=value, tolerance=0.01)


def test_intrinsic_charge_mz_and_isotopes():
    text = "PEP[Formula:CH2:z+1]TIDE/2"
    row = calculate("analyze_peptides", [text], measurements=["mz"])["records"][0]
    assert (row["charge"], row["external_charge"], row["intrinsic_charge"]) == (3, 2, 1)
    peak = calculate("isotope_envelopes", [text], axis="mz")["records"][0]
    assert peak["charge"] == 3
    assert peak["position"] == pytest.approx(pt.parse(text).isotopic_distribution()[0].mass / 3)


@pytest.mark.parametrize("target,dependency", [("pyteomics", "pyteomics"), ("psm_utils", "psm_utils"), ("alphabase_row", "alphabase")])
def test_optional_conversion(target, dependency):
    pytest.importorskip(dependency)
    row = calculate("convert_annotations", ["M[Oxidation]PEPTIDE/2"], target=target)["records"][0]
    assert row["status"] == "complete", row
    assert row["value"]


@pytest.mark.parametrize(
    "settings",
    [
        {"query_type": "mass", "query": 0.0, "tolerance": 1.0, "tolerance_unit": "ppm"},
        {"query_type": "name", "query": "Oxidation", "tolerance": 1.0},
    ],
)
def test_mass_search_validation(settings):
    with pytest.raises(ValidationError):
        c.FindModifications(**settings)


@pytest.mark.parametrize("query_type,query,name_mode", [("accession", "35", "exact"), ("name", "Oxidation", "exact"), ("name", "Oxida", "prefix")])
def test_reference_name_modes(query_type, query, name_mode):
    result = find_modifications(c.FindModifications(query_type=query_type, query=query, name_mode=name_mode))
    assert any(r["name"] == "Oxidation" for r in result)


@pytest.mark.parametrize("topic", ["ions", "notation", "conventions", "schemas", "capabilities"])
def test_all_reference_topics(topic):
    assert reference_rows(c.GetReference(topic=topic), {})


def test_all_errors_still_respect_row_budget():
    result = calculate("analyze_peptides", ["BAD!", "BAD!", "BAD!"], max_rows=1)
    assert len(result["records"]) == 1
    assert not result["computation"]["complete"]
    assert result["records"][0]["diagnostics"][0]["stage"] == "parse"


def test_fragment_filters_and_precursor_intrinsic_neutralization():
    result = calculate("fragment_peptides", ["PEPTIDE/2"], min_mz=200.0, max_mz=400.0)
    assert all(200 <= row["mz"] <= 400 for row in result["records"])
    neutralized = calculate("analyze_peptides", ["PEP[Formula:CH2:z+1]TIDE/-1"], measurements=["mz"])["records"][0]
    assert neutralized["diagnostics"][0]["code"] == "missing_charge"


@pytest.mark.parametrize(
    "edit",
    [
        {"action": "remove", "location": "internal", "index": 0, "modification": "Oxidation"},
        {"action": "clear", "location": "internal", "index": 0},
    ],
)
def test_remove_edits(edit):
    row = calculate("edit_peptides", ["M[Oxidation]PEPTIDE"], edits=[edit])["records"][0]
    assert row["proforma"] == "MPEPTIDE"

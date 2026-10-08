"""MCP contracts accept exactly what the library accepts.

Each test sends, through the real server, an input the 5.0 contracts rejected although the
library calculates with it: U/O residues, non-a-z ion types, charge-carrier charges,
non-13C isotope labels and the RESID/GNOme vocabularies.
"""

import asyncio

import pytest

pytest.importorskip("mcp")

from pydantic import ValidationError  # noqa: E402
from tacular import AA_LOOKUP, IonType  # noqa: E402

import peptacular as pt  # noqa: E402
from peptacular.constants import CV  # noqa: E402
from peptacular.mcp import contracts as c  # noqa: E402
from peptacular.mcp.server import create_server  # noqa: E402


@pytest.fixture(scope="module")
def server():
    return create_server()


def call(server, name, request):
    result = asyncio.run(server.call_tool(name, {"request": request}))
    assert not result.is_error, result.structured_content
    return result.structured_content


def test_vocabularies_follow_library_cvs():
    resolvable = {cv for cv in CV if cv not in (CV.CUSTOM, CV.OBSERVED)}
    assert len(c.VOCABULARIES) == len(resolvable)
    assert c.FindModifications.model_json_schema()["properties"]["vocabularies"]["items"]["enum"] == list(c.VOCABULARIES)
    assert {"unimod", "psimod", "xlmod", "resid", "gnome"} == set(c.VOCABULARIES)


@pytest.mark.parametrize("vocabulary,accession", [("resid", "AA0002"), ("gnome", "G00008BG")])
def test_find_modifications_in_resid_and_gnome(server, vocabulary, accession):
    result = call(server, "find_modifications", {"query_type": "accession", "query": accession, "vocabularies": [vocabulary]})
    assert [(row["vocabulary"], row["accession"]) for row in result["records"]] == [(vocabulary, accession)]


def test_fragment_accepts_every_library_ion_type(server):
    assert c.Fragment.model_json_schema()["properties"]["ion_types"]["items"]["enum"] == [str(t) for t in IonType]
    request = {"inputs": [{"annotation": "PEPTIDEK/2"}], "ion_types": [str(t) for t in IonType], "include": ["sequence"], "max_rows": 5000}
    rows = call(server, "fragment_peptides", request)["records"]
    kinds = {row["ion_type"] for row in rows}
    assert {"z.", "c-H", "i", "by", "w"} <= kinds


@pytest.mark.parametrize("ion", ["z.", "c-H", "w", "i", "by", "cy"])
def test_widened_ion_spans_match_fragment_sequence(server, ion):
    text = "PEPTIDEK"
    request = {"inputs": [{"annotation": f"{text}/1"}], "ion_types": [ion], "include": ["sequence"]}
    rows = call(server, "fragment_peptides", request)["records"]
    expected = pt.parse(f"{text}/1").fragment(ion_types=[ion], charges=[1])
    assert len(rows) == len(expected) > 0
    for row, frag in zip(rows, expected, strict=True):
        assert row["mz"] == frag.mz
        assert text[row["start"] : row["end"]] == pt.parse(row["sequence"]).sequence
        assert row["position"] is None if isinstance(frag.position, tuple) else row["position"] == frag.position


def test_fragment_isotope_label_other_than_13c(server):
    request = {"inputs": [{"annotation": "PEPTIDE/1"}], "ion_types": ["y"], "isotopes": [0, {"15N": 1}]}
    rows = call(server, "fragment_peptides", request)["records"]
    expected = pt.parse("PEPTIDE/1").fragment(ion_types=["y"], charges=[1], isotopes=[0, {"15N": 1}])
    assert [row["mz"] for row in rows] == [f.mz for f in expected]
    assert {"15N": 1} in [row["isotopes"] for row in rows]


def test_fragment_rejects_unknown_isotope_label():
    with pytest.raises(ValidationError, match="Unknown element or isotope"):
        c.Fragment.model_validate({"inputs": [{"annotation": "PEPTIDE/1"}], "isotopes": [{"N15": 1}]})


def test_analyze_with_charge_carrier(server):
    rows = call(server, "analyze_peptides", {"inputs": [{"annotation": "PEPTIDE"}], "charges": ["Na:z+1"], "measurements": ["mz"]})["records"]
    assert rows[0]["mz"] == pytest.approx(pt.parse("PEPTIDE/[Na:z+1]").mz())


def test_carrier_charge_agrees_with_encoded_carrier(server):
    request = {"inputs": [{"annotation": "PEPTIDE/[Na:z+1,H:z+1]"}], "charges": [["Na:z+1", "H:z+1"]], "measurements": ["mz"]}
    rows = call(server, "analyze_peptides", request)["records"]
    assert rows[0]["mz"] == pytest.approx(pt.parse("PEPTIDE/[Na:z+1,H:z+1]").mz())


def test_fragment_keeps_encoded_carriers(server):
    rows = call(server, "fragment_peptides", {"inputs": [{"annotation": "PEPTIDE/[Na:z+1]"}], "ion_types": ["b"]})["records"]
    expected = pt.parse("PEPTIDE").fragment(ion_types=["b"], charges=["Na:z+1"])
    assert [row["mz"] for row in rows] == [f.mz for f in expected]


def test_isotopes_and_compare_and_edit_accept_carriers(server):
    peaks = call(server, "isotope_envelopes", {"inputs": [{"annotation": "PEPTIDE"}], "charges": ["Na:z+1"], "axis": "mz"})["records"]
    assert peaks[0]["charge"] == 1
    request = {"inputs": [{"annotation": "PEPTIDE"}], "reference": {"annotation": "PEPTIDE"}, "charges": ["Na:z+1"], "measurements": ["mz"]}
    compared = call(server, "compare_peptides", request)["records"]
    assert compared[0]["mz"]["delta_input_minus_reference"] == 0
    edited = call(server, "edit_peptides", {"inputs": [{"annotation": "PEPTIDE"}], "edits": [{"action": "charge", "charge": "Na:z+1"}]})["records"]
    assert edited[0]["proforma"] == "PEPTIDE/[Na:z+1]"


@pytest.mark.parametrize("charge", ["not a carrier", ["Na:z+1", "?"], []])
def test_invalid_charge_carrier_rejected(charge):
    with pytest.raises(ValidationError):
        c.Analyze.model_validate({"inputs": [{"annotation": "PEPTIDE"}], "charges": [charge]})


def test_enumerate_rules_accept_every_library_residue(server):
    assert set(c.RESIDUES) == {str(aa.id) for aa in AA_LOOKUP}
    rule = {"residues": "UO", "modification": "Oxidation"}
    rows = call(server, "enumerate_modifications", {"inputs": [{"annotation": "PEUTIOE"}], "rules": [rule], "max_variable_mods": 1})["records"]
    assert {row["proforma"] for row in rows} == {"PEUTIOE", "PEU[Oxidation]TIOE", "PEUTIO[Oxidation]E"}

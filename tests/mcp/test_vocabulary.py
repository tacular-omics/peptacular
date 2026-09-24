"""Drift guard: MCP names follow the library and the shared tacular-omics vocabulary.

Walks every tool's input and output schema, so a renamed library argument or a new
field that drifts from it fails here instead of reaching a model.
"""

import asyncio
import inspect
import re

import pytest

pytest.importorskip("mcp")

import peptacular as pt  # noqa: E402
from peptacular.mcp import contracts as c  # noqa: E402
from peptacular.mcp.server import create_server  # noqa: E402

# Shared tacular-omics bans, plus the peptacular 4.x MCP names replaced by library names in 5.0.
BANNED = re.compile(
    r"^(unit|tolerance_type|.*_tolerance_type|retention_time.*|inverse_reduced.*|ion_mobility_.*|target_mz|scan_start_time|ce|tic|TIC|time"
    r"|one_over_k0.*|mz_begin|mz_end|window_group|monoisotopic_mz"
    r"|ion_series|isotope_offsets|max_variable_modifications|ordinal|label|losses)$"
)
ALLOWED_UNIT_NAMES: set[str] = set()

# Row fields every scientific tool adds for provenance and status.
ROW_COMMON = {"source_id", "source_key", "source_index", "original_annotation", "row_key", "status", "diagnostics", "proforma"}
# Fragment row fields that are not library names: the residue span and the intrinsic part of the charge.
FRAGMENT_MCP_ONLY = {"start", "end", "intrinsic_charge"}
# Digest rows keep MCP names (protein spans for multi-step workflows); listed so any new key is a decision.
DIGEST_KEYS = {"sequence", "length", "start", "end", "missed_cleavages", "enzyme", "specificity"}
DIGEST_MCP_ONLY = {"protein_key", "protein_id", "protein_record_index", "protein_start", "protein_end"}

MINIMAL = {
    "get_reference": {},
    "find_modifications": {"query_type": "name", "query": "Phospho"},
    "inspect_peptides": {"inputs": [{"annotation": "PEPTIDE"}]},
    "analyze_peptides": {"inputs": [{"annotation": "PEPTIDE"}]},
    "fragment_peptides": {"inputs": [{"annotation": "PEPTIDE/2"}]},
    "compare_peptides": {"inputs": [{"annotation": "PEPTIDE"}], "reference": {"annotation": "PEPTIDE"}},
    "isotope_envelopes": {"inputs": [{"annotation": "PEPTIDE"}]},
    "digest_proteins": {"inputs": [{"annotation": "PEPTIDEKAAR"}]},
    "edit_peptides": {"inputs": [{"annotation": "PEPTIDE"}], "edits": [{"action": "charge", "charge": 2}]},
    "enumerate_modifications": {"inputs": [{"annotation": "PEPTMIDE"}], "rules": [{"residues": "M", "modification": "Oxidation"}]},
    "map_peptides": {"inputs": [{"annotation": "PEP"}], "proteins": [{"annotation": "APEPA"}]},
    "convert_annotations": {"inputs": [{"annotation": "PEPTIDE"}], "target": "proforma"},
}


def _walk(schema, path, out, defs, seen):
    if not isinstance(schema, dict):
        return
    ref = schema.get("$ref")
    if ref:
        name = ref.split("/")[-1]
        if name not in seen:
            seen.add(name)
            _walk(defs.get(name, {}), f"{path}<{name}>", out, defs, seen)
    for key, value in (schema.get("properties") or {}).items():
        out.append((f"{path}.{key}", key, value))
        _walk(value, f"{path}.{key}", out, defs, seen)
    for key in ("items", "anyOf", "oneOf", "allOf", "additionalProperties"):
        value = schema.get(key)
        for item in value if isinstance(value, list) else [value]:
            _walk(item, path, out, defs, seen)


def _enum_values(schema):
    values = set()
    for item in [schema, *schema.get("anyOf", [])]:
        values |= set(item.get("enum", []))
        if "const" in item:
            values.add(item["const"])
    return values


@pytest.fixture(scope="module")
def server():
    return create_server()


@pytest.fixture(scope="module")
def properties(server):
    rows = []
    for tool in asyncio.run(server.list_tools()):
        for kind, schema in (("in", tool.input_schema), ("out", tool.output_schema or {})):
            found = []
            _walk(schema, f"{tool.name}:{kind}", found, schema.get("$defs", {}), set())
            rows.extend(found)
    return rows


def call(server, name, arguments):
    return asyncio.run(server.call_tool(name, arguments))


def test_no_banned_names(properties):
    assert sorted(path for path, name, _ in properties if BANNED.match(name)) == []


def test_tolerance_switches(properties):
    for path, name, schema in properties:
        if name.endswith("unit") and name not in ALLOWED_UNIT_NAMES:
            assert re.fullmatch(r"(\w+_)?tolerance_unit", name), path
            assert _enum_values(schema) <= {"da", "ppm"}, path


@pytest.mark.parametrize("name", sorted(MINIMAL))
def test_unknown_argument_rejected(server, name):
    assert set(MINIMAL) == set(c.REQUESTS)
    assert not call(server, name, {"request": MINIMAL[name]}).is_error
    with pytest.raises(Exception, match="__bogus__"):
        call(server, name, {"request": MINIMAL[name], "__bogus__": 1})
    with pytest.raises(Exception, match="__bogus__"):
        call(server, name, {"request": {**MINIMAL[name], "__bogus__": 1}})


def test_fragment_request_uses_library_argument_names():
    fragment_args = {"ion_types", "isotopes", "deltas", "monoisotopic", "charges"}
    assert fragment_args <= set(c.Fragment.model_fields)
    assert fragment_args <= set(inspect.signature(pt.ProFormaAnnotation.fragment).parameters)
    assert "max_variable_mods" in c.Enumerate.model_fields
    assert "max_variable_mods" in inspect.signature(pt.ProFormaAnnotation.modify).parameters


def _fragment_attributes():
    """Public ``Fragment`` attributes plus the ``fragment_records`` keys."""
    return {name for name in dir(pt.Fragment) if not name.startswith("_")} | set(pt.FRAGMENT_RECORD_KEYS)


def test_fragment_row_keys_match_library(server):
    request = {"inputs": [{"annotation": "PEPTIDE/2"}], "include": ["composition", "sequence", "mzpaf"], "max_rows": 3}
    result = call(server, "fragment_peptides", {"request": request})
    assert result.structured_content["contract_version"] == "2.0"
    rows = result.structured_content["records"]
    assert rows
    library = _fragment_attributes()
    for row in rows:
        assert set(row) - ROW_COMMON - FRAGMENT_MCP_ONLY <= library, set(row) - ROW_COMMON - FRAGMENT_MCP_ONLY - library
        for key in ("ion_type", "position", "charge_state", "mz", "mass", "neutral_mass", "mzpaf"):
            assert key in row
    from peptacular.mcp.outputs import FragmentRow

    declared = set(FragmentRow.model_fields) - ROW_COMMON - FRAGMENT_MCP_ONLY - {"sequence"}
    assert declared <= library, declared - library


def test_digest_row_keys_are_listed(server):
    result = call(server, "digest_proteins", {"request": MINIMAL["digest_proteins"]})
    for row in result.structured_content["records"]:
        assert set(row) - ROW_COMMON == DIGEST_KEYS | DIGEST_MCP_ONLY


def test_reference_ions_use_ion_type(server):
    result = call(server, "get_reference", {"request": {"topic": "ions"}})
    assert all(set(row) == {"ion_type", "kind"} for row in result.structured_content["records"])

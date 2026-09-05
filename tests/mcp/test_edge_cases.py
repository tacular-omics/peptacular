import json

import pytest
from pydantic import ValidationError

from peptacular.mcp import contracts as c
from peptacular.mcp.cli import main
from peptacular.mcp.server import call_tool


@pytest.mark.parametrize(
    "model,values",
    [
        (c.Delta, {"kind": "mass", "value": "H2O"}),
        (c.Delta, {"kind": "formula", "value": 1.0}),
        (c.Digest, {"min_length": 10, "max_length": 5}),
        (c.Fragment, {"min_mz": 500.0, "max_mz": 100.0}),
        (c.ModificationEdit, {"action": "add", "location": "internal", "modification": "Oxidation"}),
        (c.ModificationEdit, {"action": "clear", "location": "nterm", "modification": "Oxidation"}),
        (c.Rule, {"modification": "Oxidation"}),
    ],
)
def test_contradictory_requests(model, values):
    if issubclass(model, c.Scientific):
        values["inputs"] = [{"annotation": "PEPTIDE"}]
    with pytest.raises(ValidationError):
        model.model_validate(values)


@pytest.mark.parametrize(
    "field,value",
    [
        ("execution", {"mode": "job"}),
        ("preflight", True),
        ("inputs", {"kind": "reference", "reference_id": "old_result"}),
        ("inputs", {"kind": "inline", "records": [{"annotation": "PEPTIDE"}]}),
    ],
)
def test_removed_workflow_arguments_rejected(field, value):
    with pytest.raises(ValidationError):
        c.Analyze.model_validate({"inputs": [{"annotation": "PEPTIDE"}], field: value})


def test_input_budget():
    result = call_tool("inspect_peptides", c.Inspect(inputs=[{"annotation": "A" * 10000}] * 20))
    assert result.status == "error"
    assert result.diagnostics[0].code == "input_limit"


def test_result_byte_budget():
    result = call_tool("digest_proteins", c.Digest(inputs=[{"annotation": "A" * 1000}], specificity="nonspecific", max_rows=5000))
    assert result.status == "partial"
    assert result.computation.stop_reason == "byte_limit"
    assert len(json.dumps(result.records).encode()) < 250000
    assert result.total_rows is None


def test_cli_check_and_removed_workspace_options(capsys):
    assert main(["--check"]) == 0
    assert len(json.loads(capsys.readouterr().out)["tools"]) == 12
    with pytest.raises(SystemExit):
        main(["--workspace", "/unused"])


def test_unexpected_errors_do_not_leak_inputs(monkeypatch):
    import peptacular.mcp.server as server

    def broken(*args, **kwargs):
        raise RuntimeError("private sequence")

    monkeypatch.setattr(server, "run_operation", broken)
    result = call_tool("analyze_peptides", c.Analyze(inputs=[{"annotation": "PEPTIDE"}]))
    assert result.status == "error"
    assert result.diagnostics[0].code == "internal_error"
    assert "private sequence" not in result.model_dump_json()

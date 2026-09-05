import asyncio
import json
import time
from pathlib import Path

import pytest
from pydantic import ValidationError

from peptacular.mcp import contracts as c
from peptacular.mcp.cli import main
from peptacular.mcp.config import Config
from peptacular.mcp.operations import ServiceError
from peptacular.mcp.server import Service
from peptacular.mcp.storage import Store, matches


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
        (c.Dataset, {}),
        (c.Query, {"result_id": "x", "aggregate": "min"}),
        (c.Query, {"result_id": "x", "aggregate": "count", "create_view": True}),
    ],
)
def test_contradictory_requests(model, values):
    if issubclass(model, c.Scientific):
        values["inputs"] = {"kind": "inline", "records": [{"annotation": "PEPTIDE"}]}
    with pytest.raises(ValidationError):
        model.model_validate(values)


@pytest.mark.parametrize(
    "text,format",
    [
        ("PEPTIDE", "fasta"),
        (">\nPEPTIDE", "fasta"),
        ("", "fasta"),
        ("a,a\nx,y\n", "csv"),
        ("other\nPEPTIDE\n", "csv"),
        ("[]\n", "jsonl"),
    ],
)
def test_invalid_files(store, tmp_path, text, format):
    path = tmp_path / "bad"
    path.write_text(text)
    with pytest.raises(ValueError):
        store.register(c.Dataset(path=str(path), format=format))


@pytest.mark.parametrize(
    "operator,value,expected",
    [
        ("eq", 2.0, True),
        ("ne", 3.0, True),
        ("lt", 3.0, True),
        ("le", 2.0, True),
        ("gt", 1.0, True),
        ("ge", 2.0, True),
        ("in", [2.0, 3.0], True),
        ("is_null", False, True),
    ],
)
def test_filter_operators(operator, value, expected):
    assert matches(2.0, c.Filter(column="x", operator=operator, value=value)) == expected


@pytest.mark.parametrize("operator,value", [("in", 2.0), ("is_null", None), ("lt", "a"), ("eq", [1.0])])
def test_invalid_filter_values(operator, value):
    with pytest.raises(ValueError):
        matches(2, c.Filter(column="x", operator=operator, value=value))


def test_projection_invalid_columns_and_large_rows(store):
    result = store.save("result", "large", [{"proforma": "PEPTIDE", "detail": "x" * 500, "nested": {"x": 1}}])
    with pytest.raises(ServiceError, match="page"):
        store.query(c.Query(result_id=result), byte_limit=100)
    assert store.query(c.Query(result_id=result, columns=["proforma"]), byte_limit=100)["records"] == [{"proforma": "PEPTIDE"}]
    with pytest.raises(ServiceError, match="Unknown columns"):
        store.query(c.Query(result_id=result, columns=["missing"]))
    with pytest.raises(ValueError):
        store.query(c.Query(result_id=result, sort=[{"column": "nested"}]))
    with pytest.raises(ValueError):
        store.query(c.Query(result_id=result, aggregate="min", aggregate_column="proforma"))
    with pytest.raises(ValueError):
        store.query(c.Query(result_id=result, aggregate="group_count", aggregate_column="nested"))


def test_cleanup_managed_exports_preserves_user_exports(store, tmp_path):
    source = store.save("result", "test", [{"sequence": "PEPTIDE", "source_id": "=SUM(1,2)"}])
    managed = store.export(c.Export(result_id=source))
    user = store.export(c.Export(result_id=source, destination="user.csv"))
    assert "'=SUM" in Path(user["path"]).read_text()
    with store.connect() as db:
        db.execute("UPDATE objects SET expires = ?", (time.time() - 1,))
    assert store.clean() == 3
    assert not Path(managed["path"]).exists()
    assert Path(user["path"]).exists()


def test_cache_workspace_boundary(config, tmp_path):
    Store(config)
    other = tmp_path / "other"
    other.mkdir()
    with pytest.raises(ValueError, match="another workspace"):
        Store(Config(other, cache=config.cache))


def test_cli_check_cleanup_invalid_and_default_cache(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "application-cache"))
    assert main(["--workspace", str(tmp_path), "--check"]) == 0
    assert json.loads(capsys.readouterr().out)["versions"]["mcp"]
    assert main(["--workspace", str(tmp_path), "cache", "clean"]) == 0
    assert "removed_expired_objects" in capsys.readouterr().out
    assert main(["--workspace", str(tmp_path / "missing"), "--check"]) == 2
    assert "Peptacular MCP" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main(["--workspace", str(tmp_path), "cache"])


@pytest.mark.parametrize("kwargs", [{"workers": 0}, {"storage_bytes": True}])
def test_invalid_config(tmp_path, kwargs):
    with pytest.raises(ValueError):
        Config(tmp_path, **kwargs)


@pytest.mark.asyncio
async def test_preflight_and_graceful_queued_shutdown(config):
    service = Service(config)
    request = c.Analyze(inputs={"kind": "inline", "records": [{"annotation": "PEPTIDE"}]}, preflight=True)
    response = await service.call("analyze_peptides", request)
    assert response.records[0]["input_records"] == 1
    request = request.model_copy(update={"preflight": False, "execution": c.Execution(mode="job")})
    for _ in range(config.workers):
        await service.runner.semaphore.acquire()
    job = await service.call("analyze_peptides", request)
    await asyncio.sleep(0)
    await service.runner.close()
    assert service.store.get(job.job_id)["data"]["state"] == "interrupted"
    response = await service.call("analyze_peptides", request)
    assert response.diagnostics[0].code == "server_stopping"


def test_duplicate_json_fields_rejected(store, tmp_path):
    path = tmp_path / "duplicate.jsonl"
    path.write_text('{"proforma":"PEPTIDE","proforma":"OTHER"}\n')
    with pytest.raises(ValueError, match="Duplicate JSON"):
        store.register(c.Dataset(path=str(path), format="jsonl"))


@pytest.mark.parametrize("value,expected", [(2, True), (2.0, True), (True, False)])
def test_numeric_membership_keeps_booleans_distinct(value, expected):
    clause = c.Filter(column="charge", operator="in", value=[2.0, 1.0])
    assert matches(value, clause) is expected

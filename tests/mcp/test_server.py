import asyncio
import json
import multiprocessing
import os
import subprocess
import sys
import time

import pytest
from mcp import Client
from mcp.client.stdio import StdioServerParameters

from peptacular.mcp import contracts as c
from peptacular.mcp.server import Service, create_server

INPUT = {"kind": "inline", "records": [{"id": "example", "annotation": "PEPTIDE/2"}]}


@pytest.mark.asyncio
async def test_sdk_schemas_resources_and_calculation(config):
    async with Client(create_server(config)) as client:
        listed = await client.list_tools()
        assert {tool.name for tool in listed.tools} == set(c.REQUESTS)
        total_schema_bytes = len(json.dumps([tool.model_dump() for tool in listed.tools]))
        assert total_schema_bytes < 180000
        for tool in listed.tools:
            assert tool.input_schema["type"] == "object"
            assert tool.output_schema["type"] == "object"
        result = await client.call_tool("analyze_peptides", {"request": {"inputs": INPUT, "measurements": ["mz"]}})
        assert not result.is_error
        assert result.structured_content["records"][0]["mz"] > 0
        resources = await client.read_resource("peptacular://conventions")
        assert "Spectacular" in str(resources)
        invalid = await client.call_tool("analyze_peptides", {"request": {"inputs": INPUT, "unexpected": True}})
        assert invalid.is_error


@pytest.mark.asyncio
@pytest.mark.parametrize("mode", ["auto", "legacy"])
async def test_stdio_subprocess(config, mode):
    parameters = StdioServerParameters(
        command=sys.executable, args=["-m", "peptacular.mcp", "--workspace", str(config.workspace), "--cache", str(config.cache)]
    )
    async with Client(parameters, mode=mode) as client:
        tools = await client.list_tools()
        assert len(tools.tools) == 18
        result = await client.call_tool("inspect_peptides", {"request": {"inputs": INPUT}})
        assert result.structured_content["records"][0]["sequence"] == "PEPTIDE"
        assert not result.is_error


@pytest.mark.asyncio
async def test_dataset_digest_analyze_view_export_workflow(config):
    service = Service(config)
    try:
        config.workspace.joinpath("proteins.fasta").write_text(">same\nAKPEPTIDERAAK\n>same\nMPEPTIDERAAK\n")
        dataset = await service.call("register_dataset", c.Dataset(path="proteins.fasta"))
        source = c.Reference(kind="reference", reference_id=dataset.records[0]["dataset_id"])
        digest = await service.call("digest_proteins", c.Digest(inputs=source, min_length=3))
        assert digest.status == "complete", digest
        analysis = await service.call(
            "analyze_peptides", c.Analyze(inputs=c.Reference(kind="reference", reference_id=digest.result_id), charges=[2, 3], measurements=["mz", "length"])
        )
        assert analysis.status == "complete", analysis
        assert len({r["source_key"] for r in analysis.records}) == len(digest.records)
        assert all("start" in row and "parent_source_key" in row for row in analysis.records)
        view = await service.call(
            "query_result", c.Query(result_id=analysis.result_id, filters=[{"column": "mz", "operator": "ge", "value": 300.0}], create_view=True)
        )
        exported = await service.call("export_result", c.Export(result_id=view.result_id, destination="precursors.csv"))
        assert exported.status == "complete", exported
        assert config.workspace.joinpath("precursors.csv").read_text().count("\n") == view.page.total_rows + 1
    finally:
        await service.runner.close()


@pytest.mark.asyncio
async def test_job_success_idempotency_and_restart(config):
    service = Service(config)
    request = c.Analyze(inputs=INPUT, execution={"mode": "job", "idempotency_key": "one"})
    try:
        job = await service.call("analyze_peptides", request)
        repeated = await service.call("analyze_peptides", request)
        assert job.job_id == repeated.job_id
        async with asyncio.timeout(15):
            while True:
                state = await service.call("get_job", c.Job(job_id=job.job_id))
                if state.status not in ("queued", "running"):
                    break
                await asyncio.sleep(0.05)
        assert state.status == "complete", state
        assert state.records[0]["progress"]["records_consumed"] == 1
        second = Service(config)
        try:
            result = await second.call("query_result", c.Query(result_id=state.result_id))
            assert result.records[0]["length"] == 7
        finally:
            await second.runner.close()
    finally:
        await service.runner.close()


def blocking_worker(name, payload, records, proteins, output):
    time.sleep(30)


@pytest.mark.asyncio
async def test_cancellation_timeout_and_replacement(config, monkeypatch):
    import peptacular.mcp.execution as execution

    service = Service(config)
    original = execution.worker_main
    monkeypatch.setattr(execution, "worker_main", blocking_worker)
    try:
        job = await service.call("analyze_peptides", c.Analyze(inputs=INPUT, execution={"mode": "job"}))
        async with asyncio.timeout(5):
            while not service.runner.active:
                await asyncio.sleep(0.01)
        result = await service.call("cancel_job", c.Job(job_id=job.job_id))
        assert result.status == "cancelled"
        assert not service.runner.active
        assert not multiprocessing.active_children()
        result = await service.call("analyze_peptides", c.Analyze(inputs=INPUT, execution={"mode": "inline", "timeout_seconds": 1}))
        assert result.diagnostics[0].code == "timeout"
        assert not multiprocessing.active_children()
        monkeypatch.setattr(execution, "worker_main", original)
        result = await service.call("analyze_peptides", c.Analyze(inputs=INPUT))
        assert result.status == "complete", result
    finally:
        await service.runner.close()


@pytest.mark.asyncio
async def test_foreign_job_cannot_be_cancelled(config):
    service, other = Service(config), Service(config)
    try:
        job = service.store.save("job", "foreign", {"state": "running"}, {"owner": service.runner.owner})
        state = await other.call("get_job", c.Job(job_id=job))
        assert state.status == "interrupted"
        assert state.records[0]["state"] == "unavailable"
        cancel = await other.call("cancel_job", c.Job(job_id=job))
        assert cancel.diagnostics[0].code == "job_owner"
        assert service.store.get(job)["data"]["state"] == "running"
    finally:
        await service.runner.close()
        await other.runner.close()


def test_cli_check_creates_no_cache(config):
    result = subprocess.run(
        [sys.executable, "-m", "peptacular.mcp", "--workspace", str(config.workspace), "--cache", str(config.cache), "--check"],
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["ok"]
    assert not config.cache.exists()


def test_import_has_no_optional_side_effects(tmp_path):
    code = "\n".join(
        [
            "import sys",
            "import peptacular",
            "import peptacular.mcp",
            "import peptacular.mcp.cli",
            "assert 'mcp' not in sys.modules",
            "assert 'pydantic' not in sys.modules",
        ]
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, timeout=15, env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")}
    )
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "cache").exists()


@pytest.mark.asyncio
async def test_stdio_disconnect_stops_owned_job(config):
    from peptacular.mcp.storage import Store

    parameters = StdioServerParameters(
        command=sys.executable,
        args=[
            "-m",
            "peptacular.mcp",
            "--workspace",
            str(config.workspace),
            "--cache",
            str(config.cache),
        ],
    )
    async with Client(parameters) as client:
        result = await client.call_tool(
            "enumerate_modifications",
            {
                "request": {
                    "inputs": {"kind": "inline", "records": [{"annotation": "M" * 100}]},
                    "rules": [{"residues": "M", "modification": "Oxidation"}],
                    "execution": {"mode": "job"},
                }
            },
        )
        job_id = result.structured_content["job_id"]
    store = Store(config)
    assert store.get(job_id)["data"]["state"] in ("interrupted", "succeeded", "partially_succeeded")
    assert not list(config.cache.glob("worker_*"))


@pytest.mark.asyncio
async def test_job_field_failure_is_not_truncated_computation(config):
    service = Service(config)
    try:
        request = c.Analyze(
            inputs={"kind": "inline", "records": [{"annotation": "PEP[+15.5]TIDE/2"}]}, measurements=["mz", "composition"], execution={"mode": "job"}
        )
        job = await service.call("analyze_peptides", request)
        async with asyncio.timeout(15):
            while True:
                state = await service.call("get_job", c.Job(job_id=job.job_id))
                if state.status not in ("queued", "running"):
                    break
                await asyncio.sleep(0.05)
        assert state.status == "partial"
        assert state.computation.complete
        assert state.computation.stop_reason is None
    finally:
        await service.runner.close()

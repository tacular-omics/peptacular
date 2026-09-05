import json
import multiprocessing
import os
import subprocess
import sys

import pytest
from mcp import Client
from mcp.client.stdio import StdioServerParameters

from peptacular.mcp import contracts as c
from peptacular.mcp.server import call_tool, create_server

INPUT = [{"id": "example", "annotation": "PEPTIDE/2"}]


@pytest.mark.asyncio
async def test_sdk_schemas_resources_and_calculation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    async with Client(create_server()) as client:
        listed = await client.list_tools()
        assert {tool.name for tool in listed.tools} == set(c.REQUESTS)
        assert len(listed.tools) == 12
        for tool in listed.tools:
            assert tool.input_schema["type"] == "object"
            assert tool.output_schema["type"] == "object"
            assert tool.annotations.read_only_hint
            assert not tool.annotations.destructive_hint
        result = await client.call_tool("analyze_peptides", {"request": {"inputs": INPUT, "measurements": ["mz"]}})
        assert not result.is_error
        assert result.structured_content["records"][0]["mz"] > 0
        assert json.loads(result.content[0].text) == result.structured_content
        assert not {"result_id", "job_id", "page"} & result.structured_content.keys()
        resources = await client.read_resource("peptacular://conventions")
        assert "Spectacular" in str(resources)
        schema = await client.read_resource("peptacular://schemas/analyze_peptides")
        assert "measurements" in str(schema)
        invalid = await client.call_tool("analyze_peptides", {"request": {"inputs": INPUT, "unexpected": True}})
        assert invalid.is_error
        invalid_outer = await client.call_tool("analyze_peptides", {"request": {"inputs": INPUT}, "extra": True})
        assert invalid_outer.is_error
    assert not list(tmp_path.iterdir())
    assert not multiprocessing.active_children()


@pytest.mark.asyncio
@pytest.mark.parametrize("mode", ["auto", "legacy"])
async def test_stdio_subprocess(mode):
    parameters = StdioServerParameters(command=sys.executable, args=["-m", "peptacular.mcp"])
    async with Client(parameters, mode=mode) as client:
        tools = await client.list_tools()
        assert len(tools.tools) == 12
        result = await client.call_tool("inspect_peptides", {"request": {"inputs": INPUT}})
        assert result.structured_content["records"][0]["sequence"] == "PEPTIDE"
        assert not result.is_error


@pytest.mark.asyncio
async def test_sdk_diverts_dependency_stdout():
    code = "\n".join(
        [
            "import os",
            "import peptacular.mcp.server as server",
            "original = server.run_operation",
            "def noisy(*args, **kwargs):",
            "    print('dependency print', flush=True)",
            "    os.write(1, b'native dependency output\\n')",
            "    return original(*args, **kwargs)",
            "server.run_operation = noisy",
            "server.create_server().run(transport='stdio')",
        ]
    )
    parameters = StdioServerParameters(command=sys.executable, args=["-c", code])
    async with Client(parameters) as client:
        response = await client.call_tool("analyze_peptides", {"request": {"inputs": INPUT}})
        assert response.structured_content["status"] == "complete"


def test_inline_digest_to_analysis_preserves_caller_ids():
    digest = call_tool(
        "digest_proteins",
        c.Digest(
            inputs=[
                {"id": "same", "annotation": "AKPEPTIDERAAK"},
                {"id": "same", "annotation": "MPEPTIDERAAK"},
            ],
            min_length=3,
        ),
    )
    assert digest.status == "complete"
    assert {r["source_index"] for r in digest.records} == {0, 1}
    inputs = [{"id": row["row_key"], "annotation": row["proforma"]} for row in digest.records]
    analysis = call_tool("analyze_peptides", c.Analyze(inputs=inputs, charges=[2, 3], measurements=["mz", "length"]))
    assert analysis.status == "complete"
    assert len(analysis.records) == len(digest.records) * 2
    assert {row["source_id"] for row in analysis.records} == {row["row_key"] for row in digest.records}


def test_direct_calls_use_no_processes(monkeypatch):
    def unexpected_process(*args, **kwargs):
        raise AssertionError("A scalar calculation should not create a process")

    monkeypatch.setattr(multiprocessing, "get_context", unexpected_process)
    result = call_tool("analyze_peptides", c.Analyze(inputs=INPUT))
    assert result.status == "complete"
    assert result.returned_rows == result.total_rows == 1


def test_bounded_inline_response_has_no_hidden_results():
    response = call_tool("digest_proteins", c.Digest(inputs=[{"annotation": "PEPTIDE"}], specificity="nonspecific", max_rows=3))
    assert len(response.records) == response.returned_rows == 3
    assert response.status == "partial"
    assert response.computation.stop_reason == "row_limit"
    assert response.total_rows is None
    assert response.next_offset is None


def test_partial_measurement_does_not_imply_truncation():
    response = call_tool("analyze_peptides", c.Analyze(inputs=[{"annotation": "PEP[+15.5]TIDE/2"}], measurements=["mz", "composition"]))
    assert response.status == "partial"
    assert response.computation.complete
    assert response.records[0]["mz"] > 0
    assert response.records[0]["composition"] is None


def test_reference_lookup_pages():
    first = call_tool("get_reference", c.GetReference(topic="enzymes", limit=2))
    second = call_tool("get_reference", c.GetReference(topic="enzymes", limit=2, offset=first.next_offset))
    assert first.returned_rows == 2
    assert first.records != second.records
    assert first.total_rows == second.total_rows
    mods = call_tool("find_modifications", c.FindModifications(query_type="mass", query=15.9949, tolerance=0.001))
    assert any(row["name"] == "Oxidation" for row in mods.records)


def test_cli_check_creates_no_files(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "peptacular.mcp", "--check"],
        capture_output=True,
        text=True,
        timeout=15,
        cwd=tmp_path,
        env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / "cache")},
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["ok"]
    assert not list(tmp_path.iterdir())


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

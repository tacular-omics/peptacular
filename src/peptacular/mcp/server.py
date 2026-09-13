"""Stateless MCP tools using the SDK's ordinary synchronous tool execution."""

import json
import logging
import uuid
from typing import Annotated

from mcp.server import MCPServer
from mcp.server.mcpserver.tools.base import Tool
from mcp.types import CallToolResult, TextContent, ToolAnnotations

from . import contracts as c
from .operations import CONVENTIONS, LIMITS, ServiceError, diagnostic, find_modifications, reference_rows, run_operation, versions
from .outputs import OUTPUTS

DESCRIPTIONS = {
    "get_reference": "Discover supported tools, installed extras, limits, enzyme IDs, property scales, notation and request schemas. Start here.",
    "inspect_peptides": "Explain ProForma or stable JSON without mass calculation. Reports names, modifications, ambiguity and encoded charge.",
    "analyze_peptides": (
        "Calculate selected theoretical precursor mass, m/z, composition or sequence properties. Keeps per-measurement successes and diagnostics."
    ),
    "fragment_peptides": (
        "Generate theoretical backbone and precursor ions. Supply or encode a charge. Use max_rows and ion filters to bound the response. No observed spectra."
    ),
    "compare_peptides": (
        "Compare each input with one explicit reference annotation and theoretical precursor values. Deltas are input minus reference. No spectral matching."
    ),
    "isotope_envelopes": (
        "Calculate aggregated theoretical isotope distributions on an explicit mass, m/z or neutron axis. Reports normalization and truncation settings."
    ),
    "digest_proteins": (
        "Digest proteins using a known enzyme ID and bounded full, semi or nonspecific settings. Keeps source IDs and zero-based end-exclusive spans."
    ),
    "find_modifications": (
        "Find reference modifications by accession, name or mass with explicit Da/ppm tolerance. Returns candidates, never assigns a modification."
    ),
    "edit_peptides": (
        "Apply explicit modifications, charge edits, slicing or fixed-mod expansion atomically to each input copy. Internal edit indexes are zero-based."
    ),
    "enumerate_modifications": (
        "Enumerate bounded structural candidates from residue/terminal rules. Candidates are not probabilities. Inspect computation.complete for limits."
    ),
    "map_peptides": "Map peptide sequences to proteins, retaining overlapping and repeated matches. Modified inputs require an explicit ignore policy.",
    "convert_annotations": (
        "Return portable ProForma, stable JSON, AlphaBase rows or compatibility-checked Pyteomics/psm_utils text. Optional targets require their extras."
    ),
}


def input_records(records, prefix):
    return [
        {"annotation": record.annotation, "id": record.id, "source_index": index, "source_key": f"{prefix}:{index}"} for index, record in enumerate(records)
    ]


def call_tool(name, request):
    """Run one bounded request and return its complete inline response."""
    request_id = uuid.uuid4().hex
    try:
        if len(request.model_dump_json().encode()) > LIMITS["request_bytes"]:
            raise ServiceError("request_limit", "Tool arguments exceed 1 MiB. Split the batch into smaller requests.")
        if name in c.SCIENTIFIC:
            records = input_records(request.inputs, "input")
            proteins = input_records(request.proteins, "protein") if name == "map_peptides" else None
            annotations = [row["annotation"] for row in records + (proteins or [])]
            if name == "compare_peptides":
                annotations.append(request.reference.annotation)
            if sum(len(json.dumps(value)) for value in annotations) > LIMITS["annotation_characters"]:
                raise ServiceError("input_limit", "The batch exceeds 100,000 annotation characters. Split it into smaller requests.")
            result = run_operation(name, request.model_dump(), records, proteins)
            rows = result["records"]
            partial = not result["computation"]["complete"] or any(row.get("diagnostics") for row in rows)
            return c.Envelope(
                request_id=request_id,
                status="partial" if partial else "complete",
                records=rows,
                returned_rows=len(rows),
                total_rows=len(rows) if result["computation"]["complete"] else None,
                diagnostics=result["diagnostics"],
                computation=result["computation"],
                applied_settings={**request.model_dump(exclude={"inputs", "proteins", "reference"}), "versions": versions()},
            )
        rows = reference_rows(request, LIMITS) if name == "get_reference" else find_modifications(request)
        selected = rows[request.offset : request.offset + request.limit]
        next_offset = request.offset + len(selected)
        return c.Envelope(
            request_id=request_id,
            records=selected,
            returned_rows=len(selected),
            total_rows=len(rows),
            next_offset=next_offset if next_offset < len(rows) else None,
            applied_settings=request.model_dump(),
        )
    except (ValueError, KeyError, OSError) as error:
        exc = error
    except Exception as error:
        logging.getLogger(__name__).error("Request %s failed in %s with %s", request_id, name, type(error).__name__)
        exc = ServiceError("internal_error", f"Unexpected failure for request {request_id}. Check server logs.")
    return c.Envelope(
        request_id=request_id,
        status="error",
        diagnostics=[c.Diagnostic.model_validate(diagnostic(exc))],
        computation=c.Computation(complete=False, stop_reason=getattr(exc, "code", "invalid_request")),
    )


def create_server():
    tools = []

    def register(name, model):
        output_model = OUTPUTS.get(name, c.Envelope)

        def invoke(request):
            envelope = call_tool(name, request)
            data = envelope.model_dump(mode="json")
            if envelope.status != "error":
                output_model.model_validate(data)
            # Both content forms carry the same bounded records. Nothing is hidden in a cache.
            return CallToolResult(
                content=[TextContent(type="text", text=json.dumps(data, ensure_ascii=False, separators=(",", ":")))],
                structured_content=data,
                is_error=envelope.status == "error",
            )

        invoke.__annotations__ = {"request": model, "return": Annotated[CallToolResult, output_model]}
        tool = Tool.from_function(
            invoke,
            name=name,
            description=DESCRIPTIONS[name],
            structured_output=True,
            annotations=ToolAnnotations(read_only_hint=True, destructive_hint=False, idempotent_hint=True, open_world_hint=False),
        )
        # SDK argument models otherwise ignore unknown outer fields. Keep this compatibility
        # adjustment local and tested against the narrowly pinned SDK release line.
        tool.fn_metadata.arg_model.model_config["extra"] = "forbid"
        tool.fn_metadata.arg_model.model_rebuild(force=True)
        tool.parameters = tool.fn_metadata.arg_model.model_json_schema()
        tools.append(tool)

    for name, model in c.REQUESTS.items():
        register(name, model)

    server = MCPServer(
        "Peptacular",
        version="1.0",
        tools=tools,
        instructions="Use structured tools for theoretical peptide and protein calculations. Arguments are under request. "
        "Provide inputs as a list of annotation records. Each call returns results directly and retains no data. "
        "For truncated results, narrow the request or split the batch. Discover conventions with get_reference. "
        "Spectacular owns observed spectra and spectrum matching.",
    )

    @server.resource("peptacular://conventions", mime_type="application/json")
    def conventions():
        return json.dumps(CONVENTIONS)

    @server.resource("peptacular://schemas/{tool}", mime_type="application/json")
    def schema(tool: str):
        if tool not in c.REQUESTS:
            raise ValueError("Unknown tool")
        return json.dumps(c.REQUESTS[tool].model_json_schema())

    return server

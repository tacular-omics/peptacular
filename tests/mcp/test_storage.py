import gzip
import json
import time
from pathlib import Path

import pytest

from peptacular.mcp import contracts as c
from peptacular.mcp.config import Config
from peptacular.mcp.operations import ServiceError
from peptacular.mcp.storage import Store


def test_snapshot_duplicates_and_reference_keys(store, tmp_path):
    path = tmp_path / "input.fasta.gz"
    path.write_bytes(gzip.compress(b">duplicate\nAKPEPTIDER\n>duplicate\nPEPTIDE\n"))
    dataset = store.register(c.Dataset(path=str(path)))
    path.write_bytes(b"changed")
    resolved = store.resolve_inputs(c.Reference(kind="reference", reference_id=dataset["dataset_id"]))
    assert [r["annotation"] for r in resolved] == ["AKPEPTIDER", "PEPTIDE"]
    assert resolved[0]["id"] == resolved[1]["id"] == "duplicate"
    assert resolved[0]["source_key"] != resolved[1]["source_key"]


@pytest.mark.parametrize(
    "format,text",
    [
        ("csv", "id,peptide\na,PEPTIDE\n"),
        ("tsv", "id\tpeptide\na\tPEPTIDE\n"),
        ("jsonl", '{"id":"a","peptide":"PEPTIDE"}\n'),
    ],
)
def test_table_inputs(store, tmp_path, format, text):
    path = tmp_path / f"input.{format}"
    path.write_text(text)
    result = store.register(c.Dataset(path=str(path), format=format, annotation_column="peptide", id_column="id"))
    rows, _ = store.rows(result["dataset_id"])
    assert rows[0]["id"] == "a"
    assert rows[0]["proforma"] == "PEPTIDE"


def test_symlink_and_traversal_boundaries(store, tmp_path):
    outside = tmp_path.parent / f"outside-{tmp_path.name}.fasta"
    outside.write_text(">a\nPEPTIDE\n")
    link = tmp_path / "link.fasta"
    try:
        link.symlink_to(outside)
    except OSError:
        pytest.skip("Symlink creation is unavailable")
    with pytest.raises(ValueError):
        store.register(c.Dataset(path=str(link)))
    with pytest.raises(ValueError):
        store.config.destination("../outside.csv")


def test_compressed_byte_limit(tmp_path):
    path = tmp_path / "large.gz"
    path.write_bytes(gzip.compress(b">a\n" + b"A" * 2000))
    store = Store(Config(tmp_path, cache=tmp_path / "cache", input_bytes=1000))
    with pytest.raises(ServiceError, match="Decompressed"):
        store.register(c.Dataset(path=str(path)))


def test_query_views_cursors_aggregates_and_restart(store, config):
    result = store.save("result", "test", [{"proforma": "PEPTIDE", "mz": i, "group": i % 2} for i in range(10)])
    query = c.Query(result_id=result, filters=[{"column": "mz", "operator": "ge", "value": 3.0}], sort=[{"column": "mz", "descending": True}], limit=2)
    first = store.query(query)
    assert [r["mz"] for r in first["records"]] == [9, 8]
    second = store.query(query.model_copy(update={"cursor": first["page"]["next_cursor"]}))
    assert [r["mz"] for r in second["records"]] == [7, 6]
    with pytest.raises(ValueError, match="cursor"):
        store.query(c.Query(result_id=result, cursor=first["page"]["next_cursor"]))
    view = store.query(query.model_copy(update={"create_view": True}))["view_id"]
    assert len(Store(config).rows(view)[0]) == 7
    assert store.query(c.Query(result_id=view, aggregate="count"))["records"] == [{"count": 7}]
    assert store.query(c.Query(result_id=view, aggregate="min", aggregate_column="mz"))["records"] == [{"min": 3}]
    grouped = store.query(c.Query(result_id=result, aggregate="group_count", aggregate_column="group"))["records"]
    assert grouped == [{"value": 0, "count": 5}, {"value": 1, "count": 5}]


@pytest.mark.parametrize("format", ["json", "jsonl", "csv", "fasta"])
def test_exports_are_atomic_and_idempotent(store, format):
    result = store.save("result", "test", [{"sequence": "PEPTIDE", "source_id": "a", "mz": 500}])
    request = c.Export(result_id=result, format=format, destination=f"output.{format}", idempotency_key="key")
    first = store.export(request)
    assert Path(first["path"]).is_file()
    assert store.export(request) == first
    with pytest.raises(ServiceError, match="different"):
        store.export(request.model_copy(update={"destination": f"other.{format}"}))
    with pytest.raises(FileExistsError):
        store.export(request.model_copy(update={"idempotency_key": None}))


def test_failed_fasta_export_leaves_original(store, tmp_path):
    result = store.save("result", "test", [{"sequence": "M[Oxidation]"}])
    destination = tmp_path / "original.fasta"
    destination.write_text("keep")
    with pytest.raises(ValueError):
        store.export(c.Export(result_id=result, format="fasta", destination=destination.name, overwrite=True))
    assert destination.read_text() == "keep"


def test_ttl_and_quota(tmp_path):
    store = Store(Config(tmp_path, cache=tmp_path / "cache", storage_bytes=1000))
    result = store.save("result", "test", [{"value": 1}])
    with store.connect() as db:
        db.execute("UPDATE objects SET expires = ? WHERE id = ?", (time.time() - 1, result))
    with pytest.raises(ServiceError, match="expired"):
        store.get(result)
    assert store.clean() == 1
    with pytest.raises(ServiceError, match="quota"):
        store.save("dataset", "large", ["x" * 1001])


def test_numeric_reference_column_rejected(store):
    result = store.save("result", "numbers", [{"mz": 123.0}])
    with pytest.raises(ServiceError, match="annotation text"):
        store.resolve_inputs(c.Reference(kind="reference", reference_id=result, column="mz"))


def test_stable_json_dataset(store, tmp_path):
    import peptacular as pt

    path = tmp_path / "input.json"
    path.write_text(json.dumps([pt.parse("PEPTIDE").to_dict()]))
    dataset = store.register(c.Dataset(path=str(path), format="stable_json"))
    assert dataset["record_count"] == 1

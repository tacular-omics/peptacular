"""Workspace snapshots, finite queries, bounded ingestion, and atomic exports."""

import base64
import csv
import gzip
import hashlib
import io
import json
import os
import sqlite3
import tempfile
import threading
import time
import uuid
from contextlib import contextmanager
from pathlib import Path

from .contracts import Record
from .operations import ServiceError


def encode(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))


def fingerprint(value):
    return hashlib.sha256(encode(value).encode()).hexdigest()


class BoundedReader(io.RawIOBase):
    def __init__(self, stream, limit):
        self.stream = stream
        self.remaining = limit

    def readable(self):
        return True

    def readinto(self, buffer):
        chunk = self.stream.read(min(len(buffer), self.remaining + 1))
        self.remaining -= len(chunk)
        if self.remaining < 0:
            raise ServiceError("input_limit", "Decompressed input exceeds the byte budget.")
        buffer[: len(chunk)] = chunk
        return len(chunk)


def csv_value(value):
    if isinstance(value, (dict, list)):
        return encode(value)
    if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
        return "'" + value
    return value


def unique_json_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON field: {key}")
        result[key] = value
    return result


class Store:
    def __init__(self, config):
        self.config = config
        config.cache.mkdir(parents=True, exist_ok=True, mode=0o700)
        self.lock = threading.RLock()
        self.path = config.cache / "workspace.sqlite3"
        with self.connect() as db:
            db.execute("CREATE TABLE IF NOT EXISTS settings (workspace TEXT NOT NULL)")
            row = db.execute("SELECT workspace FROM settings").fetchone()
            if row and row[0] != str(config.workspace):
                raise ValueError("Cache belongs to another workspace")
            if not row:
                db.execute("INSERT INTO settings VALUES (?)", (str(config.workspace),))
            db.execute("""CREATE TABLE IF NOT EXISTS objects (
                id TEXT PRIMARY KEY, kind TEXT NOT NULL, name TEXT NOT NULL, created REAL NOT NULL,
                expires REAL NOT NULL, metadata TEXT NOT NULL, data TEXT NOT NULL, bytes INTEGER NOT NULL,
                idem TEXT UNIQUE, request_hash TEXT)""")

    @contextmanager
    def connect(self):
        connection = sqlite3.connect(self.path, timeout=30)
        connection.row_factory = sqlite3.Row
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def save(self, kind, name, data, metadata=None, *, idem=None, request_hash=None, object_id=None):
        encoded = encode(data)
        meta = encode(metadata or {})
        size = len(encoded.encode()) + len(meta.encode()) + (metadata or {}).get("managed_file_bytes", 0)
        now = time.time()
        identifier = object_id or f"{kind}_{uuid.uuid4().hex}"
        with self.lock, self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            if idem:
                prior = db.execute("SELECT * FROM objects WHERE idem = ?", (idem,)).fetchone()
                if prior:
                    if prior["request_hash"] != request_hash:
                        raise ServiceError("idempotency_conflict", "This idempotency key was used for a different request.")
                    if prior["expires"] <= now:
                        raise ServiceError("expired_reference", "The idempotent result has expired. Use a new key.")
                    return prior["id"]
            used = db.execute("SELECT COALESCE(SUM(bytes),0) FROM objects").fetchone()[0]
            if used + size > self.config.storage_bytes:
                raise ServiceError("storage_limit", "Workspace storage quota reached. Run the scoped cache clean command or increase the quota.")
            db.execute(
                "INSERT INTO objects VALUES (?,?,?,?,?,?,?,?,?,?)",
                (identifier, kind, name, now, now + self.config.ttl_seconds, meta, encoded, size, idem, request_hash),
            )
        return identifier

    def get(self, identifier):
        with self.connect() as db:
            row = db.execute("SELECT * FROM objects WHERE id = ?", (identifier,)).fetchone()
        if row is None:
            raise ServiceError("unknown_reference", "No managed object has this ID in the configured workspace.")
        if row["expires"] <= time.time():
            raise ServiceError("expired_reference", "This managed object has expired. Register the input again or use a retained export.")
        result = dict(row)
        result["data"] = json.loads(result["data"])
        result["metadata"] = json.loads(result["metadata"])
        return result

    def idempotent(self, key, request_hash):
        if not key:
            return None
        with self.connect() as db:
            row = db.execute("SELECT id,request_hash FROM objects WHERE idem = ?", (key,)).fetchone()
        if row:
            if row["request_hash"] != request_hash:
                raise ServiceError("idempotency_conflict", "This key belongs to a different request.")
            return self.get(row["id"])
        return None

    def update_job(self, identifier, owner, data):
        encoded = encode(data)
        with self.lock, self.connect() as db:
            row = db.execute("SELECT metadata FROM objects WHERE id = ? AND kind = 'job'", (identifier,)).fetchone()
            if row is None or json.loads(row[0]).get("owner") != owner:
                raise ServiceError("job_owner", "Only the owning server can change a live job.")
            db.execute("UPDATE objects SET data = ?, bytes = ? WHERE id = ?", (encoded, len(encoded.encode()), identifier))

    def list(self, request):
        with self.connect() as db:
            rows = db.execute("SELECT id,kind,name,created,expires,metadata FROM objects ORDER BY created DESC,id").fetchall()
        result = []
        for raw in rows:
            row = dict(raw)
            row["metadata"] = json.loads(row["metadata"])
            row["expired"] = row["expires"] <= time.time()
            if request.kind and row["kind"] != request.kind:
                continue
            if request.search.casefold() not in encode(row).casefold():
                continue
            result.append(row)
        return result

    def clean(self):
        # Original inputs, active jobs, and user exports are never deleted.
        with self.lock, self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            expired = db.execute("SELECT id,data,metadata FROM objects WHERE expires <= ? AND kind = 'export'", (time.time(),)).fetchall()
            for row in expired:
                if json.loads(row["metadata"]).get("managed"):
                    data = json.loads(row["data"])
                    path = Path(data["path"])
                    expected = self.config.cache / f"{row['id']}.{data['format']}"
                    if path == expected and not path.is_symlink():
                        path.unlink(missing_ok=True)
            cursor = db.execute("DELETE FROM objects WHERE expires <= ? AND kind != 'job'", (time.time(),))
            count = cursor.rowcount
            terminal = db.execute(
                "DELETE FROM objects WHERE expires <= ? AND kind = 'job' AND json_extract(data,'$.state') IN "
                "('succeeded','partially_succeeded','failed','cancelled','interrupted')",
                (time.time(),),
            )
            count += terminal.rowcount
        return count

    def rows(self, identifier):
        obj = self.get(identifier)
        if obj["kind"] not in ("dataset", "result", "view"):
            raise ServiceError("invalid_reference_type", "Expected a dataset, result, or view reference.")
        return obj["data"], obj

    def resolve_inputs(self, inputs):
        if inputs.kind == "inline":
            key = "inline_" + uuid.uuid4().hex
            rows = [{"annotation": r.annotation, "id": r.id, "source_index": i, "source_key": f"{key}:{i}"} for i, r in enumerate(inputs.records)]
        else:
            data, obj = self.rows(inputs.reference_id)
            rows = []
            for i, row in enumerate(data):
                value = row.get(inputs.column)
                if not isinstance(value, (str, dict)):
                    raise ServiceError("invalid_input_column", f"Column {inputs.column!r} must contain annotation text or stable JSON on every row.")
                rows.append(
                    {
                        "annotation": value,
                        "id": row.get("id", row.get("source_id")),
                        "source_index": i,
                        "source_key": f"{obj['id']}:{i}",
                        "context": {
                            **{
                                k: v
                                for k, v in row.items()
                                if k
                                in (
                                    "protein_id",
                                    "protein_key",
                                    "protein_start",
                                    "protein_record_index",
                                    "protein_end",
                                    "start",
                                    "end",
                                    "missed_cleavages",
                                    "enzyme",
                                    "specificity",
                                )
                            },
                            "parent_reference_id": obj["id"],
                            "parent_row_key": row.get("row_key"),
                            "parent_source_key": row.get("source_key"),
                        },
                    }
                )
        size = len(encode(rows).encode())
        residues = sum(len(encode(row["annotation"])) for row in rows)
        if size > self.config.input_bytes or len(rows) > self.config.max_records or residues > self.config.max_residues:
            raise ServiceError("input_limit", "Input exceeds the configured byte, record, or residue budget.")
        if any(len(encode(row["annotation"])) > self.config.max_sequence_length for row in rows):
            raise ServiceError("sequence_limit", "An annotation exceeds the configured length budget.")
        return rows

    def register(self, request):
        if request.records is not None:
            records = [record.model_dump() for record in request.records]
            original_path = None
        else:
            path = self.config.input_path(request.path)
            original_path = str(path)
            if path.stat().st_size > self.config.input_bytes:
                raise ServiceError("input_limit", "Compressed or source file exceeds the byte budget.")
            with gzip.open(path, "rb") if path.suffix.lower() == ".gz" else path.open("rb") as stream:
                with io.TextIOWrapper(io.BufferedReader(BoundedReader(stream, self.config.input_bytes)), encoding="utf-8-sig") as text:
                    records = self._parse_file(text, request)
        if len(records) > self.config.max_records:
            raise ServiceError("input_limit", "Dataset exceeds the record budget.")
        rows = []
        for index, record in enumerate(records):
            header = record.pop("_header", None)
            validated = Record.model_validate(record)
            rows.append({"id": validated.id, "proforma": validated.annotation, "record_index": index, "row_key": f"row_{index}"})
            if header is not None:
                rows[-1]["header"] = header
        if not rows:
            raise ServiceError("empty_dataset", "Input contains no records.")
        if sum(len(encode(row["proforma"])) for row in rows) > self.config.max_residues:
            raise ServiceError("input_limit", "Dataset exceeds the annotation character budget.")
        identifier = self.save("dataset", request.name, rows, {"source_path": original_path, "format": request.format, "record_count": len(rows)})
        return {"dataset_id": identifier, "record_count": len(rows), "annotation_column": "proforma", "snapshot": True}

    def _parse_file(self, text, request):
        if request.format == "fasta":
            records = []
            name, sequence = None, []
            for line in text:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        records.append({"id": name.split()[0], "annotation": "".join(sequence), "_header": name})
                    name, sequence = line[1:], []
                    if not name:
                        raise ValueError("FASTA headers must not be empty")
                elif name is None:
                    raise ValueError("FASTA sequence precedes the first header")
                else:
                    sequence.append(line)
                if len(records) > self.config.max_records:
                    raise ServiceError("input_limit", "Too many FASTA records.")
            if name is not None:
                records.append({"id": name.split()[0], "annotation": "".join(sequence), "_header": name})
            return records
        if request.format == "stable_json":
            data = json.load(text, object_pairs_hook=unique_json_object)
            return [{"annotation": item} for item in (data if isinstance(data, list) else [data])]
        if request.format in ("csv", "tsv"):
            reader = csv.DictReader(text, delimiter="," if request.format == "csv" else "\t")
            if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError("Table must have unique column headers")
        else:
            reader = (json.loads(line, object_pairs_hook=unique_json_object) for line in text if line.strip())
        result = []
        for row in reader:
            if not isinstance(row, dict) or request.annotation_column not in row:
                raise ValueError("The declared annotation column is missing")
            if request.id_column and request.id_column not in row:
                raise ValueError("The declared ID column is missing")
            result.append({"id": row[request.id_column] if request.id_column else None, "annotation": row[request.annotation_column]})
            if len(result) > self.config.max_records:
                raise ServiceError("input_limit", "Too many table records.")
        return result

    def query(self, request, *, byte_limit=None):
        rows, obj = self.rows(request.result_id)
        available = set().union(*(row.keys() for row in rows)) if rows else set()
        requested = set(request.columns or []) | {f.column for f in request.filters} | {s.column for s in request.sort}
        for clause in request.filters:
            # Validate even when the result is empty or every row is filtered out.
            if clause.operator == "in" and not isinstance(clause.value, list):
                raise ValueError("in requires a list value")
            if clause.operator == "is_null" and not isinstance(clause.value, bool):
                raise ValueError("is_null requires a boolean value")
        if request.aggregate_column:
            requested.add(request.aggregate_column)
        if requested - available:
            raise ServiceError("unknown_column", f"Unknown columns: {', '.join(sorted(requested - available))}")
        rows = [row for row in rows if all(matches(row.get(f.column), f) for f in request.filters)]
        for sort in reversed(request.sort):
            values = [row.get(sort.column) for row in rows if row.get(sort.column) is not None]
            if any(isinstance(v, (dict, list)) for v in values) or len({type(v) for v in values}) > 1:
                raise ValueError("Sort columns must contain consistently typed scalar values")
            rows.sort(key=lambda row: (row.get(sort.column) is None, row.get(sort.column)), reverse=sort.descending)
        if request.aggregate:
            values = [row.get(request.aggregate_column) for row in rows if row.get(request.aggregate_column) is not None]
            if request.aggregate == "count":
                rows = [{"count": len(rows)}]
            elif request.aggregate == "group_count":
                counts = {}
                for row in rows:
                    value = row.get(request.aggregate_column)
                    if isinstance(value, (dict, list)):
                        raise ValueError("Group counts require scalar values")
                    key = encode(value)
                    counts[key] = counts.get(key, 0) + 1
                rows = [{"value": json.loads(key), "count": value} for key, value in sorted(counts.items())]
            else:
                if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in values):
                    raise ValueError("Min/max require a numeric column")
                rows = [{request.aggregate: (min(values) if request.aggregate == "min" else max(values)) if values else None}]
        elif request.columns:
            rows = [{column: row.get(column) for column in request.columns} for row in rows]
        shape = request.model_dump(exclude={"cursor", "limit", "create_view"})
        query_hash = fingerprint(shape)
        offset = 0
        if request.cursor:
            try:
                cursor = json.loads(base64.urlsafe_b64decode(request.cursor))
                if cursor["query"] != query_hash or type(cursor["offset"]) is not int or cursor["offset"] < 0:
                    raise ValueError("Cursor does not match this query")
                offset = cursor["offset"]
            except (ValueError, KeyError, TypeError) as exc:
                raise ValueError("Invalid or mismatched query cursor") from exc
        page, next_cursor = page_rows(rows, offset, request.limit, byte_limit or self.config.page_bytes, query_hash)
        view = None
        if request.create_view:
            view = self.save(
                "view", f"View of {obj['name']}", rows, {"parent_id": obj["id"], "query": shape, "computation": obj["metadata"].get("computation")}
            )
        return {
            "records": page,
            "page": {"returned_rows": len(page), "total_rows": len(rows), "next_cursor": next_cursor},
            "view_id": view,
            "computation": obj["metadata"].get("computation") or {"complete": True, "stop_reason": None},
        }

    def export(self, request):
        request_hash = fingerprint(request.model_dump(exclude={"idempotency_key"}))
        key = f"export:{request.idempotency_key}" if request.idempotency_key else None
        with self.lock:
            prior = self.idempotent(key, request_hash)
            if prior:
                if not Path(prior["data"]["path"]).is_file():
                    raise ServiceError("export_missing", "The earlier export file is no longer present.")
                return prior["data"]
            rows, obj = self.rows(request.result_id)
            identifier = "export_" + uuid.uuid4().hex
            destination = self.config.destination(request.destination) if request.destination else self.config.cache / f"{identifier}.{request.format}"
            if destination.exists() and not request.overwrite:
                raise FileExistsError("Destination exists. Choose another name or explicitly enable overwrite.")
            fd, temp_name = tempfile.mkstemp(prefix=".peptacular-export-", dir=destination.parent)
            temp = Path(temp_name)
            try:
                with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
                    if request.format == "json":
                        stream.write(encode({"contract_version": "1.0", "metadata": obj["metadata"], "records": rows}))
                    elif request.format == "jsonl":
                        for row in rows:
                            stream.write(encode(row) + "\n")
                    elif request.format == "csv":
                        columns = sorted(set().union(*(row.keys() for row in rows))) if rows else []
                        writer = csv.DictWriter(stream, fieldnames=columns)
                        writer.writeheader()
                        for row in rows:
                            writer.writerow({k: csv_value(v) for k, v in row.items()})
                    else:
                        for index, row in enumerate(rows):
                            sequence = row.get(request.sequence_column)
                            if not isinstance(sequence, str) or not sequence or not sequence.isalpha() or not sequence.isascii():
                                raise ValueError("FASTA export requires a plain residue sequence column on every row")
                            name = str(row.get("source_id") or row.get("id") or index).replace("\n", " ").replace("\r", " ")
                            stream.write(f">{name}\n{sequence}\n")
                    stream.flush()
                    os.fsync(stream.fileno())
                if temp.stat().st_size > self.config.storage_bytes:
                    raise ServiceError("storage_limit", "Export exceeds the byte quota.")
                result = {
                    "export_id": identifier,
                    "path": str(destination),
                    "format": request.format,
                    "record_count": len(rows),
                    "resource_uri": f"peptacular://exports/{identifier}",
                    "sha256": hashlib.sha256(temp.read_bytes()).hexdigest(),
                }
                if request.format == "csv":
                    result["formula_like_text"] = "escaped_with_leading_apostrophe"
                # Reserve metadata before publication. Files outside cache are always user-owned.
                saved_id = self.save(
                    "export",
                    destination.name,
                    result,
                    {
                        "parent_id": obj["id"],
                        "managed": request.destination is None,
                        "managed_file_bytes": temp.stat().st_size if request.destination is None else 0,
                    },
                    object_id=identifier,
                    idem=key,
                    request_hash=request_hash,
                )
                if saved_id != identifier:
                    return self.get(saved_id)["data"]
                try:
                    if request.destination and self.config.destination(request.destination) != destination:
                        raise ValueError("Destination changed during export")
                    if request.overwrite:
                        os.replace(temp, destination)
                    else:
                        os.link(temp, destination)
                except BaseException:
                    with self.connect() as db:
                        db.execute("DELETE FROM objects WHERE id = ?", (identifier,))
                    raise
            finally:
                temp.unlink(missing_ok=True)
            return result


def matches(value, clause):
    expected, op = clause.value, clause.operator
    if op == "is_null":
        if not isinstance(expected, bool):
            raise ValueError("is_null requires a boolean value")
        return (value is None) == expected
    if op == "in":
        if not isinstance(expected, list):
            raise ValueError("in requires a list value")
        return any(value == item and (isinstance(value, bool) == isinstance(item, bool)) for item in expected)
    if isinstance(expected, list):
        raise ValueError("Only in accepts a list")
    if op in ("eq", "ne"):
        equal = value == expected and (isinstance(value, bool) == isinstance(expected, bool))
        return equal if op == "eq" else not equal
    if value is None:
        return False
    if isinstance(value, bool) or isinstance(expected, bool) or not isinstance(value, (float, int)) or not isinstance(expected, (float, int)):
        raise ValueError("Range comparisons require numeric values")
    return {"lt": value < expected, "le": value <= expected, "gt": value > expected, "ge": value >= expected}[op]


def page_rows(rows, offset, limit, byte_limit, query_hash):
    page, size = [], 0
    for row in rows[offset : offset + limit]:
        row_size = len(encode(row).encode())
        if size + row_size > byte_limit:
            if not page:
                raise ServiceError("page_limit", "One row exceeds the page byte budget. Project fewer columns or export the result.")
            break
        size += row_size
        page.append(row)
    next_offset = offset + len(page)
    cursor = base64.urlsafe_b64encode(encode({"query": query_hash, "offset": next_offset}).encode()).decode() if next_offset < len(rows) else None
    return page, cursor

"""Spawned workers with enforced deadlines and instance-owned job lifetimes."""

import asyncio
import contextlib
import json
import logging
import multiprocessing
import os
import sys
import time
import uuid
from typing import Any, Literal, cast

from .contracts import Computation, Envelope, Page, Query
from .operations import ServiceError, diagnostic, run_operation, versions
from .storage import encode, fingerprint

_LOG = logging.getLogger(__name__)


def worker_main(name, payload, records, proteins, output):
    # A dependency print or warning must never become protocol traffic.
    sys.stdout = sys.stderr
    last_report = 0.0

    def report(progress):
        nonlocal last_report
        if time.monotonic() - last_report < 0.2:
            return
        with open(output + ".progress.tmp", "w", encoding="utf-8") as stream:
            stream.write(encode(progress))
        os.replace(output + ".progress.tmp", output + ".progress")
        last_report = time.monotonic()

    try:
        result = run_operation(name, payload, records, proteins, progress=report)
    except BaseException:
        _LOG.error("Scientific worker failed for operation %s", name)
        result = {"worker_error": "Unexpected calculation failure. Use the request ID to locate server logs."}
    with open(output, "w", encoding="utf-8") as stream:
        stream.write(encode(result))


class Runner:
    def __init__(self, store):
        self.store = store
        self.config = store.config
        self.owner = uuid.uuid4().hex
        self.semaphore = asyncio.Semaphore(self.config.workers)
        self.tasks = {}
        self.closed = False
        self.active = {}
        self.progress_paths = {}
        self.idempotency_lock = asyncio.Lock()

    async def submit(self, name, request, request_id):
        if self.closed:
            raise ServiceError("server_stopping", "Server is shutting down.")
        records = await asyncio.to_thread(self.store.resolve_inputs, request.inputs)
        proteins = await asyncio.to_thread(self.store.resolve_inputs, request.proteins) if name == "map_peptides" else None
        characters = sum(len(encode(row["annotation"])) for row in records)
        large = len(records) > 25 or characters > 5000 or name in ("enumerate_modifications", "isotope_envelopes")
        if proteins:
            large = large or len(proteins) * len(records) > 100
            if len(proteins) * len(records) > 100000:
                raise ServiceError("resource_limit", "Mapping exceeds 100,000 peptide/protein pairs. Split the inputs.")
        mode = "job" if request.execution.mode == "job" or (request.execution.mode == "auto" and large) else "inline"
        settings = request.model_dump(exclude={"inputs", "execution", "proteins", "reference"})
        settings["versions"] = versions()
        if request.preflight:
            return Envelope(
                request_id=request_id,
                applied_settings=settings,
                records=[
                    {
                        "mode": mode,
                        "input_records": len(records),
                        "input_characters": characters,
                        "max_rows": request.max_rows,
                        "worker_limit": self.config.workers,
                    }
                ],
                page=Page(returned_rows=1, total_rows=1),
            )
        if request.execution.mode == "inline" and (len(records) > 100 or characters > 1000000):
            raise ServiceError("inline_limit", "Inline input exceeds 100 records or 1 MiB. Select auto or job execution.")
        payload = request.model_dump()
        request_hash = fingerprint({"tool": name, "request": payload, "versions": versions()})
        key = f"scientific:{request.execution.idempotency_key}" if request.execution.idempotency_key else None
        async with self.idempotency_lock:
            prior = await asyncio.to_thread(self.store.idempotent, key, request_hash)
            if prior:
                if prior["kind"] == "job":
                    return self.job_envelope(prior, request_id)
                return await self.result_envelope(prior["id"], request_id)
            if len(self.tasks) >= 32:
                raise ServiceError("queue_limit", "Server already has 32 pending requests. Wait for a job to complete.")
            if mode == "job":
                job_id = await asyncio.to_thread(
                    self.store.save,
                    "job",
                    name,
                    {"state": "queued", "created": time.time(), "request_id": request_id},
                    {"owner": self.owner, "operation": name},
                    idem=key,
                    request_hash=request_hash,
                )
                task = asyncio.create_task(self._job(job_id, name, payload, records, proteins, settings, request_id))
                self.tasks[job_id] = task
                task.add_done_callback(lambda _: self.tasks.pop(job_id, None))
                return Envelope(
                    request_id=request_id,
                    status="queued",
                    job_id=job_id,
                    applied_settings=settings,
                    computation=Computation(complete=False, stop_reason="queued"),
                    page=Page(total_rows=None),
                )
        # Inline requests use the same isolated worker path, with a shorter deadline.
        task = asyncio.current_task()
        self.tasks[request_id] = task
        try:
            result = await self._compute(request_id, name, payload, records, proteins, min(request.execution.timeout_seconds, self.config.inline_seconds))
            result_id = await asyncio.to_thread(self._save_result, name, result, settings, idem=key, request_hash=request_hash)
            return await self.result_envelope(result_id, request_id)
        finally:
            self.tasks.pop(request_id, None)

    async def _compute(self, identifier, name, payload, records, proteins, timeout):
        # Include queue time in the deadline. No work runs on the protocol loop.
        async with asyncio.timeout(timeout):
            async with self.semaphore:
                output = self.config.cache / f"worker_{self.owner}_{uuid.uuid4().hex}.json"
                if identifier.startswith("job_"):
                    job = self.store.get(identifier)
                    self.store.update_job(identifier, self.owner, {**job["data"], "state": "running", "started": time.time()})
                process = multiprocessing.get_context("spawn").Process(
                    target=worker_main,
                    args=(name, payload, records, proteins, str(output)),
                    daemon=True,
                )
                process.start()
                self.active[identifier] = process
                self.progress_paths[identifier] = output.with_suffix(output.suffix + ".progress")
                try:
                    while process.is_alive():
                        await asyncio.sleep(0.05)
                    process.join()
                    if process.exitcode != 0 or not output.exists():
                        raise ServiceError("worker_failure", "Calculation worker exited without a result.")
                    if output.stat().st_size > 34000000:
                        raise ServiceError("byte_limit", "Worker output exceeded the byte budget.")
                    result = await asyncio.to_thread(lambda: json.loads(output.read_text(encoding="utf-8")))
                    if "worker_error" in result:
                        raise ServiceError("internal_error", result["worker_error"])
                    return result
                finally:
                    if process.is_alive():
                        process.terminate()
                        await asyncio.to_thread(process.join, 1)
                    if process.is_alive():
                        process.kill()
                        await asyncio.to_thread(process.join, 1)
                    process.close()
                    self.active.pop(identifier, None)
                    self.progress_paths.pop(identifier, None)
                    output.unlink(missing_ok=True)
                    output.with_suffix(output.suffix + ".progress").unlink(missing_ok=True)
                    output.with_suffix(output.suffix + ".progress.tmp").unlink(missing_ok=True)

    def _save_result(self, name, result, settings, **kwargs):
        metadata = {
            "operation": name,
            "settings": settings,
            "computation": result["computation"],
            "diagnostics": result["diagnostics"],
            "progress": result["progress"],
        }
        return self.store.save("result", name, result["records"], metadata, **kwargs)

    async def _job(self, job_id, name, payload, records, proteins, settings, request_id):
        created = time.time()
        data = {"state": "queued", "created": created, "request_id": request_id, "progress": {"records_consumed": 0, "rows_generated": 0}}
        self.store.update_job(job_id, self.owner, data)
        try:
            result = await self._compute(job_id, name, payload, records, proteins, payload["execution"]["timeout_seconds"])
            identifier = await asyncio.to_thread(self._save_result, name, result, settings)
            partial = not result["computation"]["complete"] or any(row.get("diagnostics") for row in result["records"])
            data.update(
                state="partially_succeeded" if partial else "succeeded", result_id=identifier, progress=result["progress"], computation=result["computation"]
            )
        except asyncio.CancelledError:
            data.update(state="interrupted" if self.closed else "cancelled")
        except TimeoutError:
            data.update(state="failed", diagnostics=[diagnostic(ServiceError("timeout", "The job deadline expired and its worker was stopped."))])
        except Exception as exc:
            if isinstance(exc, (ValueError, OSError)):
                data.update(state="failed", diagnostics=[diagnostic(exc)])
            else:
                _LOG.error("Job %s failed with %s", job_id, type(exc).__name__)
                data.update(state="failed", diagnostics=[diagnostic(ServiceError("internal_error", f"Unexpected failure for request {request_id}."))])
        finally:
            data["elapsed_seconds"] = time.time() - created
            self.store.update_job(job_id, self.owner, data)

    async def result_envelope(self, identifier, request_id):
        obj = await asyncio.to_thread(self.store.get, identifier)
        metadata = obj["metadata"]
        result: dict[str, Any]
        try:
            result = await asyncio.to_thread(self.store.query, Query(result_id=identifier), byte_limit=self.config.preview_bytes)
        except ServiceError as exc:
            if exc.code != "page_limit":
                raise
            result = {"records": [], "page": {"returned_rows": 0, "total_rows": len(obj["data"]), "next_cursor": None}}
        partial = not metadata["computation"]["complete"] or any(row.get("diagnostics") for row in obj["data"])
        return Envelope(
            request_id=request_id,
            result_id=identifier,
            status="partial" if partial else "complete",
            records=result["records"],
            page=Page.model_validate(result["page"]),
            applied_settings=metadata["settings"],
            computation=metadata["computation"],
            diagnostics=metadata.get("diagnostics", [])[:25],
        )

    def job_envelope(self, obj, request_id):
        data = dict(obj["data"])
        state = data["state"]
        if obj["metadata"].get("owner") != self.owner and state in ("queued", "running"):
            data.update(state="unavailable", explanation="The job belongs to another server instance. Its live state cannot be established here.")
            state = "interrupted"
        status: Any = {"succeeded": "complete", "partially_succeeded": "partial", "failed": "error"}.get(state, state)
        return Envelope(
            request_id=request_id,
            job_id=obj["id"],
            result_id=data.get("result_id"),
            status=cast(Literal["complete", "partial", "error", "queued", "running", "cancelled", "interrupted"], status),
            records=[data],
            page=Page(returned_rows=1, total_rows=1),
            computation=Computation.model_validate(
                data.get(
                    "computation",
                    {
                        "complete": state == "succeeded",
                        "stop_reason": None if state == "succeeded" else state,
                    },
                )
            ),
        )

    async def get_job(self, identifier, request_id):
        obj = await asyncio.to_thread(self.store.get, identifier)
        if obj["kind"] != "job":
            raise ServiceError("invalid_reference_type", "Expected a job ID.")
        if identifier in self.progress_paths:
            with contextlib.suppress(FileNotFoundError):
                progress = await asyncio.to_thread(self.progress_paths[identifier].read_text, encoding="utf-8")
                obj["data"]["progress"] = json.loads(progress)
        if "created" in obj["data"] and obj["data"]["state"] in ("queued", "running"):
            obj["data"]["elapsed_seconds"] = time.time() - obj["data"]["created"]
        return self.job_envelope(obj, request_id)

    async def cancel(self, identifier, request_id):
        obj = await asyncio.to_thread(self.store.get, identifier)
        if obj["kind"] != "job" or obj["metadata"].get("owner") != self.owner:
            raise ServiceError("job_owner", "Only this server's own jobs can be cancelled.")
        task = self.tasks.get(identifier)
        if task:
            task.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await task
            # A queued coroutine may be cancelled before its first instruction.
            latest = self.store.get(identifier)
            if latest["data"]["state"] == "queued":
                self.store.update_job(identifier, self.owner, {**latest["data"], "state": "cancelled"})
        return await self.get_job(identifier, request_id)

    async def close(self):
        self.closed = True
        identifiers = list(self.tasks)
        tasks = list(self.tasks.values())
        for task in tasks:
            task.cancel()
        await asyncio.gather(*tasks, return_exceptions=True)
        for identifier in identifiers:
            if identifier.startswith("job_"):
                job = self.store.get(identifier)
                if job["data"]["state"] in ("queued", "running"):
                    self.store.update_job(identifier, self.owner, {**job["data"], "state": "interrupted"})
        self.tasks.clear()

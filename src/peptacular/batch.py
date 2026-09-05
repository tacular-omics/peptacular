"""Bounded sequence batches with optional collection of input diagnostics."""

import inspect
import multiprocessing as mp
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from functools import partial
from itertools import islice
from typing import Any, Literal

from .annotation import ProFormaAnnotation
from .constants import parallelMethod, parallelMethodLiteral
from .diagnostics import Diagnostic, diagnostic_from_exception
from .sequence.parallel import AUTO_PARALLEL_MIN_ITEMS, _get_optimal_method, _validate_positive_int

__all__ = ["BatchResult", "BatchOperation", "iter_batch", "batch", "diagnose"]

BatchOperation = Literal["parse", "mass", "mz", "comp", "fragment", "fast_fragment", "digest", "isotopic_distribution"]
_OPERATIONS = frozenset({"parse", "mass", "mz", "comp", "fragment", "fast_fragment", "digest", "isotopic_distribution"})


@dataclass(frozen=True)
class BatchResult:
    """One ordered batch result. A successful value may itself be ``None``.

    :param index: Zero-based position in the original iterable.
    :param input: Original sequence string or annotation.
    :param value: Operation result, or ``None`` on input failure.
    :param error: Input diagnostic, or ``None`` on success.
    """

    index: int
    input: str | ProFormaAnnotation
    value: Any = None
    error: Diagnostic | None = None

    @property
    def ok(self) -> bool:
        """Whether the operation succeeded."""
        return self.error is None


def _validate_operation(operation: BatchOperation, kwargs: dict[str, Any]) -> None:
    if operation not in _OPERATIONS:
        raise ValueError(f"Unknown batch operation {operation!r}. Choose from {', '.join(sorted(_OPERATIONS))}.")
    # Bad keywords and missing required arguments are configuration errors,
    # not failures to repeat for every sequence in a database.
    inspect.signature(getattr(ProFormaAnnotation, operation)).bind(None, **kwargs)


def _run_item(
    item: tuple[int, str | ProFormaAnnotation],
    *,
    operation: BatchOperation,
    kwargs: dict[str, Any],
    errors: Literal["raise", "collect"],
) -> BatchResult:
    index, sequence = item
    stage: Literal["parse", "validate", "calculate"] = "parse"
    if not isinstance(sequence, (str, ProFormaAnnotation)):
        if errors == "raise":
            raise TypeError("Batch inputs must be sequence strings or ProFormaAnnotation objects")
        return BatchResult(index, sequence, error=Diagnostic("invalid_input", "parse", "Expected a sequence string or annotation", "TypeError"))
    try:
        annotation = ProFormaAnnotation.parse(sequence) if isinstance(sequence, str) else sequence.copy()
        stage = "validate"
        if operation != "parse":
            annotation.validate_sequence()
            annotation.validate_ambiguous_labels()
        if operation == "parse":
            if kwargs.get("validate", False):
                annotation.validate_annotation()
            value = annotation
        else:
            stage = "calculate"
            value = getattr(annotation, operation)(**kwargs)
            if operation == "digest":
                # Consume lazy failures here and return a process-safe value.
                value = list(value)
        return BatchResult(index, sequence, value=value)
    except (ValueError, KeyError) as exc:
        if errors == "raise":
            raise
        return BatchResult(index, sequence, error=diagnostic_from_exception(exc, stage))


def diagnose(sequence: str | ProFormaAnnotation, operation: BatchOperation = "mass", **kwargs: Any) -> Diagnostic | None:
    """Run one operation and return its input diagnostic, or ``None`` on success.

    Parsing, annotation validation, and calculation failures are distinguished.
    This executes the requested operation, so fragment or isotope diagnostics
    have the same computational cost as the corresponding calculation.
    Invalid operation names, bad keyword arguments, and unexpected internal or
    infrastructure errors still raise. The input annotation is never mutated.
    """
    _validate_operation(operation, kwargs)
    return _run_item((0, sequence), operation=operation, kwargs=kwargs, errors="collect").error


def iter_batch(
    operation: BatchOperation,
    sequences: Iterable[str | ProFormaAnnotation],
    *,
    errors: Literal["raise", "collect"] = "raise",
    batch_size: int = 1000,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    start_method: Literal["spawn", "fork", "forkserver"] | None = None,
    **kwargs: Any,
) -> Iterator[BatchResult]:
    """Yield indexed results while holding at most ``batch_size`` inputs at once.

    ``errors='raise'`` propagates input errors. ``errors='collect'`` records
    expected input failures and continues. Configuration, iterator, unexpected
    programming, and worker failures propagate in either mode. Each input is
    checked for valid residues and ambiguous labels before calculation. Parse
    operations retain the usual ``validate=False`` default. Results retain input order.

    Automatic small batches run sequentially. Explicit process/thread settings
    override that heuristic. A pool is reused across bounded batches and closed
    on exhaustion, error, or iterator ``close()``. Closing waits for already
    running work. Use ``contextlib.closing`` for early termination.
    Process execution needs a guarded script entry point on spawn platforms.
    ``start_method`` selects a local process context without changing global
    multiprocessing state. Each result may itself contain a large value, such
    as a full fragment list, so ``batch_size`` also controls retained results.

    :param operation: Annotation operation, such as ``mass``, ``comp``, or ``fragment``.
    :param sequences: An iterable of sequence strings or annotations.
    :param errors: Raise on an input failure or collect a diagnostic.
    :param batch_size: Maximum number of inputs collected for one execution batch.
    :param n_workers: Positive worker limit. Defaults to available CPU count.
    :param chunksize: Positive process-map chunk size. Defaults to 1.
    :param method: Sequential, thread, process, or automatic execution.
    :param start_method: Optional multiprocessing context for process execution.
    :param kwargs: Keyword arguments for the annotation operation.
    :return: Ordered results with zero-based input indexes.
    """
    _validate_operation(operation, kwargs)
    if errors not in ("raise", "collect"):
        raise ValueError("errors must be 'raise' or 'collect'")
    if isinstance(sequences, (str, ProFormaAnnotation)):
        raise TypeError("sequences must be an iterable of inputs, not a single sequence")
    _validate_positive_int(batch_size, "batch_size")
    if batch_size is None:
        raise ValueError("batch_size must be a positive integer")
    _validate_positive_int(n_workers, "n_workers")
    _validate_positive_int(chunksize, "chunksize")
    selected = parallelMethod(method) if method is not None else parallelMethod(_get_optimal_method())
    context = mp.get_context(start_method) if start_method is not None else None
    if context is not None and selected != parallelMethod.PROCESS:
        raise ValueError("start_method requires process execution")
    source = enumerate(sequences)
    execute = partial(_run_item, operation=operation, kwargs=kwargs, errors=errors)
    executor: ProcessPoolExecutor | ThreadPoolExecutor | None = None
    try:
        while items := list(islice(source, batch_size)):
            small_auto = method is None and n_workers is None and start_method is None and len(items) < AUTO_PARALLEL_MIN_ITEMS
            if selected == parallelMethod.SEQUENTIAL or (small_auto and executor is None):
                yield from map(execute, items)
                continue
            if executor is None:
                workers = min(n_workers or _available_cpus(), len(items))
                if selected == parallelMethod.THREAD:
                    executor = ThreadPoolExecutor(max_workers=workers)
                else:
                    executor = ProcessPoolExecutor(max_workers=workers, mp_context=context)
            yield from executor.map(execute, items, chunksize=chunksize or 1)
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)


def _available_cpus() -> int:
    import os

    return getattr(os, "process_cpu_count", os.cpu_count)() or 1


def batch(
    operation: BatchOperation,
    sequences: Iterable[str | ProFormaAnnotation],
    *,
    errors: Literal["raise", "collect"] = "raise",
    batch_size: int = 1000,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    start_method: Literal["spawn", "fork", "forkserver"] | None = None,
    **kwargs: Any,
) -> list[BatchResult]:
    """Collect :func:`iter_batch` into a list. Defaults to raising on input errors."""
    return list(
        iter_batch(
            operation,
            sequences,
            errors=errors,
            batch_size=batch_size,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            start_method=start_method,
            **kwargs,
        )
    )

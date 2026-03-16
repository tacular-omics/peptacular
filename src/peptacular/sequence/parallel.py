from __future__ import annotations

import logging
import multiprocessing as mp
import sys
from collections.abc import Callable, Sequence
from functools import partial
from multiprocessing.pool import Pool, ThreadPool
from typing import Any, Literal, TypeVar

logger = logging.getLogger(__name__)

from ..constants import parallelMethod, parallelMethodLiteral

T = TypeVar("T")

# Global pool cache removed


def set_start_method(method: Literal["fork", "spawn", "forkserver"] | None = None):
    """
    Set the multiprocessing start method.

    :param method: 'fork', 'spawn', or 'forkserver'. If None, uses platform default.

    Notes:
        - fork: Fast but unsafe with threads/some libraries (Unix only)
        - spawn: Slow but safest (all platforms)
        - forkserver: Balanced (Unix only)

    Must be called before creating any pools.
    """
    if method is not None:
        try:
            mp.set_start_method(method, force=True)
        except RuntimeError as e:
            logger.warning(f"Warning: Could not set start method to '{method}': {e}")


def get_start_method() -> str:
    """
    Get the current multiprocessing start method.

    :return: Current start method ('fork', 'spawn', or 'forkserver')
    """
    return mp.get_start_method()


def get_available_start_methods() -> list[str]:
    """
    Get all available start methods on this platform.

    :return: List of available methods
    """
    return mp.get_all_start_methods()


def _is_gil_disabled() -> bool:
    """
    Detect if the GIL is disabled (Python 3.13+ with free-threading).
    :return: True if GIL is disabled, False otherwise
    """
    # Python 3.13+ with free-threading
    if sys.version_info >= (3, 13):
        try:
            if hasattr(sys, "_is_gil_enabled"):
                result = not sys._is_gil_enabled()
                return result
        except AttributeError:
            pass
    return False


def _get_optimal_method() -> Literal["process", "thread"]:
    """
    Automatically determine the best parallelization method.
    :return: 'thread' if GIL is disabled, 'process' otherwise
    """
    return "thread" if _is_gil_disabled() else "process"


def _apply_wrapper[T](item: Any, func: Callable[..., T], func_kwargs: dict[str, Any]) -> T:
    """
    Wrapper function at module level so it can be pickled.
    :param item: Single item to process
    :param func: Function to apply
    :param func_kwargs: Keyword arguments for the function
    :return: Result of func(item, **func_kwargs)
    """
    return func(item, **func_kwargs)


def _create_pool(method: parallelMethod, n_workers: int) -> Pool | ThreadPool:
    """
    Create a new pool.

    :param method: Process or thread pool
    :param n_workers: Number of workers
    :return: Pool instance
    """
    if method == parallelMethod.THREAD:
        return ThreadPool(processes=n_workers)
    else:
        return Pool(processes=n_workers)


def parallel_apply_internal[T](
    func: Callable[..., T],
    items: Sequence[Any],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = False,
    verbose: bool = False,
    **func_kwargs: Any,
) -> list[T]:
    """
    Internal function for parallel processing.
    Automatically detects if GIL is disabled and uses threading for better performance.
    Otherwise defaults to multiprocessing.
    Results are returned in the same order as the input items.

    :param func: Function to apply to each item
    :param items: Sequence of items to process
    :param n_workers: Number of worker processes/threads. If None, uses CPU count
    :param chunksize: Number of items per chunk. If None, auto-calculated
    :param method: 'process', 'thread', 'sequential', or None (auto-detect). Default is None.
    :param reuse_pool: Ignored. Kept for backward compatibility.
    :param verbose: If True, print worker and chunksize info
    :param func_kwargs: Keyword arguments to pass to the function
    :return: List of results in the same order as input items
    """
    # Convert to list if needed
    items_list = list(items)

    # Handle empty input
    if not items_list:
        return []

    method_enum = parallelMethod(method) if method is not None else parallelMethod.PROCESS

    # Handle sequential execution
    if method_enum == parallelMethod.SEQUENTIAL:
        return [func(item, **func_kwargs) for item in items_list]

    # Handle parallel execution
    if n_workers is None:
        n_workers = mp.cpu_count()

    if verbose:
        logger.debug(f"Using n_workers: {n_workers}")

    # Optimize chunksize calculation
    if chunksize is None:
        # Better heuristic: balance overhead vs parallelism
        total_items = len(items_list)
        if total_items < n_workers:
            chunksize = 1
        else:
            # Aim for ~4 chunks per worker to balance load
            chunksize = max(1, total_items // (n_workers * 4))

    if verbose:
        logger.debug(f"Using chunksize: {chunksize}")

    wrapper = partial(_apply_wrapper, func=func, func_kwargs=func_kwargs)

    # Create temporary pool
    pool = _create_pool(method_enum, n_workers)
    with pool:
        return pool.map(wrapper, items_list, chunksize=chunksize)

#!/usr/bin/env python3
"""Benchmark parallelization strategies for peptacular."""

import csv
import multiprocessing as mp
import random
import statistics
import sys
import time
from pathlib import Path
from typing import Literal

# import matplotlib.pyplot as plt
import peptacular as pt


def generate_test_data(n_sequences: int = 10_000):
    """Generate random ProForma sequences."""
    random.seed(42)
    sequences = [pt.ProFormaAnnotation.random(6, random.randint(10, 30)).serialize() for _ in range(n_sequences)]
    annotations = [pt.parse(seq) for seq in sequences]
    return sequences, annotations


def benchmark_method(
    data,
    n_workers: int,
    method: Literal["sequential", "thread", "process"],
    n_runs: int = 3,
) -> tuple[float, float]:
    """Run benchmark multiple times and return mean and stdev."""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        _ = pt.mass(data, n_workers=n_workers, method=method)
        times.append(time.perf_counter() - start)

    return statistics.mean(times), statistics.stdev(times) if len(times) > 1 else 0.0


def check_gil_status():
    """Check if GIL is disabled."""
    try:
        return not sys._is_gil_enabled()
    except AttributeError:
        return False


def save_results_to_tsv(results: dict, output_path: Path):
    """Save benchmark results to TSV file."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["data_type", "workers", "method", "mean_time", "std_time", "speedup"])

        for data_type, data in results.items():
            baseline = data["baseline"]
            for result in data["results"]:
                workers = result["workers"]
                for method in ["sequential", "thread", "process"]:
                    mean_time, std_time, speedup = result[method]
                    writer.writerow([data_type, workers, method, f"{mean_time:.4f}", f"{std_time:.4f}", f"{speedup:.2f}"])


def main():
    pt.set_start_method("fork")

    # Config
    NUM_SEQUENCES = 10_000
    NUM_RUNS = 2
    MAX_WORKERS = mp.cpu_count()
    WORKER_COUNTS = list(range(1, MAX_WORKERS * 2))

    OUTPUT_DIR = Path(__file__).parent / "benchmark_results"
    OUTPUT_DIR.mkdir(exist_ok=True)
    timestamp = time.strftime("%Y%m%d_%H%M%S")

    # System info
    print("PEPTACULAR PARALLELIZATION BENCHMARK")
    print(f"Python: {sys.version.split()[0]}")
    print(f"GIL: {'disabled' if check_gil_status() else 'enabled'}")
    print(f"CPUs: {MAX_WORKERS}")
    print(f"Sequences: {NUM_SEQUENCES:,}")
    print(f"Runs per test: {NUM_RUNS}")
    print()

    # Generate test data
    sequences, annotations = generate_test_data(NUM_SEQUENCES)

    # Warmup
    _ = pt.mass(sequences, n_workers=1, method="sequential")
    _ = pt.mass(sequences, n_workers=1, method="process")
    _ = pt.mass(sequences, n_workers=1, method="thread")

    # Baseline
    baseline_seq, _ = benchmark_method(sequences, 1, "sequential", NUM_RUNS)
    baseline_ann, _ = benchmark_method(annotations, 1, "sequential", NUM_RUNS)

    print("Baselines (1 worker, sequential):")
    print(f"  Strings: {baseline_seq:.3f}s")
    print(f"  Annotations: {baseline_ann:.3f}s ({baseline_ann / baseline_seq:.2f}x overhead)")
    print()

    # Store results
    all_results = {}

    # Test both data types
    for data, data_type, baseline in [
        (sequences, "strings", baseline_seq),
        (annotations, "annotations", baseline_ann),
    ]:
        print(f"{data_type.capitalize()}:")
        print(f"{'Workers':<10} {'Sequential':<15} {'Thread':<15} {'Process':<15}")

        results = []
        for n in WORKER_COUNTS:
            seq_time, seq_std = benchmark_method(data, n, "sequential", NUM_RUNS)
            thread_time, thread_std = benchmark_method(data, n, "thread", NUM_RUNS)
            process_time, process_std = benchmark_method(data, n, "process", NUM_RUNS)

            seq_speedup = baseline / seq_time
            thread_speedup = baseline / thread_time
            process_speedup = baseline / process_time

            results.append(
                {
                    "workers": n,
                    "sequential": (seq_time, seq_std, seq_speedup),
                    "thread": (thread_time, thread_std, thread_speedup),
                    "process": (process_time, process_std, process_speedup),
                }
            )

            print(f"{n:<10} {seq_time:.3f}s ({seq_speedup:.1f}x)  {thread_time:.3f}s ({thread_speedup:.1f}x)  {process_time:.3f}s ({process_speedup:.1f}x)")

        all_results[data_type] = {"baseline": baseline, "results": results}
        print()

    # Save results
    tsv_path = OUTPUT_DIR / f"scaling_benchmark_{timestamp}.tsv"
    save_results_to_tsv(all_results, tsv_path)
    print(f"Results saved to: {tsv_path}")


if __name__ == "__main__":
    main()

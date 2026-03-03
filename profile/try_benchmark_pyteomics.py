"""
Mass Calculation Benchmark: Peptacular vs Pyteomics
Comprehensive testing with multiple workload patterns
"""

import random
import statistics
import time

from pyteomics import mass

import peptacular as pt

# Configure parallel processing
pt.set_start_method("fork")


def generate_random_proforma(length: int = 20, mod_probability: float = 0.2) -> str:
    """Generate a random ProForma sequence with modifications."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    mods = ["[-15.995]", "[+57.021]", "[Oxidation]", "[+42.011]", "[unimod:1]"]

    sequence_parts: list[str] = []

    if random.random() < mod_probability:
        sequence_parts.append(random.choice(mods))
        sequence_parts.append("-")

    for _ in range(length):
        aa = random.choice(amino_acids)
        sequence_parts.append(aa)
        if random.random() < mod_probability:
            sequence_parts.append(random.choice(mods))

    if random.random() < mod_probability:
        sequence_parts.append("-")
        sequence_parts.append(random.choice(mods))

    return "".join(sequence_parts)


def time_function(func, iterations: int, warmup: int = 3) -> tuple[list[float], list]:
    """Time a function with warmup runs."""
    # Warmup
    for _ in range(warmup):
        result = func()

    # Actual timing
    times = []
    results = []
    for _ in range(iterations):
        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
        times.append(elapsed * 1_000_000)  # Convert to microseconds
        results.append(result)

    return times, results


def benchmark_single_sequence_workload():
    """Benchmark: Individual sequence processing (typical API usage)"""

    num_sequences = 1000
    iterations_per_seq = 10

    print(f"\n{'=' * 80}")
    print("BENCHMARK 1: Single Sequence Workload (per-sequence timing)")
    print(f"{'=' * 80}")
    print(f"Sequences: {num_sequences}")
    print(f"Iterations per sequence: {iterations_per_seq}")
    print(f"Total calculations: {num_sequences * iterations_per_seq:,}\n")

    # Generate test data
    sequences = [generate_random_proforma(length=random.randint(10, 30)) for _ in range(num_sequences)]

    print("Example sequences:")
    for i, seq in enumerate(sequences[:3]):
        print(f"  {seq}")

    # Benchmark Peptacular (sequential)
    print("\nTiming Peptacular (sequential)...")
    pept_times = []
    pept_masses = []

    for seq in sequences:
        seq_times, seq_masses = time_function(lambda s=seq: pt.mass(s), iterations=iterations_per_seq, warmup=1)
        pept_times.extend(seq_times)
        pept_masses.append(seq_masses[0])

    # Benchmark Pyteomics
    print("Timing Pyteomics...")
    pyto_times = []
    pyto_masses = []

    for seq in sequences:
        seq_times, seq_masses = time_function(lambda s=seq: mass.calculate_mass(proforma=s), iterations=iterations_per_seq, warmup=1)
        pyto_times.extend(seq_times)
        pyto_masses.append(seq_masses[0])

    # Statistics
    pept_mean = statistics.mean(pept_times)
    pept_median = statistics.median(pept_times)
    pept_stdev = statistics.stdev(pept_times) if len(pept_times) > 1 else 0

    pyto_mean = statistics.mean(pyto_times)
    pyto_median = statistics.median(pyto_times)
    pyto_stdev = statistics.stdev(pyto_times) if len(pyto_times) > 1 else 0

    speedup = pyto_mean / pept_mean

    # Mass agreement
    mass_diffs = [abs(m1 - m2) for m1, m2 in zip(pept_masses, pyto_masses)]
    max_diff = max(mass_diffs)
    avg_diff = statistics.mean(mass_diffs)

    # Results
    print(f"\n{'Method':<20} {'Mean (µs)':<12} {'Median (µs)':<12} {'StdDev (µs)':<12}")
    print(f"{'-' * 56}")
    print(f"{'Peptacular':<20} {pept_mean:>10.2f}   {pept_median:>10.2f}   {pept_stdev:>10.2f}")
    print(f"{'Pyteomics':<20} {pyto_mean:>10.2f}   {pyto_median:>10.2f}   {pyto_stdev:>10.2f}")
    print(f"{'-' * 56}")

    if speedup > 1:
        print(f"✅ Peptacular is {speedup:.2f}x FASTER")
    else:
        print(f"⚠️  Pyteomics is {1 / speedup:.2f}x FASTER")

    print(f"\nMass Agreement: Max diff = {max_diff:.6f} Da, Avg diff = {avg_diff:.6f} Da")

    return pept_times, pyto_times


def benchmark_batch_workload():
    """Benchmark: Batch processing (parallel vs sequential)"""

    batch_sizes = [100, 500, 1000, 5000]
    iterations = 10

    print(f"\n{'=' * 80}")
    print("BENCHMARK 2: Batch Workload (total batch time)")
    print(f"{'=' * 80}")
    print(f"Iterations per batch size: {iterations}")
    print(f"Testing batch sizes: {batch_sizes}\n")

    for batch_size in batch_sizes:
        print(f"\n--- Batch Size: {batch_size} sequences ---")

        # Generate batch
        sequences = [generate_random_proforma(length=random.randint(10, 30)) for _ in range(batch_size)]

        # Peptacular Sequential
        print("  Peptacular (sequential)...", end=" ", flush=True)
        pept_seq_times, _ = time_function(lambda: [pt.mass(s) for s in sequences], iterations=iterations, warmup=2)
        pept_seq_mean = statistics.mean(pept_seq_times)
        print(f"{pept_seq_mean / 1000:.2f} ms")

        # Peptacular Parallel
        print("  Peptacular (parallel)...  ", end=" ", flush=True)
        # Extra warmup to spawn processes
        for _ in range(3):
            _ = pt.mass(sequences)

        pept_par_times, pept_masses = time_function(
            lambda: pt.mass(sequences),
            iterations=iterations,
            warmup=0,  # Already warmed up
        )
        pept_par_mean = statistics.mean(pept_par_times)
        print(f"{pept_par_mean / 1000:.2f} ms")

        # Pyteomics (sequential only)
        print("  Pyteomics (sequential)... ", end=" ", flush=True)
        pyto_times, pyto_masses = time_function(lambda: [mass.calculate_mass(proforma=s) for s in sequences], iterations=iterations, warmup=2)
        pyto_mean = statistics.mean(pyto_times)
        print(f"{pyto_mean / 1000:.2f} ms")

        # Mass agreement
        mass_diffs = [abs(float(m1) - m2) for m1, m2 in zip(pept_masses[0], pyto_masses[0])]
        max_diff = max(mass_diffs)

        # Speedups
        speedup_seq = pyto_mean / pept_seq_mean
        speedup_par = pyto_mean / pept_par_mean
        parallel_efficiency = pept_seq_mean / pept_par_mean

        print(f"\n  Results:")
        print(f"    Peptacular seq vs Pyteomics:  {speedup_seq:>5.2f}x")
        print(f"    Peptacular par vs Pyteomics:  {speedup_par:>5.2f}x")
        print(f"    Parallel speedup (internal):  {parallel_efficiency:>5.2f}x")
        print(f"    Mass agreement: max Δ = {max_diff:.6f} Da")


def benchmark_per_sequence_in_batch():
    """Benchmark: Per-sequence cost within batches (amortized overhead)"""

    batch_size = 1000
    iterations = 20

    print(f"\n{'=' * 80}")
    print("BENCHMARK 3: Per-Sequence Cost in Batch Processing")
    print(f"{'=' * 80}")
    print(f"Batch size: {batch_size} sequences")
    print(f"Iterations: {iterations}\n")

    sequences = [generate_random_proforma(length=random.randint(10, 30)) for _ in range(batch_size)]

    # Peptacular parallel
    print("Warming up parallel processing...")
    for _ in range(5):
        _ = pt.mass(sequences)

    print("Timing parallel batch processing...")
    batch_times, _ = time_function(lambda: pt.mass(sequences), iterations=iterations, warmup=0)

    batch_mean = statistics.mean(batch_times)
    per_sequence_cost = batch_mean / batch_size
    throughput = (batch_size * 1_000_000) / batch_mean  # sequences per second

    print(f"\nBatch processing {batch_size} sequences:")
    print(f"  Total time:        {batch_mean / 1000:>8.2f} ms")
    print(f"  Per-sequence cost: {per_sequence_cost:>8.2f} µs")
    print(f"  Throughput:        {throughput:>8,.0f} sequences/sec")


def benchmark_sequence_length_scaling():
    """Benchmark: How performance scales with sequence length"""

    lengths = [10, 20, 50, 100, 200]
    iterations = 100

    print(f"\n{'=' * 80}")
    print("BENCHMARK 4: Scaling with Sequence Length")
    print(f"{'=' * 80}")
    print(f"Iterations per length: {iterations}\n")

    print(f"{'Length':<10} {'Peptacular (µs)':<18} {'Pyteomics (µs)':<18} {'Speedup':<10}")
    print(f"{'-' * 56}")

    for length in lengths:
        # Generate sequences of specific length
        sequences = [generate_random_proforma(length=length, mod_probability=0.2) for _ in range(iterations)]

        # Time both implementations
        pept_times, _ = time_function(lambda i=0: pt.mass(sequences[i % len(sequences)]), iterations=iterations, warmup=5)

        pyto_times, _ = time_function(lambda i=0: mass.calculate_mass(proforma=sequences[i % len(sequences)]), iterations=iterations, warmup=5)

        pept_mean = statistics.mean(pept_times)
        pyto_mean = statistics.mean(pyto_times)
        speedup = pyto_mean / pept_mean

        print(f"{length:<10} {pept_mean:>16.2f}   {pyto_mean:>16.2f}   {speedup:>8.2f}x")


if __name__ == "__main__":
    print("Mass Calculation Benchmarks")
    print(f"Parallel backend: {pt.get_start_method()}")

    # Run all benchmarks
    benchmark_single_sequence_workload()
    benchmark_batch_workload()
    benchmark_per_sequence_in_batch()
    benchmark_sequence_length_scaling()

    print(f"\n{'=' * 80}")
    print("All benchmarks complete!")
    print(f"{'=' * 80}\n")

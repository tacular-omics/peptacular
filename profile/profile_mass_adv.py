import cProfile
import pstats
import random
import time
from io import StringIO
from pathlib import Path
from typing import Any, Callable

import peptacular as pt


def generate_random_protein(length: int = 500) -> str:
    """Generate a random protein sequence."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(length))


class Profiler:
    """Robust profiling context manager with multiple output formats."""

    def __init__(self, name: str = "profile", output_dir: Path | None = None):
        self.name = name
        self.output_dir = output_dir or Path("profile_results")
        self.output_dir.mkdir(exist_ok=True)
        self.profiler = cProfile.Profile()
        self.start_time = None
        self.end_time = None

    def __enter__(self):
        self.start_time = time.perf_counter()
        self.profiler.enable()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.profiler.disable()
        self.end_time = time.perf_counter()
        return False

    @property
    def elapsed(self) -> float:
        """Get elapsed time in seconds."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return 0.0

    def save_stats(self):
        """Save profiling stats to file."""
        stats_file = self.output_dir / f"{self.name}.prof"
        self.profiler.dump_stats(str(stats_file))
        print(f"\nSaved raw stats to: {stats_file}")

    def print_stats(self, sort_by: str = "cumulative", top_n: int = 30):
        """Print profiling statistics."""
        s = StringIO()
        ps = pstats.Stats(self.profiler, stream=s)
        ps.strip_dirs()
        ps.sort_stats(sort_by)
        ps.print_stats(top_n)

        output = s.getvalue()
        print(output)

        # Save to file
        txt_file = self.output_dir / f"{self.name}_{sort_by}.txt"
        txt_file.write_text(output)
        print(f"Saved {sort_by} stats to: {txt_file}")

        return output

    def print_callers(self, top_n: int = 20):
        """Print caller statistics."""
        s = StringIO()
        ps = pstats.Stats(self.profiler, stream=s)
        ps.strip_dirs()
        ps.sort_stats("cumulative")
        ps.print_callers(top_n)

        output = s.getvalue()
        print("\n" + "=" * 80)
        print("CALLER STATISTICS")
        print("=" * 80)
        print(output)

        txt_file = self.output_dir / f"{self.name}_callers.txt"
        txt_file.write_text(output)

        return output

    def summary(self, iterations: int = 1):
        """Print summary statistics."""
        print("\n" + "=" * 80)
        print("PROFILING SUMMARY")
        print("=" * 80)
        print(f"Total time: {self.elapsed:.3f}s")
        if iterations > 1:
            print(f"Iterations: {iterations}")
            print(f"Time per iteration: {self.elapsed / iterations * 1000:.3f}ms")
        print("=" * 80 + "\n")


def benchmark(func: Callable, iterations: int = 1, warmup: int = 0) -> dict[str, Any]:
    """
    Benchmark a function with optional warmup runs.

    Returns:
        dict with timing statistics
    """
    # Warmup
    for _ in range(warmup):
        func()

    # Actual timing
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    times.sort()
    n = len(times)

    return {
        "result": result,
        "iterations": iterations,
        "total": sum(times),
        "mean": sum(times) / n,
        "median": times[n // 2],
        "min": times[0],
        "max": times[-1],
        "std": (sum((t - sum(times) / n) ** 2 for t in times) / n) ** 0.5 if n > 1 else 0,
    }


def print_benchmark_results(name: str, stats: dict[str, Any], units: str = "ms"):
    """Pretty print benchmark results."""
    multiplier = 1000 if units == "ms" else 1

    print(f"\n{'=' * 60}")
    print(f"BENCHMARK: {name}")
    print(f"{'=' * 60}")
    print(f"Iterations: {stats['iterations']}")
    print(f"Total time: {stats['total']:.3f}s")
    print(f"Mean: {stats['mean'] * multiplier:.3f}{units}")
    print(f"Median: {stats['median'] * multiplier:.3f}{units}")
    print(f"Min: {stats['min'] * multiplier:.3f}{units}")
    print(f"Max: {stats['max'] * multiplier:.3f}{units}")
    print(f"Std Dev: {stats['std'] * multiplier:.3f}{units}")
    print(f"{'=' * 60}\n")


def run_mass_single(protein: str) -> float:
    """Run mass calculation on a single protein."""
    protein = "{Oxidation}M[+100]" + protein + "K-[Oxidation]"
    masses = pt.mass(protein, charge=0, ion_type="p")
    return masses


def profile_mass_single(num_proteins: int = 10, protein_length: int = 20, iterations_per_protein: int = 500, sort_by: str = "cumulative"):
    """Profile mass calculation for single protein peptides."""

    # Generate test data
    proteins = [generate_random_protein(protein_length) for _ in range(num_proteins)]
    total_calls = num_proteins * iterations_per_protein

    print(f"Testing with {num_proteins} proteins, {iterations_per_protein} iterations each")
    print(f"Total calls: {total_calls}")

    # Profile
    with Profiler(name="mass_single") as prof:
        for protein in proteins:
            for _ in range(iterations_per_protein):
                result = run_mass_single(protein)

    prof.summary(iterations=total_calls)
    prof.print_stats(sort_by=sort_by, top_n=40)
    prof.print_stats(sort_by="tottime", top_n=20)
    prof.print_callers(top_n=15)
    prof.save_stats()

    return result


def run_mass_batch() -> list[float]:
    """Run mass calculation on many peptides."""
    # Generate proteins and digest them
    proteins = [generate_random_protein(1_000) for _ in range(100)]
    all_peptides: list[str] = []

    for protein in proteins:
        peptides = pt.simple_digest(protein, cleave_on="KR", missed_cleavages=2, semi=False, min_len=6, max_len=50)
        all_peptides.extend(list(peptides))

    all_peptides = pt.modify(all_peptides, nterm_static="[+100]", cterm_static="[-50]", internal_variable={"M": "[+15.995]", "C": "[+57.021]"})

    print(f"Generated {len(all_peptides)} peptides for mass calculation")
    print(f"Sample peptides: {all_peptides[:3]}")

    # Calculate masses
    masses = []
    start_time = time.perf_counter()
    for peptide in all_peptides:
        mass = pt.mass(peptide, charge=0, ion_type="p")
        masses.append(mass)
    elapsed = time.perf_counter() - start_time

    print(f"\nCalculated masses for {len(all_peptides)} peptides")
    print(f"Took {elapsed:.2f}s ({elapsed / len(all_peptides) * 1000:.3f}ms per peptide)")

    return masses


def profile_mass_batch():
    """Profile mass calculation for many peptides."""
    with Profiler(name="mass_batch") as prof:
        result = run_mass_batch()

    prof.summary()
    prof.print_stats(sort_by="cumulative", top_n=40)
    prof.print_stats(sort_by="tottime", top_n=20)
    prof.save_stats()

    return result


def compare_implementations():
    """Compare different implementations or settings."""
    print("\n" + "=" * 80)
    print("COMPARATIVE BENCHMARKING")
    print("=" * 80 + "\n")

    # Test different protein lengths
    lengths = [10, 20, 50, 100]

    for length in lengths:
        protein = generate_random_protein(length)

        stats = benchmark(lambda: run_mass_single(protein), iterations=100, warmup=10)

        print_benchmark_results(f"Protein length {length}", stats, units="ms")


def analyze_hotspots():
    """Identify and analyze performance hotspots."""
    print("\n" + "=" * 80)
    print("HOTSPOT ANALYSIS")
    print("=" * 80 + "\n")

    with Profiler(name="hotspot_analysis") as prof:
        proteins = [generate_random_protein(20) for _ in range(10)]
        for protein in proteins:
            for _ in range(500):
                run_mass_single(protein)

    # Get stats object for detailed analysis
    stats = pstats.Stats(prof.profiler)
    stats.strip_dirs()
    stats.sort_stats("tottime")

    # Find functions taking >1% of total time
    print("Functions taking >1% of total time:")
    print("-" * 80)

    total_time = stats.total_tt
    for func, data in stats.stats.items():
        cc, nc, tt, ct, callers = data
        if tt / total_time > 0.01:  # >1% of time
            filename, line, func_name = func
            print(f"{func_name:50s} {tt:8.3f}s ({tt / total_time * 100:5.1f}%)")

    prof.save_stats()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Profile peptacular mass calculations")
    parser.add_argument("--mode", choices=["single", "batch", "compare", "hotspots", "all"], default="single", help="Profiling mode")
    parser.add_argument("--iterations", type=int, default=500, help="Iterations per protein (single mode)")
    parser.add_argument("--proteins", type=int, default=10, help="Number of proteins (single mode)")
    parser.add_argument("--length", type=int, default=20, help="Protein length (single mode)")

    args = parser.parse_args()

    if args.mode == "single" or args.mode == "all":
        print("=" * 80)
        print("SINGLE PROTEIN MASS CALCULATION")
        print("=" * 80)
        profile_mass_single(num_proteins=args.proteins, protein_length=args.length, iterations_per_protein=args.iterations)

    if args.mode == "batch" or args.mode == "all":
        print("\n" + "=" * 80)
        print("BATCH MASS CALCULATION")
        print("=" * 80)
        profile_mass_batch()

    if args.mode == "compare" or args.mode == "all":
        compare_implementations()

    if args.mode == "hotspots" or args.mode == "all":
        analyze_hotspots()

    print("\n✓ Profiling complete! Check the 'profile_results' directory for detailed outputs.")

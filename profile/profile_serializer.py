import cProfile
import pstats
import random
import statistics
import time
from io import StringIO
from typing import Callable

import peptacular as pt


def generate_random_protein(length: int = 500) -> str:
    """Generate a random protein sequence."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(length))


def benchmark(func: Callable, iterations: int = 1000, warmup: int = 100) -> dict:
    """Benchmark a function with warmup and multiple iterations."""
    # Warmup
    for _ in range(warmup):
        func()

    # Actual benchmark
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        func()
        end = time.perf_counter()
        times.append(end - start)

    return {
        "mean": statistics.mean(times),
        "median": statistics.median(times),
        "stdev": statistics.stdev(times) if len(times) > 1 else 0,
        "min": min(times),
        "max": max(times),
        "total": sum(times),
    }


def format_time(seconds: float) -> str:
    """Format time in appropriate units."""
    if seconds < 1e-6:
        return f"{seconds * 1e9:.2f} ns"
    elif seconds < 1e-3:
        return f"{seconds * 1e6:.2f} µs"
    elif seconds < 1:
        return f"{seconds * 1e3:.2f} ms"
    else:
        return f"{seconds:.2f} s"


def print_results(name: str, results: dict):
    """Print benchmark results nicely."""
    print(f"\n{name}:")
    print(f"  Mean:   {format_time(results['mean'])}")
    print(f"  Median: {format_time(results['median'])}")
    print(f"  Stdev:  {format_time(results['stdev'])}")
    print(f"  Min:    {format_time(results['min'])}")
    print(f"  Max:    {format_time(results['max'])}")


def profile_serialization():
    """Profile serialization for various peptide complexities."""

    # Test cases with increasing complexity
    test_cases = {
        "Simple peptide": "PEPTIDE",
        "With C-term mod": "PEPTIDE-[Amidated]",
        "With N-term mod": "[Acetyl]-PEPTIDE",
        "With internal mods": "PEP[+15.99]TI[+0.98]DE",
        "With multiple mods": "[Acetyl]-PEP[+15.99]TI[+0.98]DE-[Amidated]",
        "Complex with charge": "[Acetyl]-M[+15.99]PEPTIDEWITHM[+15.99]ETHIONINE-[Amidated]/2",
        "Long peptide": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL",
        "Long with mods": "[Acetyl]-M[+15.99]KTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEM[+15.99]PQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL-[Amidated]/3",
    }

    print("=" * 80)
    print("SERIALIZATION BENCHMARK")
    print("=" * 80)

    results = {}

    for name, sequence in test_cases.items():
        print(f"\nTesting: {name}")
        print(f"  Input: {sequence[:80]}{'...' if len(sequence) > 80 else ''}")

        # Parse once to get annotation object
        annotation = pt.parse(sequence)

        def test_serialize():
            return pt.serialize(annotation)

        result = benchmark(test_serialize, iterations=5000, warmup=500)
        results[name] = result
        print_results(f"  Serialization", result)

        # Verify output
        serialized = pt.serialize(annotation)
        print(f"  Output: {serialized[:80]}{'...' if len(serialized) > 80 else ''}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY - Mean Serialization Times")
    print("=" * 80)

    sorted_results = sorted(results.items(), key=lambda x: x[1]["mean"])
    for name, result in sorted_results:
        print(f"{name:<30} {format_time(result['mean'])}")

    return results


def profile_batch_serialization():
    """Profile serialization of many peptides (like in mod building)."""

    print("\n" + "=" * 80)
    print("BATCH SERIALIZATION BENCHMARK")
    print("=" * 80)

    # Generate realistic test data - peptides with modifications
    protein = generate_random_protein(1000)
    peptides = list(pt.simple_digest(protein, cleave_on="KR", missed_cleavages=2, semi=False, min_len=6, max_len=50))
    sequence_spans = peptides[:500]  # Limit to first 500 for speed
    sequences = [seq for seq, s in sequence_spans]

    print(f"\nGenerated {len(peptides)} peptides from digestion")
    # Add modifications to some peptides
    modified_peptides = []
    for peptide in sequences[:100]:  # Take first 100 for speed
        variants = list(
            pt.modify(
                peptide,
                internal_variable={"M": [15.995], "C": [57.021]},
                max_variable_mods=2,
            )
        )
        modified_peptides.extend(variants)

    print(f"Generated {len(modified_peptides)} modified peptide variants")

    # Parse all to get annotation objects
    annotations = []
    parse_start = time.time()
    for peptide in modified_peptides:
        ann = pt.parse(peptide)
        annotations.append(ann)
    parse_time = time.time() - parse_start

    print(f"Parsed {len(annotations)} annotations in {parse_time:.3f}s")
    print(f"Average parse time: {format_time(parse_time / len(annotations))}")

    # Benchmark serialization
    serialize_start = time.time()
    serialized = []
    for ann in annotations:
        serialized.append(pt.serialize(ann))
    serialize_time = time.time() - serialize_start

    print(f"\nSerialized {len(annotations)} annotations in {serialize_time:.3f}s")
    print(f"Average serialize time: {format_time(serialize_time / len(annotations))}")
    print(f"Throughput: {len(annotations) / serialize_time:.0f} serializations/sec")

    return annotations, serialized


def profile_detailed_single():
    """Detailed profiling of a single complex serialization."""

    print("\n" + "=" * 80)
    print("DETAILED PROFILING - Complex Peptide Serialization")
    print("=" * 80)

    # Create a complex peptide with many features
    complex_peptide = "[Acetyl]-M[+15.99]PEPTIDEWITHM[+15.99]ETHIONINEANDC[+57.02]YSTEINE-[Amidated]/2"

    print(f"\nProfiling: {complex_peptide}")

    # Parse once
    annotation = pt.parse(complex_peptide)

    # Profile serialization
    profiler = cProfile.Profile()
    profiler.enable()

    # Run many iterations to get good profiling data
    for _ in range(10000):
        pt.serialize(annotation)

    profiler.disable()

    # Print stats
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
    ps.print_stats(40)
    print(s.getvalue())


def profile_edge_cases():
    """Test serialization of edge cases."""

    print("\n" + "=" * 80)
    print("EDGE CASES BENCHMARK")
    print("=" * 80)

    edge_cases = {
        "Empty modifications": "PEPTIDE",
        "Only N-term": "[Acetyl]-PEPTIDE",
        "Only C-term": "PEPTIDE-[Amidated]",
        "Only internal": "PEPT[+15.99]IDE",
        "Multiple same position": "M[+15.99][+0.98]PEPTIDE",
        "All mod types": "[Acetyl]-M[+15.99]PEPTIDE-[Amidated]/2",
        "Very long sequence": "A" * 1000,
        "Many mods": "M[+15.99]" * 50,
    }

    for name, sequence in edge_cases.items():
        try:
            annotation = pt.parse(sequence)

            def test_serialize():
                return pt.serialize(annotation)

            result = benchmark(test_serialize, iterations=1000, warmup=100)
            print(f"\n{name}:")
            print(f"  Length: {len(sequence)}")
            print(f"  Time: {format_time(result['mean'])}")

            # Verify round-trip
            serialized = pt.serialize(annotation)
            reparsed = pt.parse(serialized)
            reserialized = pt.serialize(reparsed)

            if serialized == reserialized:
                print(f"  ✓ Round-trip successful")
            else:
                print(f"  ✗ Round-trip FAILED")
                print(f"    Original:  {serialized}")
                print(f"    Reparsed:  {reserialized}")
        except Exception as e:
            print(f"\n{name}: ERROR - {e}")


if __name__ == "__main__":
    # Run all benchmarks
    profile_serialization()
    profile_batch_serialization()
    profile_edge_cases()
    profile_detailed_single()

    print("\n" + "=" * 80)
    print("PROFILING COMPLETE")
    print("=" * 80)

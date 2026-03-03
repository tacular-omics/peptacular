import cProfile
import pstats
import random
import time
from io import StringIO

import peptacular as pt


def generate_random_protein(length: int = 500) -> str:
    """Generate a random protein sequence."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(length))


def run_mass():
    """Run mass calculation on many peptides."""
    # Generate proteins and digest them
    proteins = [generate_random_protein(1_000) for _ in range(100)]

    all_peptides: list[str] = []
    for protein in proteins:
        peptides = pt.simple_digest(protein, cleave_on="KR", missed_cleavages=2, semi=False, min_len=6, max_len=50)
        all_peptides.extend(list(peptides))

    all_peptides = pt.modify(all_peptides, nterm_static="[+100]", cterm_static="[-50]", internal_variable={"M": "[+15.995]", "C": "[+57.021]"})

    print(f"Generated {len(all_peptides)} peptides for mass calculation")

    print(all_peptides[:5])  # Show first 5 peptides

    # Now profile the mass calculation
    start_time = time.time()
    for peptide in all_peptides:
        masses = pt.mass(peptide, charge=0, ion_type="p")
    elapsed = time.time() - start_time

    print(f"Calculated masses for {len(all_peptides)} peptides")
    print(f"Took {elapsed:.2f} seconds ({elapsed / len(all_peptides) * 1000:.3f} ms per peptide)")
    return masses


def run_mass_single(protein):
    """Run mass calculation on peptides from a single protein."""

    protein = "{Oxidation}M[+100]" + protein + "K-[Oxidation]"  # Ensure cleavable ends

    start_time = time.time()
    masses = pt.mass(protein, charge=0, ion_type="p")
    elapsed = time.time() - start_time
    return masses


def profile_mass_single():
    """Profile mass calculation for single protein peptides (no multiprocessing)."""
    profiler = cProfile.Profile()
    profiler.enable()

    proteins = [generate_random_protein(20) for _ in range(10)]

    for protein in proteins:
        for i in range(500):
            result = run_mass_single(protein)
    profiler.disable()

    # Print stats
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
    ps.print_stats(30)  # Top 30 functions
    print(s.getvalue())
    return result


def profile_mass_multi():
    """Profile mass calculation for many peptides (with multiprocessing)."""
    profiler = cProfile.Profile()
    profiler.enable()
    result = run_mass()
    profiler.disable()

    # Print stats
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
    ps.print_stats(30)  # Top 30 functions
    print(s.getvalue())
    return result


if __name__ == "__main__":
    print("=" * 60)
    print("SINGLE PROTEIN MASS CALCULATION (no multiprocessing)")
    print("=" * 60)
    profile_mass_single()

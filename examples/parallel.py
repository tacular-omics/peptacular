"""
Parallel Execution Examples
============================
Examples of using parallel processing with peptacular sequence functions.
Most functions in the sequence module accept lists and automatically use
multiprocessing for better performance on large datasets.
"""

import peptacular as pt
import time

# ============================================================================
# BASIC PARALLEL EXECUTION (must be run through __main__)
# ============================================================================

def run():
    print("=" * 60)
    print("BASIC PARALLEL EXECUTION")
    print("=" * 60)

    # Create a list of peptide sequences
    peptides = [
        'PEPTIDE',
        'PATRICK',
        'TYLER',
        'GARRETT'
    ]

    # Parse multiple sequences in parallel
    # Pass a list of strings - automatically uses multiprocessing
    annotations = pt.parse(peptides)
    print(f"\nParsed {len(annotations)} peptides:")
    for annot in annotations:
        print(f"  {annot.serialize()}")

    # ============================================================================
    # MASS CALCULATIONS IN PARALLEL
    # ============================================================================

    print("\n" + "=" * 60)
    print("MASS CALCULATIONS")
    print("=" * 60)

    # Calculate masses for multiple peptides
    masses = pt.mass(peptides)
    print("\nMasses:")
    for pep, mass in zip(peptides, masses):
        print(f"  {pep}: {mass:.4f} Da")

    # With charge states
    masses_charged = pt.mass(peptides, charge=2)
    print("\nMasses at +2 charge:")
    for pep, mass in zip(peptides, masses_charged):
        print(f"  {pep}: {mass:.4f} Da")

    # m/z calculations
    mz_values = pt.mz(peptides, charge=2)
    print("\nm/z at +2 charge:")
    for pep, mz in zip(peptides, mz_values):
        print(f"  {pep}: {mz:.4f}")

    # ============================================================================
    # COMPOSITION IN PARALLEL
    # ============================================================================

    print("\n" + "=" * 60)
    print("COMPOSITION CALCULATIONS")
    print("=" * 60)

    # Get compositions for multiple peptides
    compositions = pt.comp(peptides)
    print("\nCompositions:")
    for pep, comp in zip(peptides, compositions):
        comp_str = {str(elem): count for elem, count in comp.items()}
        print(f"  {pep}: {comp_str}")

    # ============================================================================
    # PHYSICOCHEMICAL PROPERTIES
    # ============================================================================

    print("\n" + "=" * 60)
    print("PHYSICOCHEMICAL PROPERTIES")
    print("=" * 60)

    # Calculate hydrophobicity for all peptides
    hydrophobicity = pt.hydrophobicity(peptides)
    print("\nHydrophobicity:")
    for pep, hydro in zip(peptides, hydrophobicity):
        print(f"  {pep}: {hydro:.3f}")

    # Calculate pI values
    pi_values = pt.pi(peptides)
    print("\nIsoelectric points:")
    for pep, pi in zip(peptides, pi_values):
        print(f"  {pep}: {pi:.2f}")

    # Calculate aromaticity
    aromaticity = pt.aromaticity(peptides)
    print("\nAromaticity:")
    for pep, arom in zip(peptides, aromaticity):
        print(f"  {pep}: {arom:.3f}")

    # ============================================================================
    # USING ANNOTATION OBJECTS
    # ============================================================================

    print("\n" + "=" * 60)
    print("USING ANNOTATION OBJECTS")
    print("=" * 60)

    # Can also pass lists of ProFormaAnnotation objects
    modified_peptides = [
        pt.parse('[Acetyl]-PEPTIDE'),
        pt.parse('PEM[Oxidation]TIDE'),
        pt.parse('SEQS[Phospho]UENCE/2')
    ]

    # Calculate masses from annotations
    annot_masses = pt.mass(modified_peptides)
    print("\nModified peptide masses:")
    for annot, mass in zip(modified_peptides, annot_masses):
        print(f"  {annot.serialize()}: {mass:.4f} Da")

    # ============================================================================
    # CUSTOM PROPERTY CALCULATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("CUSTOM PROPERTY CALCULATIONS")
    print("=" * 60)

    # Calculate custom properties in parallel
    custom_props = pt.calc_property(
        peptides,
        scale=pt.HydrophobicityScale.KYTE_DOOLITTLE,
        aggregation_method='avg'
    )
    print("\nKyte-Doolittle hydrophobicity (average):")
    for pep, prop in zip(peptides, custom_props):
        print(f"  {pep}: {prop:.3f}")

    # ============================================================================
    # PERFORMANCE COMPARISON
    # ============================================================================

    print("\n" + "=" * 60)
    print("PERFORMANCE COMPARISON")
    print("=" * 60)

    # Create larger dataset for timing
    large_dataset = peptides * 10000  # 40,000 peptides

    # Time parallel execution
    start = time.time()
    _ = pt.mass(large_dataset)
    parallel_time = time.time() - start

    print(f"\nProcessed {len(large_dataset)} peptides:")
    print(f"  Parallel execution: {parallel_time:.3f} seconds")
    print(f"  Average per peptide: {parallel_time/len(large_dataset)*1000:.2f} ms")

    # Time serial execution for comparison
    start = time.time()
    _ = [pt.mass(pep) for pep in large_dataset]
    serial_time = time.time() - start
    print(f"  Serial execution: {serial_time:.3f} seconds")
    print(f"  Average per peptide: {serial_time/len(large_dataset)*1000:.2f} ms")

    print("\n" + "=" * 60)

if __name__ == "__main__":
    run()
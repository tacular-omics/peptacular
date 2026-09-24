"""
Sequence Property Calculations
===============================
Examples of calculating physicochemical and structural properties of peptides.
Note: These calculations use only the amino acid sequence; modifications are not considered.
"""

import peptacular as pt


def run():
    # Parse a test peptide
    annot = pt.parse('PEPTIDE')

    # ============================================================================
    # SIMPLE PHYSICOCHEMICAL PROPERTIES
    # ============================================================================

    print("=" * 60)
    print("PHYSICOCHEMICAL PROPERTIES")
    print("=" * 60)

    # These properties return single float values
    print(f"Sequence: {annot}")
    print(f"Hydrophobicity: {annot.prop.hydrophobicity:.3f}")
    print(f"Flexibility: {annot.prop.flexibility:.3f}")
    print(f"Hydrophilicity: {annot.prop.hydrophilicity:.3f}")
    print(f"Surface accessibility: {annot.prop.surface_accessibility:.3f}")
    print(f"Polarity: {annot.prop.polarity:.3f}")
    print(f"Aromaticity: {annot.prop.aromaticity:.3f}")
    print(f"Isoelectric point (pI): {annot.prop.pi:.2f}")
    print(f"HPLC retention: {annot.prop.hplc:.3f}")
    print(f"Refractivity: {annot.prop.refractivity:.3f}")

    # ============================================================================
    # STRUCTURAL PROPERTIES
    # ============================================================================

    print("\n" + "=" * 60)
    print("STRUCTURAL PROPERTIES")
    print("=" * 60)

    # Secondary structure percentages
    print(f"Alpha helix: {annot.prop.alpha_helix_percent:.1f}%")
    print(f"Beta sheet: {annot.prop.beta_sheet_percent:.1f}%")
    print(f"Beta turn: {annot.prop.beta_turn_percent:.1f}%")
    print(f"Coil: {annot.prop.coil_percent:.1f}%")

    # Predicted secondary structure using different methods
    ss_dr = annot.prop.secondary_structure(pt.SecondaryStructureMethod.DELEAGE_ROUX)
    print("\nSecondary structure (Deleage-Roux method):")
    print(f"  Alpha helix: {ss_dr['alpha_helix']:.1f}%")
    print(f"  Beta sheet: {ss_dr['beta_sheet']:.1f}%")
    print(f"  Beta turn: {ss_dr['beta_turn']:.1f}%")
    print(f"  Coil: {ss_dr['coil']:.1f}%")

    # ============================================================================
    # COMPOSITION-BASED PROPERTIES
    # ============================================================================

    print("\n" + "=" * 60)
    print("COMPOSITION PROPERTIES")
    print("=" * 60)

    # Amino acid composition
    proline_pct = annot.prop.aa_property_percentage('P')
    acidic_pct = annot.prop.aa_property_percentage('DE')  # D and E
    basic_pct = annot.prop.aa_property_percentage('KR')   # K and R
    print(f"Proline content: {proline_pct:.1f}%")
    print(f"Acidic residues (D, E): {acidic_pct:.1f}%")
    print(f"Basic residues (K, R): {basic_pct:.1f}%")

    # Charge at different pH values
    print(f"\nNet charge at pH 7.0: {annot.prop.charge_at_ph(7.0):.2f}")
    print(f"Net charge at pH 3.0: {annot.prop.charge_at_ph(3.0):.2f}")
    print(f"Net charge at pH 11.0: {annot.prop.charge_at_ph(11.0):.2f}")

    # ============================================================================
    # CUSTOM PROPERTY CALCULATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("CUSTOM PROPERTY CALCULATIONS")
    print("=" * 60)

    # --- Basic calculation with default options ---
    prop = annot.prop.calc_property(
        scale=pt.HydrophobicityScale.ABRAHAM_LEO,
        missing_aa_handling=pt.MissingAAHandling.ERROR,  # default
        aggregation_method=pt.AggregationMethod.SUM,     # default
        normalize=False,                                  # default
        weighting_scheme=pt.WeightingMethods.UNIFORM,    # default
        min_weight=0.0,                                   # default
        max_weight=1.0,                                   # default
    )
    print(f"Abraham-Leo hydrophobicity (sum): {prop:.2f}")

    # --- Using string identifiers ---
    prop_avg = annot.prop.calc_property(
        scale="deleage_roux_alpha_helix",
        missing_aa_handling="avg",
        aggregation_method="avg"
    )
    print(f"Alpha helix propensity (avg): {prop_avg:.3f}")

    # --- Custom scale dictionary ---
    custom_scale = {
        'A': 1.0, 'C': 2.0, 'D': 3.0, 'E': 4.0,
        'F': 5.0, 'G': 6.0, 'H': 7.0, 'I': 8.0,
        'K': 9.0, 'L': 10.0, 'M': 11.0, 'N': 12.0,
        'P': 13.0, 'Q': 14.0, 'R': 15.0, 'S': 16.0,
        'T': 17.0, 'V': 18.0, 'W': 19.0, 'Y': 20.0
    }
    custom_prop = annot.prop.calc_property(scale=custom_scale, missing_aa_handling="avg")
    print(f"Custom scale (sum): {custom_prop:.2f}")

    # ============================================================================
    # AVAILABLE OPTIONS FOR calc_property()
    # ============================================================================

    print("\n" + "=" * 60)
    print("CALC_PROPERTY OPTIONS")
    print("=" * 60)

    """
    [Scale]
    - Use built-in scale enums (e.g., HydrophobicityScale.ABRAHAM_LEO)
    - Use scale name as string (e.g., "abraham_leo")
    - Provide custom dict (e.g., {'A': 1.0, 'C': 2.0, ...})
    - ~50 built-in scales available

    [missing_aa_handling]
    - 'avg': Use average of known values
    - 'min': Use minimum of known values
    - 'max': Use maximum of known values
    - 'median': Use median of known values
    - 'error': Raise error (default)
    - 'zero': Use 0.0
    - 'skip': Skip missing amino acids

    [aggregation_method]
    - 'sum': Sum of amino acid values (default)
    - 'avg': Average of amino acid values

    [normalize]
    - True: Normalize each AA's property value to [0, 1] before aggregation
    - False: Use raw values (default)

    [weighting_scheme]
    - 'uniform': All positions weighted equally (default)
    - 'linear': Linear weighting across sequence
    - 'exponential': Exponential weighting
    - 'gaussian': Gaussian weighting
    - 'sigmoid': Sigmoid weighting
    - 'cosine': Cosine weighting
    - 'sinusoidal': Sinusoidal weighting

    [min_weight, max_weight]
    - Define weight range (default: 0.0 to 1.0)
    """

    # ============================================================================
    # SLIDING WINDOW CALCULATIONS
    # ============================================================================

    print("=" * 60)
    print("SLIDING WINDOW CALCULATIONS")
    print("=" * 60)

    # Calculate property over sliding windows
    windows = annot.prop.property_windows(
        scale=pt.HydrophobicityScale.ABRAHAM_LEO,
        window_size=4,
        missing_aa_handling=pt.MissingAAHandling.ERROR,
        aggregation_method=pt.AggregationMethod.SUM,
        normalize=False,
        weighting_scheme=pt.WeightingMethods.UNIFORM,
        min_weight=0.0,
        max_weight=1.0,
    )
    print("\nWindow size 4 (overlapping):")
    print(f"  Values: {[f'{v:.2f}' for v in windows]}")
    print(f"  Number of windows: {len(windows)}")

    # Different window size
    windows_large = annot.prop.property_windows(
        scale=pt.HydrophobicityScale.ABRAHAM_LEO,
        window_size=3
    )
    print("\nWindow size 3:")
    print(f"  Values: {[f'{v:.2f}' for v in windows_large]}")

    # ============================================================================
    # PARTITIONED WINDOW CALCULATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("PARTITIONED WINDOW CALCULATIONS")
    print("=" * 60)

    # Divide sequence into fixed number of non-overlapping partitions
    partitions = annot.prop.property_partitions(
        scale=pt.HydrophobicityScale.ABRAHAM_LEO,
        num_windows=3,
        aa_overlap=0,
        missing_aa_handling=pt.MissingAAHandling.ERROR,
        aggregation_method=pt.AggregationMethod.SUM,
        normalize=False,
        weighting_scheme=pt.WeightingMethods.UNIFORM,
        min_weight=0.0,
        max_weight=1.0,
    )
    print("\n3 partitions (no overlap):")
    print(f"  Values: {[f'{v:.2f}' for v in partitions]}")

    # With overlap between partitions
    partitions_overlap = annot.prop.property_partitions(
        scale=pt.HydrophobicityScale.ABRAHAM_LEO,
        num_windows=3,
        aa_overlap=1
    )
    print("\n3 partitions (1 AA overlap):")
    print(f"  Values: {[f'{v:.2f}' for v in partitions_overlap]}")

    # ============================================================================
    # PRACTICAL EXAMPLES
    # ============================================================================

    print("\n" + "=" * 60)
    print("PRACTICAL EXAMPLES")
    print("=" * 60)

    # Example: Hydrophobicity profile for transmembrane prediction
    tm_peptide = pt.parse('LFGAIAGFIENGWEGMIDG')
    tm_windows = tm_peptide.prop.property_windows(
        scale=pt.HydrophobicityScale.KYTE_DOOLITTLE,
        window_size=9
    )
    print(f"\nTransmembrane peptide: {tm_peptide}")
    print("Kyte-Doolittle hydrophobicity profile (window=9):")
    for i, val in enumerate(tm_windows):
        print(f"  Position {i+1}: {val:.2f}")

    # Example: Charge distribution analysis
    charged_peptide = pt.parse('PKDEPKDE')
    charge_partitions = charged_peptide.prop.property_partitions(
        scale={'K': 1, 'R': 1, 'D': -1, 'E': -1},  # Simple charge scale
        num_windows=4,
        aa_overlap=0,
        missing_aa_handling='zero'
    )
    print(f"\nCharged peptide: {charged_peptide}")
    print("Charge distribution (4 regions):")
    for i, val in enumerate(charge_partitions):
        print(f"  Region {i+1}: {val:+.1f}")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    run()

"""
Isotopic Distribution Calculations
===================================
Examples of calculating isotopic distributions from ProForma annotations.
"""

import peptacular as pt


def run():
    # Parse a simple peptide sequence
    annot = pt.parse("PEPTIDE")

    # ============================================================================
    # BASIC ISOTOPIC DISTRIBUTION
    # ============================================================================

    print("=" * 60)
    print("BASIC ISOTOPIC DISTRIBUTION")
    print("=" * 60)

    # --- Default Distribution ---
    # Returns list of IsotopicData with mass, neutron_count, and abundance
    # Abundances normalized so max peak = 1.0
    dist = annot.isotopic_distribution()
    print(f"\nPeptide: {annot.serialize()}")
    print(f"Monoisotopic mass: {annot.mass():.3f} Da")
    print("Default isotopic distribution:")
    for iso in dist:
        print(
            f"  mass: {iso.mass:>8.3f} Da, abundance: {iso.abundance:>6.3f}, neutrons: {iso.neutron_count}"
        )

    # --- Control Number of Isotopes ---
    dist_limited = annot.isotopic_distribution(max_isotopes=3)
    print("\nLimited to the first 3 nominal isotope positions:")
    for iso in dist_limited:
        print(f"  mass: {iso.mass:>8.3f} Da, abundance: {iso.abundance:>6.3f}")

    # --- Abundance Threshold ---
    # Retain peaks through the last one meeting the threshold, including weaker leading peaks.
    dist_filtered = annot.isotopic_distribution(min_abundance_threshold=0.05)
    print("\nEnvelope through the last peak at least 5% of the maximum:")
    for iso in dist_filtered:
        print(f"  mass: {iso.mass:>8.3f} Da, abundance: {iso.abundance:>6.3f}")

    # --- Neutron Offsets ---
    # Every peak provides both its center mass and nominal neutron count.
    print("\nNominal neutron offsets:")
    for iso in dist:
        print(f"  neutron offset: {iso.neutron_count:>3}, abundance: {iso.abundance:>6.3f}")

    # ============================================================================
    # AGGREGATED CENTER MASSES
    # ============================================================================

    print("\n" + "=" * 60)
    print("AGGREGATED CENTER MASSES")
    print("=" * 60)

    # BRAIN produces a probability-weighted center mass per nominal isotope peak.
    # Formatting controls displayed precision, without changing the calculation.
    print("\nCenter masses displayed to 5 decimal places:")
    for iso in dist[:3]:
        print(f"  mass: {iso.mass:.5f} Da, abundance: {iso.abundance:>6.3f}")

    # ============================================================================
    # COMBINING WITH COMP PARAMETERS
    # ============================================================================

    print("\n" + "=" * 60)
    print("COMBINING WITH COMP PARAMETERS")
    print("=" * 60)

    # isotopic_distribution() accepts same parameters as comp()
    # Combine charge, isotopes, losses, and ion type
    dist_combined = annot.isotopic_distribution(
        ion_type="y", charge=2, isotopes=1, deltas={"H2O": 1}
    )
    print("\ny-ion, +2 charge, +1 13C, -H2O:")
    for iso in dist_combined[:4]:
        print(f"  m/z: {iso.mass / 2:>8.3f}, abundance: {iso.abundance:>6.3f}")


if __name__ == "__main__":
    run()

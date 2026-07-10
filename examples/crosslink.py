"""Cross-linked and Branched Peptidoform Ions
==============================================
Examples of parsing and analysing cross-linked peptides (ProForma 2.1 §9.2.2, §9.3).

Chains joined by ``//`` share a single charge state and parse into a
``MultiProFormaAnnotation`` container; a single peptidoform (including an intra-chain
cross-link) parses into an ordinary ``ProFormaAnnotation``.
"""

import peptacular as pt


def run():
    # ============================================================================
    # INTRA-CHAIN CROSS-LINK (both ends on one chain -> single annotation)
    # ============================================================================

    print("=" * 60)
    print("INTRA-CHAIN CROSS-LINK (§9.2.1)")
    print("=" * 60)

    intra = pt.parse("EVTK[XLMOD:02001#XL1]LEK[#XL1]SEFD")
    print(f"Type:  {type(intra).__name__}")
    print(f"Mass:  {intra.mass():.4f} Da  (DSBU linker mass counted once)")

    # ============================================================================
    # INTER-CHAIN CROSS-LINK (ends span two chains -> MultiProFormaAnnotation)
    # ============================================================================

    print("\n" + "=" * 60)
    print("INTER-CHAIN CROSS-LINK (§9.2.2)")
    print("=" * 60)

    ion = pt.parse("EVTK[XLMOD:02001#XL1]LE//AK[#XL1]ENLYFQ/3")
    print(f"Type:         {type(ion).__name__}")
    print(f"Chains:       {len(ion)}")
    for i, chain in enumerate(ion):
        print(f"  chain {i}:    {chain.serialize()}")
    print(f"Shared charge: {ion.charge}")
    print(f"Neutral mass:  {ion.neutral_mass():.4f} Da")
    print(f"Ion mass:      {ion.mass():.4f} Da")
    print(f"m/z:           {ion.mz():.4f}")
    print(f"Serialize:     {ion.serialize()}")

    # ============================================================================
    # BRANCHED PEPTIDE
    # ============================================================================

    print("\n" + "=" * 60)
    print("BRANCHED PEPTIDE (§9.3)")
    print("=" * 60)

    branch = pt.parse("ED[MOD:00093#BRANCH]//D[#BRANCH]ATR/1")
    print(f"Serialize: {branch.serialize()}")
    print(f"Mass:      {branch.mass():.4f} Da")

    # ============================================================================
    # CROSS-LINK LABEL VALIDATION
    # ============================================================================

    print("\n" + "=" * 60)
    print("CROSS-LINK VALIDATION (validate=True)")
    print("=" * 60)

    # A dangling back-reference (#XL1 referenced but never defined) is rejected.
    try:
        pt.parse("PEK[#XL1]TIDE//AKENLYFQ", validate=True)
    except ValueError as e:
        print(f"Rejected: {e}")


if __name__ == "__main__":
    run()

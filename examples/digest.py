"""
Protein Digestion Examples
===========================
Simple examples of in-silico enzymatic digestion using ProForma.
All digestion methods return Span objects (start, end, missed_cleavages).
Use annotation[span] to get the actual peptide.
"""

import re

import peptacular as pt


def run():
    # ============================================================================
    # SIMPLE DIGESTION (AA Based)
    # ============================================================================

    protein = pt.parse("[Amidated]-PEPTIDEKPEPTIDERPEPT[Phospho]IDER-[+57]")

    print("=" * 60)
    print("SIMPLE DIGESTION (AA BASED)")
    print("=" * 60)
    print(f"Protein: {protein}\n")

    # Basic trypsin-like digestion
    print("Trypsin-like (cleave after K/R):")
    for span in protein.simple_digest_spans(cleave_on="KR"):
        peptide = protein[span]
        print(f"  {peptide.serialize()} - span: {span}")

    # With restrictions
    print("\nWith restrictions (cleave after K/R, but not before N or after P):")
    for span in protein.simple_digest_spans(cleave_on="KR", restrict_before="N", restrict_after="P", cterminal=True):
        print(f"  {protein[span].serialize()}")

    # ============================================================================
    # DIGESTION (REGEX BASED)
    # ============================================================================

    print("\n" + "=" * 60)
    print("DIGESTION (REGEX BASED)")
    print("=" * 60)

    # Using predefined enzyme enum
    print("\nUsing Proteases enum:")
    for span in protein.digest_spans(pt.Proteases.TRYPSIN):
        print(f"  {protein[span].serialize()}")

    # Using enzyme string
    print("\nUsing enzyme string 'trypsin':")
    for span in protein.digest_spans("trypsin"):
        print(f"  {protein[span].serialize()}")

    # Custom regex: pass a compiled pattern (a plain string is only a protease name)
    print("\nCustom regex (cleave after A or E):")
    for span in protein.digest_spans(re.compile("(?<=[AE])")):
        print(f"  {protein[span].serialize()}")

    # ============================================================================
    # CLEAVAGE SITES
    # ============================================================================

    print("\n" + "=" * 60)
    print("CLEAVAGE SITES")
    print("=" * 60)

    print("\nCleavage positions for trypsin (after K/R):")
    sites = list(
        protein.simple_cleavage_sites(
            cleave_on="KR",
            restrict_after="P",
            restrict_before="N",
            cterminal=True,
        )
    )
    print(f"  Sites: {sites}")
    print(f"  Sequence: {protein.sequence}")
    print(f"            {''.join('^' if i in sites else ' ' for i in range(len(protein.sequence)))}")

    print("\nCleavage positions for included trypsin regex:")
    # can also use Proteases.TRYPSIN or custom regex
    sites_regex = list(protein.cleavage_sites("trypsin"))
    print(f"  Sites: {sites_regex}")

    # ============================================================================
    # MISSED CLEAVAGES
    # ============================================================================

    print("\n" + "=" * 60)
    print("MISSED CLEAVAGES")
    print("=" * 60)

    print("\nWith 1 missed cleavage:")
    for span in protein.digest_spans("trypsin", missed_cleavages=1):
        print(f"  {protein[span].serialize()}")

    # ============================================================================
    # LENGTH FILTERING
    # ============================================================================

    print("\n" + "=" * 60)
    print("LENGTH FILTERING")
    print("=" * 60)

    print("\nPeptides between 7-15 amino acids:")
    for span in protein.digest_spans("trypsin", min_len=7, max_len=15):
        peptide = protein[span]
        print(f"  {peptide.serialize()} (length: {len(peptide)})")

    # ============================================================================
    # SEMI-ENZYMATIC DIGESTION
    # ============================================================================

    print("\n" + "=" * 60)
    print("SEMI-ENZYMATIC")
    print("=" * 60)

    print("\nSemi-enzymatic (one end must be enzymatic):")
    for span in protein.digest_spans("trypsin", semi=True, min_len=5, max_len=10):
        print(f"  {protein[span].serialize()}")


if __name__ == "__main__":
    run()

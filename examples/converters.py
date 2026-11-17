"""
Sequence Format Conversion Examples
====================================
Examples of converting peptide sequences from other tools (IP2, DIANN, MaxQuant,
Casanovo, MS2PIP) to ProForma 2.1 format. All conversion functions support parallel execution.
"""

import peptacular as pt


def run():
    # ============================================================================
    # IP2 SEQUENCE CONVERSION
    # ============================================================================

    print("=" * 60)
    print("IP2 SEQUENCE CONVERSION")
    print("=" * 60)

    # Basic IP2 format: K.SEQUENCE.K
    ip2_seq = "K.PEPTIDE.K"
    proforma = pt.convert_ip2_sequence(ip2_seq)
    print(f"IP2: {ip2_seq}")
    print(f"ProForma: {proforma}\n")

    # ============================================================================
    # DIANN SEQUENCE CONVERSION
    # ============================================================================

    print("\n" + "=" * 60)
    print("DIANN SEQUENCE CONVERSION")
    print("=" * 60)

    # With modification
    diann_mod = "_YMGTLRGC[Carbamidomethyl]LLRLYHD_"
    proforma_mod = pt.convert_diann_sequence(diann_mod)
    print(f"DIANN with mod: {diann_mod}")
    print(f"ProForma: {proforma_mod}\n")

    # ============================================================================
    # MAXQUANT SEQUENCE CONVERSION
    # ============================================================================

    print("\n" + "=" * 60)
    print("MAXQUANT SEQUENCE CONVERSION")
    print("=" * 60)

    # The Modified sequence column of evidence.txt and msms.txt: the sequence inside
    # underscores, each modification in round brackets after its residue.
    maxquant_mod = "_(ac)SM(ox)K_"
    print(f"MaxQuant with mods: {maxquant_mod}")
    print(f"ProForma: {pt.convert_maxquant_sequence(maxquant_mod)}")

    # Full modification names, nested parentheses and all (MaxQuant 1.6.7.0 onward)
    maxquant_full = "_ALTDM(Oxidation (M))PQM(Oxidation (M))R_"
    print(f"MaxQuant full names: {maxquant_full}")
    print(f"ProForma: {pt.convert_maxquant_sequence(maxquant_full)}")

    # Writing a MaxQuant sequence back out
    maxquant_annot = pt.parse("PEP[Phospho]TIDE-[Amidated]")
    print(f"\nProForma: {maxquant_annot.serialize()}")
    print(f"MaxQuant: {maxquant_annot.to_maxquant()}\n")

    # ============================================================================
    # CASANOVO SEQUENCE CONVERSION
    # ============================================================================

    print("\n" + "=" * 60)
    print("CASANOVO SEQUENCE CONVERSION")
    print("=" * 60)

    # Complex example
    casanovo_complex = "+43.006P+100EPTIDE"
    proforma_complex = pt.convert_casanovo_sequence(casanovo_complex)
    print(f"Casanovo complex: {casanovo_complex}")
    print(f"ProForma: {proforma_complex}")

    # Parse Casanovo format using annotation method
    casanovo_annot = pt.ProFormaAnnotation.from_casanovo("+43.006PEPTIDE")
    print(f"\nCasanovo (annotation method): {casanovo_annot.serialize()}")
    print(f"  Mass: {casanovo_annot.mass():.4f} Da")

    # ============================================================================
    # MS2PIP FORMAT CONVERSION
    # ============================================================================

    print("\n" + "=" * 60)
    print("MS2PIP FORMAT CONVERSION")
    print("=" * 60)

    # Convert TO MS2PIP format
    pf_annot = pt.parse("[Acetyl]-PEM[Oxidation]TIDE")
    unmod_seq, mod_str = pf_annot.to_ms2_pip()
    print(f"\nProForma: {pf_annot.serialize()}")
    print(f"MS2PIP sequence: {unmod_seq}")
    print(f"MS2PIP mods: {mod_str}")

    # Convert FROM MS2PIP format
    ms2pip_annot = pt.ProFormaAnnotation.from_ms2_pip(
        sequence="PEPTIDE", mod_str="0|Acetyl|3|Oxidation"
    )
    print(f"\nMS2PIP -> ProForma: {ms2pip_annot.serialize()}")

    # With static modifications
    ms2pip_static = pt.ProFormaAnnotation.from_ms2_pip(
        sequence="PEPTIDE", mod_str="0|Acetyl", static_mods={"C": "Carbamidomethyl"}
    )
    print(f"MS2PIP with static mods: {ms2pip_static.serialize()}")


if __name__ == "__main__":
    run()

"""
Sequence Format Conversion Examples
====================================
Examples of converting peptide sequences from other tools (IP2, DIANN, MaxQuant, Skyline,
Casanovo, MS2PIP, MSFragger) to ProForma 2.1 format. All conversion functions support
parallel execution.
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
    # SKYLINE SEQUENCE CONVERSION
    # ============================================================================

    print("\n" + "=" * 60)
    print("SKYLINE SEQUENCE CONVERSION")
    print("=" * 60)

    # Skyline writes each modification in square brackets after its residue, as a signed
    # mass difference, a name, or (unimod:N) in round brackets.
    skyline_mass = "AAM[15.9949]Q[-0.98]R"
    print(f"Skyline with deltas: {skyline_mass}")
    print(f"ProForma: {pt.convert_skyline_sequence(skyline_mass)}")

    skyline_unimod = "GSMAFVPK(unimod:259)"
    print(f"Skyline with UniMod: {skyline_unimod}")
    print(f"ProForma: {pt.convert_skyline_sequence(skyline_unimod)}")

    # Writing a Skyline sequence back out
    skyline_annot = pt.parse("PEPC[Carbamidomethyl]IDE")
    print(f"\nProForma: {skyline_annot.serialize()}")
    print(f"Skyline: {skyline_annot.to_skyline()}\n")

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
    # MSFRAGGER / FRAGPIPE CONVERSION
    # ============================================================================

    print("\n" + "=" * 60)
    print("MSFRAGGER / FRAGPIPE CONVERSION")
    print("=" * 60)

    # FragPipe's psm.tsv carries two columns: the stripped Peptide and the entries in
    # Assigned Modifications. This row is from compomics/psm_utils' example file.
    peptide, assigned = "LHISPSNMTNQNTPEYMEK", "8M(15.9949), C-term(-0.9840)"
    print(f"Peptide: {peptide}")
    print(f"Assigned Modifications: {assigned}")
    print(f"ProForma: {pt.from_msfragger((peptide, assigned))}")

    # Writing the pair back out. The entries are sorted as strings and every mass gets four
    # decimals, which is what Philosopher writes.
    annot = pt.parse("[+42.0106]-PEPM[+15.9949]IDE")
    print(f"\nProForma: {annot.serialize()}")
    print(f"MSFragger: {annot.to_msfragger()}\n")

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

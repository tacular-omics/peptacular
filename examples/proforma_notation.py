"""
ProForma Notation Examples
===========================
Comprehensive examples of supported ProForma 2.1 notation in peptacular.
Demonstrates parsing and serialization of various modification types and features.
"""

import peptacular as pt


def run():
    # ============================================================================
    # BASIC SEQUENCES
    # ============================================================================

    print("=" * 60)
    print("BASIC SEQUENCES")
    print("=" * 60)

    # Simple unmodified peptide
    simple = pt.parse("PEPTIDE")
    print(f"Simple sequence: {simple.serialize()}")

    # ============================================================================
    # TERMINAL MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("TERMINAL MODIFICATIONS")
    print("=" * 60)

    # Both terminals modified
    both = pt.parse("[Acetyl]-PEPTIDE-[Amidated]")
    print(f"Both terminals: {both.serialize()}")

    # Multiple N-terminal modifications
    multi_nterm = pt.parse("[Acetyl][Formyl]-PEPTIDE")
    print(f"Multiple N-term mods: {multi_nterm.serialize()}")

    # ============================================================================
    # INTERNAL MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("INTERNAL MODIFICATIONS")
    print("=" * 60)

    # Multiple different modifications
    multi_internal = pt.parse("PEM[Oxidation]TIS[Phospho]DE")
    print(f"Multiple modifications: {multi_internal.serialize()}")

    # Multiple modifications on same residue
    same_residue = pt.parse("PEM[Oxidation][Dioxidation]TIDE")
    print(f"Multiple mods on M: {same_residue.serialize()}")

    # ============================================================================
    # MODIFICATION NOTATION TYPES
    # ============================================================================

    print("\n" + "=" * 60)
    print("MODIFICATION NOTATION TYPES")
    print("=" * 60)

    # By name (Unimod/PSI-MOD)
    by_name = pt.parse("PEM[Oxidation]TIDE")
    print(f"By name: {by_name.serialize()}")

    # By accession number requires the UNIMOD: or MOD: prefix for Unimod/PSI-MOD respectively
    by_accession = pt.parse("PEM[UNIMOD:35]TIDE")
    print(f"By Unimod accession: {by_accession.serialize()}")

    # By mass (delta mass). requires sign (+/-)
    by_mass = pt.parse("PEM[+15.995]TIDE")
    print(f"By mass shift: {by_mass.serialize()}")
    neg_mass = pt.parse("PEPTIDE[-18.011]")
    print(f"Negative mass shift: {neg_mass.serialize()}")


    # By formula (requires Formula: prefix)
    by_formula = pt.parse("PEM[Formula:O]TIDE")
    print(f"By formula: {by_formula.serialize()}")

    # by glycan composition (requires Glycan: prefix)
    by_glycan = pt.parse("NEEYN[Glycan:Hex5HexNAc4]K")
    print(f"By glycan composition: {by_glycan.serialize()}")


    # ============================================================================
    # CHARGE STATES
    # ============================================================================

    print("\n" + "=" * 60)
    print("CHARGE STATES")
    print("=" * 60)

    # Positive charge
    charged_pos = pt.parse("PEPTIDE/2")
    print(f"Charge +2: {charged_pos.serialize()}")

    # Negative charge
    charged_neg = pt.parse("PEPTIDE/-2")
    print(f"Charge -2: {charged_neg.serialize()}")

    # ============================================================================
    # CHARGE ADDUCTS
    # ============================================================================

    print("\n" + "=" * 60)
    print("CHARGE ADDUCTS")
    print("=" * 60)

    # Single adduct (Total charge = +1)
    na_adduct = pt.parse("PEPTIDE/[Na:z+1]")
    print(f"Sodium adduct: {na_adduct.serialize()}")

    # Multiple copies of same adduct (Total charge = +2)
    multi_adduct = pt.parse("PEPTIDE/[Na:z+1^2]")
    print(f"Two sodium adducts: {multi_adduct.serialize()}")

    # Multiple different adducts (separated by commas) (Total charge = +3)
    mixed_adducts = pt.parse("PEPTIDE/[Na:z+1^2,H:z+1]")
    print(f"Mixed adducts: {mixed_adducts.serialize()}")

    # Metal adduct with charge (Total charge = +2)
    zn_adduct = pt.parse("PEPTIDE/[Zn:z+2]")
    print(f"Zinc adduct (+2): {zn_adduct.serialize()}")

    # ============================================================================
    # LABILE MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("LABILE MODIFICATIONS")
    print("=" * 60)

    labile = pt.parse("{Glycan:Hex}PEPTIDE")
    print(f"Labile glycan: {labile.serialize()}")

    multi_labile = pt.parse("{Phospho}PEPTIDE")
    print(f"Multiple labile: {multi_labile.serialize()}")

    # ============================================================================
    # GLYCAN NOTATION
    # ============================================================================

    print("\n" + "=" * 60)
    print("GLYCAN NOTATION")
    print("=" * 60)

    # Simple glycan
    simple_glycan = pt.parse("NEEYN[Glycan:Hex5HexNAc4]K")
    print(f"N-glycan: {simple_glycan.serialize()}")

    # ============================================================================
    # FIXED/STATIC MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("FIXED/STATIC MODIFICATIONS")
    print("=" * 60)

    # Fixed modification applied to all matching residues (M and T on all positions)
    fixed_mod = pt.parse("<[Oxidation]@M,T>MEMTIMDE")
    print(f"Fixed oxidation on all M and T: {fixed_mod.serialize()}")

    # Multiple fixed modifications
    multi_fixed = pt.parse("<[Oxidation]@M><[Phospho]@S>MSPETIDE")
    print(f"Multiple fixed mods: {multi_fixed.serialize()}")

    # Fixed modification with position rules (N-term Proline)
    fixed_nterm = pt.parse("<[Acetyl]@N-term:P>PEPTIDE")
    print(f"Fixed N-term mod: {fixed_nterm.serialize()}")

    # Fixed modification with position rules (Any C-term)
    fixed_cterm = pt.parse("<[Amidated]@C-term>PEPTIDE")
    print(f"Fixed C-term mod: {fixed_cterm.serialize()}")

    # ============================================================================
    # ISOTOPE LABELING
    # ============================================================================

    print("\n" + "=" * 60)
    print("ISOTOPE LABELING")
    print("=" * 60)

    # C13 labeling (all carbons)
    c13 = pt.parse("<13C>PEPTIDE")
    print(f"C13 labeled: {c13.serialize()}")

    # N15 labeling
    n15 = pt.parse("<15N>PEPTIDE")
    print(f"N15 labeled: {n15.serialize()}")

    # Multiple isotope labels
    multi_isotope = pt.parse("<13C><15N>PEPTIDE")
    print(f"C13 and N15 labeled: {multi_isotope.serialize()}")

    # Deuterium labeling
    deuterium = pt.parse("<2H>PEP[Oxidation]TIDE") 
    print(f"Deuterium labeled: {deuterium.serialize()}")

    # ============================================================================
    # AMBIGUOUS MODIFICATIONS (UNKNOWN POSITION)
    # ============================================================================

    print("\n" + "=" * 60)
    print("AMBIGUOUS MODIFICATIONS (UNKNOWN POSITION)")
    print("=" * 60)

    # Unknown position
    unknown_pos = pt.parse("[Phospho]?PEPTIDE")
    print(f"Phospho somewhere: {unknown_pos.serialize()}")

    # Multiple unknown modifications (Support caret for specifying multiple occurrences)
    multi_unknown = pt.parse("[Phospho]^2[Acetyl]?PEPTIDE")
    print(f"Multiple unknown: {multi_unknown.serialize()}")

    # ============================================================================
    # INTERVAL NOTATION (AMBIGUOUS LOCALIZATION)
    # ============================================================================

    print("\n" + "=" * 60)
    print("INTERVAL NOTATION (LOCALIZATION RANGES)")
    print("=" * 60)

    # Modification in a range (1-indexed, inclusive)
    interval = pt.parse("P(EP)[Phospho]TIDE")
    print(f"Phospho in positions 1-3: {interval.serialize()}")

    # Ambiguous interval (EP or PT or something with similar mass)
    ambiguous_interval = pt.parse("P(?EP)[Phospho]TIDE")
    print(f"Ambiguous intervals: {ambiguous_interval.serialize()}")


    # ============================================================================
    # INFO TAGS
    # ============================================================================

    print("\n" + "=" * 60)
    print("INFO TAGS (NON-MODIFICATION ANNOTATIONS)")
    print("=" * 60)

    # Info tag (no mass contribution)
    info_tag = pt.parse("PEPT[INFO:test]IDE")
    print(f"Info tag: {info_tag.serialize()}")


    # ============================================================================
    # PEPTIDE NAMING
    # ============================================================================

    print("\n" + "=" * 60)
    print("PEPTIDE NAMING")
    print("=" * 60)

    # Peptidoform name
    peptide_name = pt.parse("(>MyPeptide)PEPTIDE")
    print(f"Peptide name: {peptide_name.serialize()}")


    # ============================================================================
    # Multiple FEATURES COMBINED
    # ============================================================================


    print("\n" + "=" * 60)
    print("MULTIPLE FEATURES COMBINED")
    print("=" * 60)

    # Combined info tag and modification
    multi_info = pt.parse("PEPT[Phospho|INFO:quality=high]IDE")
    print(f"Info + modification: {multi_info.serialize()}")

    # Technically this is valid but no reason to do this. Peptacular only looks at the first modification in such cases.
    multi_annot2 = pt.parse("PEPT[Phospho|Oxidation|+76.0]IDE")
    print(f"Info + modification: {multi_annot2.serialize()}")


if __name__ == "__main__":
    run()
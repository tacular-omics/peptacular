"""Tests for ProFormaAnnotation.from_ms2_pip method"""

import unittest

from peptacular.annotation import ProFormaAnnotation


class TestFromMS2PIP(unittest.TestCase):
    """Test cases for MS2PIP format parsing"""

    def test_empty_modification_string(self):
        """Test parsing with no modifications"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="", static_mods=None)
        self.assertEqual(annot.serialize(), "PEPTIDE")

    def test_single_nterm_modification(self):
        """Test N-terminal modification (loc=0)"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="0|Acetyl")
        self.assertEqual(annot.serialize(), "[Acetyl]-PEPTIDE")

    def test_single_cterm_modification(self):
        """Test C-terminal modification (loc=-1)"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="-1|Amidated")
        self.assertEqual(annot.serialize(), "PEPTIDE-[Amidated]")

    def test_single_internal_modification(self):
        """Test internal modification (1-indexed)"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="2|Phospho")
        self.assertEqual(annot.serialize(), "PE[Phospho]PTIDE")

    def test_multiple_internal_modifications_same_position(self):
        """Test multiple modifications at same position"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="2|Phospho|2|Oxidation")
        serialized = annot.serialize()
        self.assertIn("Phospho", serialized)
        self.assertIn("Oxidation", serialized)

    def test_multiple_internal_modifications_different_positions(self):
        """Test modifications at different positions"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="2|Phospho|5|Oxidation")
        self.assertEqual(annot.serialize(), "PE[Phospho]PTI[Oxidation]DE")

    def test_mixed_modification_locations(self):
        """Test N-term, C-term, and internal modifications"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="0|Acetyl|3|Phospho|-1|Amidated")
        self.assertEqual(annot.serialize(), "[Acetyl]-PEP[Phospho]TIDE-[Amidated]")

    def test_static_modifications(self):
        """Test static and variable modifications"""
        annot = ProFormaAnnotation.from_ms2_pip(
            sequence="PEPTCIDE",
            mod_str="2|Phospho",
            static_mods={"C": "Carbamidomethyl"},
        )
        self.assertEqual(annot.serialize(), "<[Carbamidomethyl]@C>PE[Phospho]PTCIDE")

    def test_invalid_modification_format(self):
        """Test that invalid format raises ValueError"""
        with self.assertRaises(ValueError):
            ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="2|Phospho|3")

    def test_modification_out_of_range(self):
        """Test that out of range position raises ValueError"""
        with self.assertRaises(ValueError):
            ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="10|Phospho")
        with self.assertRaises(ValueError):
            ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="-2|Phospho")

    def test_mass_modification(self):
        """Test mass-based modification"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="2|+79.966331")
        self.assertEqual(annot.serialize(), "PE[+79.966331]PTIDE")

    def test_unimod_accession(self):
        """Test Unimod accession"""
        annot = ProFormaAnnotation.from_ms2_pip(sequence="PEPTIDE", mod_str="2|UNIMOD:21")
        self.assertEqual(annot.serialize(), "PE[UNIMOD:21]PTIDE")

    def test_complex_example(self):
        """Test complex real-world example"""
        annot = ProFormaAnnotation.from_ms2_pip(
            sequence="PEPTCMIDEK",
            mod_str="0|Acetyl|5|Oxidation|10|Phospho",
            static_mods={"C": 57.02146},
        )
        self.assertEqual(annot.serialize(), "<[+57.02146]@C>[Acetyl]-PEPTC[Oxidation]MIDEK[Phospho]")


class TestToMS2PIP(unittest.TestCase):
    """Test cases for MS2PIP format serialization"""

    def test_unmodified_peptide(self):
        """Test unmodified peptide"""
        seq, mod_str = ProFormaAnnotation(sequence="PEPTIDE").to_ms2_pip()
        self.assertEqual((seq, mod_str), ("PEPTIDE", ""))

    def test_single_nterm_modification(self):
        """Test N-terminal modification"""
        annot = ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"Acetyl": 1})
        seq, mod_str = annot.to_ms2_pip()
        self.assertEqual((seq, mod_str), ("PEPTIDE", "0|Acetyl"))

    def test_single_cterm_modification(self):
        """Test C-terminal modification"""
        annot = ProFormaAnnotation(sequence="PEPTIDE", cterm_mods={"Amidated": 1})
        seq, mod_str = annot.to_ms2_pip()
        self.assertEqual((seq, mod_str), ("PEPTIDE", "-1|Amidated"))

    def test_single_internal_modification(self):
        """Test internal modification"""
        annot = ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"Phospho": 1}})
        seq, mod_str = annot.to_ms2_pip()
        self.assertEqual((seq, mod_str), ("PEPTIDE", "2|Phospho"))

    def test_mixed_modifications(self):
        """Test mixed N-term, C-term, and internal modifications"""
        annot = ProFormaAnnotation(
            sequence="PEPTIDE",
            nterm_mods={"Acetyl": 1},
            cterm_mods={"Amidated": 1},
            internal_mods={1: {"Phospho": 1}},
        )
        seq, mod_str = annot.to_ms2_pip()
        self.assertEqual(seq, "PEPTIDE")
        parts = set(mod_str.split("|"))
        self.assertEqual(parts, {"0", "Acetyl", "-1", "Amidated", "2", "Phospho"})

    def test_multiple_internal_modifications(self):
        """Test multiple internal modifications"""
        annot = ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"Phospho": 1}, 4: {"Oxidation": 1}})
        seq, mod_str = annot.to_ms2_pip()
        self.assertEqual(seq, "PEPTIDE")
        parts = set(mod_str.split("|"))
        self.assertEqual(parts, {"2", "Phospho", "5", "Oxidation"})

    def test_mass_modification(self):
        """Test mass-based modification"""
        annot = ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"+79.966331": 1}})
        seq, mod_str = annot.to_ms2_pip()
        self.assertEqual((seq, mod_str), ("PEPTIDE", "2|+79.966331"))

    def test_unimod_accession(self):
        """Test Unimod accession"""
        annot = ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"UNIMOD:21": 1}})
        seq, mod_str = annot.to_ms2_pip()
        self.assertEqual((seq, mod_str), ("PEPTIDE", "2|UNIMOD:21"))

    def test_rejects_unsupported_features(self):
        """Test that unsupported features raise ValueError"""
        with self.assertRaises(ValueError):
            ProFormaAnnotation(sequence="PEPTIDE", labile_mods={"Hex": 1}).to_ms2_pip()
        with self.assertRaises(ValueError):
            ProFormaAnnotation(sequence="PEPTIDE", unknown_mods={"Phospho": 1}).to_ms2_pip()
        with self.assertRaises(ValueError):
            ProFormaAnnotation(sequence="PEPTIDE", charge=2).to_ms2_pip()

    def test_rejects_modification_multipliers(self):
        """Test that modification multipliers raise ValueError"""
        with self.assertRaises(ValueError):
            ProFormaAnnotation(sequence="PEPTIDE", nterm_mods={"Acetyl": 2}).to_ms2_pip()

    def test_rejects_multiple_mods_same_site(self):
        """Test that multiple mods at same site raise ValueError"""
        with self.assertRaises(ValueError):
            ProFormaAnnotation(sequence="PEPTIDE", internal_mods={1: {"Phospho": 1, "Oxidation": 1}}).to_ms2_pip()

    def test_roundtrip(self):
        """Test roundtrip: from_ms2_pip -> to_ms2_pip"""
        original_seq = "PEPTIDE"
        original_mod_str = "0|Acetyl|2|Phospho|7|Oxidation|-1|Amidated"

        annot = ProFormaAnnotation.from_ms2_pip(original_seq, original_mod_str)
        seq, mod_str = annot.to_ms2_pip()

        self.assertEqual(seq, original_seq)
        self.assertEqual(set(mod_str.split("|")), set(original_mod_str.split("|")))


if __name__ == "__main__":
    unittest.main()

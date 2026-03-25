import unittest

import peptacular as pt


class TestSlice(unittest.TestCase):
    def test_basic_slice(self):
        """Test basic slicing functionality"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")

        # Test normal slicing
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.sequence, "PEP")
        self.assertEqual(annotation.sequence, "PEPTIDE")  # Original unchanged

        # Test in-place slicing
        annotation.slice(0, 3, inplace=True)
        self.assertEqual(annotation.sequence, "PEP")

    def test_slice_preserves_nterm_modification(self):
        """Test slicing from start preserves N-terminal modification"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "[Acetyl]-PEP")

    def test_slice_preserves_cterm_modification(self):
        """Test slicing to end preserves C-terminal modification"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amidated]")
        sliced = annotation.slice(3, 7)
        self.assertEqual(sliced.serialize(), "TIDE-[Amidated]")

    def test_slice_removes_both_terminal_modifications_when_middle_slice(self):
        """Test middle slice removes both terminal modifications"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amidated]")
        sliced = annotation.slice(1, 5)
        self.assertEqual(sliced.serialize(), "EPTI")

    def test_slice_preserves_complete_interval(self):
        """Test slicing preserves interval when completely within slice"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        sliced = annotation.slice(0, 5)
        self.assertEqual(sliced.serialize(), "P(EP)[Phosphorylated]TI")

    def test_slice_preserves_complete_interval_full_sequence(self):
        """Test slicing preserves interval when taking full sequence"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        sliced = annotation.slice(0, 7)
        self.assertEqual(sliced.serialize(), "P(EP)[Phosphorylated]TIDE")

    def test_slice_raises_error_when_cutting_interval_start(self):
        """Test slicing raises ValueError when interval start gets cut"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(2, 5)

    def test_slice_raises_error_when_cutting_interval_end(self):
        """Test slicing raises ValueError when interval end gets cut"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(0, 2)

    def test_slice_raises_error_when_cutting_both_interval_ends(self):
        """Test slicing raises ValueError when both interval ends get cut"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(2, 4)

    def test_slice_preserves_complete_ambiguous_interval(self):
        """Test slicing preserves ambiguous interval when completely within slice"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPTI)[Phosphorylated]DE")
        self.assertEqual(annotation.serialize(), "P(?EPTI)[Phosphorylated]DE")

        sliced = annotation.slice(0, 7)
        self.assertEqual(sliced.serialize(), "P(?EPTI)[Phosphorylated]DE")

    def test_slice_raises_error_when_cutting_ambiguous_interval(self):
        """Test slicing raises ValueError when ambiguous interval gets cut"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPTI)[Phosphorylated]DE")

        with self.assertRaises(ValueError):
            annotation.slice(2, 4)

    def test_slice_preserves_labile_modifications(self):
        """Test slicing preserves labile modifications"""
        annotation = pt.ProFormaAnnotation.parse("{LabileMod}PEPTIDE")
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "{LabileMod}PEP")

    def test_slice_preserves_static_modifications(self):
        """Test slicing preserves static modifications"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"57@C": 1})
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "<57@C>PEP")

    def test_slice_preserves_charge_state(self):
        """Test slicing preserves charge state"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", charge=2)
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "PEP/2")

    def test_slice_drops_nterm_only_static_mod_when_nterm_removed(self):
        """N-term-only static mod is dropped when start > 0"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@N-term": 1})
        sliced = annotation.slice(1, 7)
        self.assertFalse(sliced.has_static_mods)

    def test_slice_drops_cterm_only_static_mod_when_cterm_removed(self):
        """C-term-only static mod is dropped when stop < seq_len"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@C-term": 1})
        sliced = annotation.slice(0, 6)
        self.assertFalse(sliced.has_static_mods)

    def test_slice_strips_nterm_rule_from_mixed_static_mod(self):
        """Mixed static mod keeps only non-N-term rules when N-term is removed"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@K,N-term": 1})
        sliced = annotation.slice(1, 7)
        self.assertTrue(sliced.has_static_mods)
        self.assertIn("[+1]@K", sliced.serialize())
        self.assertNotIn("N-term", sliced.serialize())

    def test_slice_preserves_anywhere_static_mod(self):
        """ANYWHERE-rule static mod is preserved through any slice"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@M": 1})
        sliced = annotation.slice(2, 5)
        self.assertTrue(sliced.has_static_mods)
        self.assertIn("[+1]@M", sliced.serialize())

    def test_slice_preserves_global_static_mod_no_rules(self):
        """Static mod with no position rules is preserved through any slice"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]": 1})
        sliced = annotation.slice(2, 5)
        self.assertTrue(sliced.has_static_mods)
        self.assertIn("[+1]", sliced.serialize())

    def test_slice_strips_both_terminal_rules_from_mixed_static_mod(self):
        """Static mod with K, N-term, and C-term rules keeps only K when both terminals removed"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@K,N-term,C-term": 1})
        sliced = annotation.slice(1, 6)
        self.assertTrue(sliced.has_static_mods)
        self.assertIn("[+1]@K", sliced.serialize())
        self.assertNotIn("N-term", sliced.serialize())
        self.assertNotIn("C-term", sliced.serialize())

    def test_slice_drops_static_mod_with_only_terminal_rules_when_both_removed(self):
        """Static mod with only N-term and C-term rules is dropped when both terminals removed"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@N-term,C-term": 1})
        sliced = annotation.slice(1, 6)
        self.assertFalse(sliced.has_static_mods)

    def test_slice_preserves_nterm_static_mod_when_nterm_kept(self):
        """N-term-only static mod is preserved when slice starts at position 0"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@N-term": 1})
        sliced = annotation.slice(0, 4)
        self.assertTrue(sliced.has_static_mods)
        self.assertIn("N-term", sliced.serialize())

    def test_slice_preserves_cterm_static_mod_when_cterm_kept(self):
        """C-term-only static mod is preserved when slice ends at seq_len"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", static_mods={"[+1]@C-term": 1})
        sliced = annotation.slice(3, 7)
        self.assertTrue(sliced.has_static_mods)
        self.assertIn("C-term", sliced.serialize())


if __name__ == "__main__":
    unittest.main()

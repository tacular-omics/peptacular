import unittest

import peptacular as pt


class TestShift(unittest.TestCase):
    def test_basic_shift(self):
        """Test basic shifting functionality"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")

        # Test positive shift
        shifted = annotation.shift(2)
        self.assertEqual(shifted.sequence, "PTIDEPE")
        self.assertEqual(annotation.sequence, "PEPTIDE")  # Original unchanged

        # Test in-place modification
        annotation.shift(2, inplace=True)
        self.assertEqual(annotation.sequence, "PTIDEPE")

    def test_zero_shift_returns_original(self):
        """Test zero shift returns original sequence"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = annotation.shift(0)
        self.assertEqual(result.sequence, "PEPTIDE")

    def test_full_length_shift_returns_original(self):
        """Test full length shift returns original sequence"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = annotation.shift(7)
        self.assertEqual(result.sequence, "PEPTIDE")

    def test_negative_shift(self):
        """Test negative shift moves sequence backwards"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = annotation.shift(-2)
        self.assertEqual(result.serialize(), "DEPEPTI")

    def test_large_shift_wraps_around(self):
        """Test large shift number wraps around correctly"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = annotation.shift(9)
        self.assertEqual(result.serialize(), "PTIDEPE")

    def test_empty_sequence_shift(self):
        """Test shifting empty sequence returns empty"""
        empty = pt.ProFormaAnnotation(sequence="")
        result = empty.shift(5)
        self.assertEqual(result.sequence, "")

    def test_single_amino_acid_shift(self):
        """Test shifting single amino acid returns same"""
        single = pt.ProFormaAnnotation(sequence="A")
        result = single.shift(5)
        self.assertEqual(result.sequence, "A")

    def test_shift_with_single_internal_modification(self):
        """Test shifting with single internal modification"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTIDE")
        result = annotation.shift(2)
        self.assertEqual(result.serialize(), "PTIDEPE[Phospho]")

    def test_shift_with_multiple_internal_modifications(self):
        """Test shifting with multiple internal modifications"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        result = annotation.shift(2)
        self.assertEqual(result.serialize(), "PTI[Methyl]DEPE[Phospho]")

    def test_shift_with_stacked_internal_modifications(self):
        """Test shifting with stacked internal modifications"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho][Methyl]PTIDE")
        result = annotation.shift(1)
        self.assertEqual(result.serialize(), "E[Phospho][Methyl]PTIDEP")

    def test_shift_preserves_nterm_modification(self):
        """Test shifting preserves N-terminal modification"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        result = annotation.shift(2)
        self.assertEqual(result.serialize(), "[Acetyl]-PTIDEPE")

    def test_shift_preserves_cterm_modification(self):
        """Test shifting preserves C-terminal modification"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        result = annotation.shift(3)
        self.assertEqual(result.serialize(), "TIDEPEP-[Amide]")

    def test_shift_preserves_both_terminal_modifications(self):
        """Test shifting preserves both terminal modifications"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        result = annotation.shift(1)
        self.assertEqual(result.serialize(), "[Acetyl]-EPTIDEP-[Amide]")

    def test_shift_preserves_labile_modification(self):
        """Test shifting preserves labile modification"""
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        result = annotation.shift(2)
        self.assertEqual(result.serialize(), "{Glycan}PTIDEPE")

    def test_shift_preserves_static_modification(self):
        """Test shifting preserves static modification"""
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        result = annotation.shift(2)
        self.assertEqual(result.serialize(), "<57@C>PTIDEPE")

    def test_shift_preserves_charge_state(self):
        """Test shifting preserves charge state"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        result = annotation.shift(2)
        self.assertEqual(result.serialize(), "PTIDEPE/2")

    def test_shift_preserves_unknown_modification(self):
        """Test shifting preserves unknown modification"""
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        result = annotation.shift(2)
        self.assertEqual(result.serialize(), "[Unknown]?PTIDEPE")

    def test_shift_raises_error_when_interval_wraps_around_case1(self):
        """Test shifting raises ValueError when interval wraps around - case 1"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        with self.assertRaises(ValueError):
            annotation.shift(2)

    def test_shift_raises_error_when_interval_wraps_around_case2(self):
        """Test shifting raises ValueError when interval wraps around - case 2"""
        annotation = pt.ProFormaAnnotation.parse("PE(?PTI)DE")
        with self.assertRaises(ValueError):
            annotation.shift(3)

    def test_shift_preserves_ambiguous_interval_when_no_wrap(self):
        """Test shifting preserves ambiguous interval when it doesn't wrap around"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPT)[Phospho]IDE")
        result = annotation.shift(1)
        self.assertEqual(result.serialize(), "(?EPT)[Phospho]IDEP")

    def test_shift_preserves_regular_interval_when_no_wrap(self):
        """Test shifting preserves regular interval when it doesn't wrap around"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        result = annotation.shift(1)
        self.assertEqual(result.serialize(), "(EP)[Phospho]TIDEP")

    def test_shift_keep_terminals_negative_and_large_shift(self):
        """Negative and >= seq_len shifts reduce modulo middle length, not full length"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE")
        # Middle is "PTIDE" (len 5): n=-1 shifts right by 1 -> "EPTID"
        self.assertEqual(annot.shift(-1, keep_nterm=2).serialize(), "PEEPTID")
        # Round-trip inverse with kept terminals
        self.assertEqual(annot.shift(1, keep_nterm=2).shift(-1, keep_nterm=2).serialize(), "PEPTIDE")
        # n == seq_len (7) still shifts the 5-residue middle by 7 % 5 == 2
        self.assertEqual(annot.shift(7, keep_nterm=2).serialize(), "PEIDEPT")
        # Both terminals kept: middle is "EPTID" (len 5), n=-1 -> "DEPTI"
        self.assertEqual(annot.shift(-1, keep_nterm=1, keep_cterm=1).serialize(), "PDEPTIE")
        # Internal modification tracks its residue under negative shift with keep_nterm
        mod_annot = pt.ProFormaAnnotation.parse("PEPT[Phospho]IDE")
        self.assertEqual(mod_annot.shift(-1, keep_nterm=2).serialize(), "PEEPT[Phospho]ID")

    def test_shift_interval_ending_at_boundary_does_not_wrap(self):
        """Intervals that land flush with the end of the shifted region must not raise"""
        annot = pt.ProFormaAnnotation.parse("A(BC)[+10]D")
        self.assertEqual(annot.shift(-1).serialize(), "DA(BC)[+10]")
        with self.assertRaises(ValueError):
            annot.shift(-2)

        kept = pt.ProFormaAnnotation.parse("PE(PT)[+10]IDE")
        self.assertEqual(kept.shift(-3, keep_nterm=2).serialize(), "PEIDE(PT)[+10]")

        full_interval = pt.ProFormaAnnotation.parse("(ABCD)[+10]")
        self.assertEqual(full_interval.shift(4).serialize(), "(ABCD)[+10]")
        with self.assertRaises(ValueError):
            full_interval.shift(1)


if __name__ == "__main__":
    unittest.main()

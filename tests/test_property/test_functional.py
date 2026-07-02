import unittest

import peptacular as pt

# manly for test coverage, as the test_mixin.py already tests the property methods via the ProFormaAnnotation.prop interface


class TestFunctionalProperties(unittest.TestCase):
    """Test the functional property calculations in peptacular"""

    def test_hydrophobicity_functional(self):
        # Single string
        val = pt.hydrophobicity("LIVM")
        self.assertIsInstance(val, float)

        # Single Annotation
        annot = pt.ProFormaAnnotation("LIVM")
        val2 = pt.hydrophobicity(annot)
        self.assertIsInstance(val2, float)
        self.assertEqual(val, val2)

        # List of strings
        vals = pt.hydrophobicity(["LIVM", "RKDE"])
        self.assertIsInstance(vals, list)
        self.assertEqual(len(vals), 2)
        self.assertEqual(vals[0], val)

    def test_pi_functional(self):
        val = pt.pi("DE")
        self.assertLess(val, 7.0)

        vals = pt.pi(["DE", "RK"])
        self.assertIsInstance(vals, list)
        self.assertLess(vals[0], 7.0)
        self.assertGreater(vals[1], 7.0)

    def test_pi_not_clamped(self):
        """pI must span the full pH range, not the old [4.05, 12.0] clamp.

        Regression: the bisection was bracketed to [4.05, 12.0], so strongly acidic
        peptides bottomed out at 4.05 and strongly basic ones topped out at 12.0.
        """
        self.assertLess(pt.pi("DDDDDDDDDD"), 4.05)
        self.assertGreater(pt.pi("RRRRRRRRRR"), 12.0)

    def test_pi_functional_matches_property(self):
        """The functional API and the annotation property must agree."""
        for seq in ("DDDDDDDDDD", "RRRRRRRRRR", "PEPTIDE"):
            self.assertAlmostEqual(pt.pi(seq), pt.parse(seq).prop.pi, places=6)

    def test_charge_at_ph_functional(self):
        # Single
        c = pt.charge_at_ph("K", pH=2.0)
        self.assertGreater(c, 0.0)

        # List
        cs = pt.charge_at_ph(["K", "D"], pH=7.0)
        self.assertIsInstance(cs, list)
        self.assertEqual(len(cs), 2)

    def test_aromaticity_functional(self):
        val = pt.aromaticity("WFY")
        self.assertIsInstance(val, float)

        vals = pt.aromaticity(["WFY", "AAA"])
        self.assertIsInstance(vals, list)
        self.assertGreater(vals[0], vals[1])

    def test_all_property_functions(self):
        """Verify all other property functions work on strings"""
        seq = "ACDEFGHIKLMNPQRSTVWY"

        # List of functions that return float for single sequence
        funcs = [
            pt.flexibility,
            pt.hydrophilicity,
            pt.surface_accessibility,
            pt.polarity,
            pt.mutability,
            pt.codons,
            pt.bulkiness,
            pt.recognition_factors,
            pt.transmembrane_tendency,
            pt.average_buried_area,
            pt.hplc,
            pt.refractivity,
        ]

        for func in funcs:
            val = func(seq)
            self.assertIsInstance(val, float, f"{func.__name__} failed to return float")

            vals = func([seq, seq])
            self.assertIsInstance(vals, list, f"{func.__name__} failed to return list")
            self.assertEqual(len(vals), 2)

    def test_secondary_structure_functional(self):
        # The main function returns dict
        res = pt.secondary_structure("AAAA")
        self.assertIsInstance(res, dict)
        self.assertIn(pt.SecondaryStructureType.ALPHA_HELIX, res)

        # Individual helpers
        self.assertIsInstance(pt.alpha_helix_percent("AAAA"), float)
        self.assertIsInstance(pt.beta_sheet_percent("AAAA"), float)
        self.assertIsInstance(pt.beta_turn_percent("AAAA"), float)
        self.assertIsInstance(pt.coil_percent("AAAA"), float)

    def test_calc_property_functional(self):
        seq = "A"
        scale = {"A": 10.0}
        val = pt.calc_property(seq, scale=scale)
        self.assertEqual(val, 10.0)

        vals = pt.calc_property([seq, seq], scale=scale)
        self.assertEqual(vals, [10.0, 10.0])

    def test_calc_window_property_functional(self):
        # NOTE: calc_window_property usually returns list[float] for a single sequence (windows)
        # It does NOT appear to support list of sequences based on the provided signature in attachment,
        # or maybe it does but signature provided was:
        # def calc_window_property(sequence: str | ProFormaAnnotation, ...) -> list[float]:

        seq = "AAAAA"
        scale = {"A": 1.0}
        windows = pt.calc_window_property(seq, scale=scale, window_size=3)
        self.assertIsInstance(windows, list)
        # Length depends on padding/implementation, just checking type and non-empty
        self.assertTrue(len(windows) > 0)
        self.assertIsInstance(windows[0], float)

    def test_aa_property_percentage_functional(self):
        seq = "ACDEF"
        residues = ["A", "C"]
        val = pt.aa_property_percentage(seq, residues=residues)
        self.assertIsInstance(val, float)
        # 2 out of 5 -> 0.4
        self.assertAlmostEqual(val, 0.4)


if __name__ == "__main__":
    unittest.main()

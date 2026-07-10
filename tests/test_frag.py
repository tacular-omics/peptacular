import unittest

import pytest
from tacular import IonType

import peptacular as pt
from peptacular.annotation.frag import Fragment


class TestFragmentClass(unittest.TestCase):
    def test_init_simple(self):
        f = Fragment(ion_type=IonType.Y, position=2, mass=100.0, monoisotopic=True, charge_state=1, parent_sequence="AB")
        self.assertEqual(f.ion_type, IonType.Y)
        self.assertEqual(f.position, 2)
        self.assertEqual(f.mass, 100.0)
        self.assertEqual(f.charge_state, 1)
        self.assertEqual(f.parent_sequence, "AB")
        self.assertTrue(f.monoisotopic)

        # mz = mass / charge
        self.assertEqual(f.mz, 100.0)
        self.assertTrue(f.is_protonated)

        # By default (charge_state=1, no explicit adducts), it assumes protonation
        # neutral mass = mass - proton_mass
        # We can check it's less than mass
        self.assertLess(f.neutral_mass, 100.0)

    def test_charge_adducts(self):
        # Na+ adduct, charge 1
        # Need to format adduct dictionary correctly: string -> int count
        f = Fragment(ion_type=IonType.B, position=3, mass=200.0, monoisotopic=True, charge_state=1, charge_adducts=("Na:z+1",))
        self.assertEqual(f.mz, 200.0)
        self.assertFalse(f.is_protonated)

        adducts = f.charge_adducts
        self.assertEqual(len(adducts), 1)
        # Check properties of the adduct object returned if possible, or just string conversion

    def test_duplicate_charge_adducts_not_collapsed(self):
        # Two identical adducts must both be counted (regression: {k: 1} dict deduped them,
        # so neutral_mass subtracted only one carrier's mass and composition undercounted).
        one = Fragment(ion_type=IonType.B, position=3, mass=200.0, monoisotopic=True, charge_state=1, charge_adducts=("Na:z+1",))
        two = Fragment(ion_type=IonType.B, position=3, mass=200.0, monoisotopic=True, charge_state=2, charge_adducts=("Na:z+1", "Na:z+1"))
        self.assertEqual(two.charge_adducts.get_charge(), 2)
        self.assertEqual(dict(two.charge_adducts._mods), {"Na:z+1": 2})
        # a second sodium removes ~23 Da more when deriving the neutral mass
        self.assertLess(two.neutral_mass, one.neutral_mass)

    def test_isotopes_int(self):
        f = Fragment(
            ion_type=IonType.Y,
            position=1,
            mass=150.0,
            monoisotopic=True,
            charge_state=1,
            isotopes=2,  # 2 C13 atoms
        )
        self.assertTrue(f.is_c13)
        isos = f.isotopes
        # key should be C13 element info, value 2
        self.assertEqual(len(isos), 1)
        # Check value is 2
        val = next(iter(isos.values()))
        self.assertEqual(val, 2)

    def test_isotopes_dict(self):
        f = Fragment(ion_type=IonType.Y, position=1, mass=150.0, monoisotopic=True, charge_state=1, isotopes={"13C": 1, "15N": 1})
        # is_c13 returns True if isotopes is int, else False (based on implementation reading)
        self.assertFalse(f.is_c13)
        isos = f.isotopes
        self.assertEqual(len(isos), 2)

    def test_losses(self):
        # Loss of water
        f = Fragment(ion_type=IonType.Y, position=1, mass=150.0, monoisotopic=True, charge_state=1, deltas={"H2O": 1})
        losses = f.losses
        self.assertEqual(len(losses), 1)
        # Key should be ChargedFormula for H2O

        # Loss by mass
        f2 = Fragment(ion_type=IonType.Y, position=1, mass=150.0, monoisotopic=True, charge_state=1, deltas={18.01: 1})
        losses2 = f2.losses
        self.assertEqual(losses2[18.01], 1)

    def test_position(self):
        annot = pt.parse("PEPTIDE")

        frag = annot.frag(position=3, ion_type=IonType.B)
        self.assertEqual(frag.position, 3)
        self.assertEqual(frag.parent_sequence, "PEPTIDE")
        self.assertEqual(frag.sequence, "PEP")

    def test_frag_charge(self):
        annot = pt.parse("PEPTIDE/2")

        frag = annot.frag(position=3, ion_type=IonType.IMMONIUM)
        self.assertEqual(frag.position, 3)
        self.assertEqual(frag.parent_sequence, "PEPTIDE/2")
        self.assertEqual(frag.charge_state, 2)
        self.assertEqual(frag.sequence, "P/2")

    def test_frag_set_charge(self):
        annot = pt.parse("PEPTIDE/2")

        frag = annot.frag(position=(2, 5), charge=3, ion_type=IonType.BY)
        self.assertEqual(frag.position, (2, 5))
        self.assertEqual(frag.parent_sequence, "PEPTIDE/3")
        self.assertEqual(frag.charge_state, 3)
        self.assertEqual(frag.sequence, "EPTI/3")

    def test_frag_set_charge_adduct(self):
        annot = pt.parse("PEPTIDE/2")

        frag = annot.frag(position=3, charge="Na:z+1^2", ion_type=IonType.Z)
        self.assertEqual(frag.position, 3)
        self.assertEqual(frag.parent_sequence, "PEPTIDE/[Na:z+1^2]")
        self.assertEqual(frag.charge_state, 2)

        frag = annot.frag(position=3, charge="Na:z+1", ion_type=IonType.Z)
        self.assertEqual(frag.position, 3)
        self.assertEqual(frag.parent_sequence, "PEPTIDE/[Na:z+1]")
        self.assertEqual(frag.charge_state, 1)
        self.assertEqual(frag._charge_adducts, ("Na:z+1",))
        self.assertEqual(frag.charge_adducts.serialize(), "[Na:z+1]")

    def test_frag_set_bad_charge_adduct(self):
        annot = pt.parse("PEPTIDE/2")

        with pytest.raises(ValueError):
            frag = annot.frag(position=3, charge="Na", ion_type=IonType.Z)
            self.assertEqual(frag.position, 3)
            self.assertEqual(frag.parent_sequence, "PEPTIDE/[Na]")

    def test_frag_position(self):
        annot = pt.parse("PEPTIDE")

        frag = annot.frag(position=4, ion_type=IonType.Y)
        self.assertEqual(frag.position, 4)
        self.assertEqual(frag.parent_sequence, "PEPTIDE")
        self.assertEqual(frag.sequence, "TIDE")

    def test_frag_position_none(self):
        annot = pt.parse("PEPTIDE")

        frag = annot.frag(position=None, ion_type=IonType.Y)
        self.assertEqual(frag.position, 7)
        self.assertEqual(frag.parent_sequence, "PEPTIDE")
        self.assertEqual(frag.sequence, "PEPTIDE")


if __name__ == "__main__":
    unittest.main()

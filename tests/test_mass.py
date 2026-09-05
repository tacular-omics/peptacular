import pytest

import peptacular as pt

PEPTIDE_Y1_MONO_MASS = 800.36728
PEPTIDE_Y1_AVG_MASS = 800.831


class TestApplyMods:
    def test_basic_mass(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(781.349)

    def test_basic_average_mass(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass = seq._base_mass(monoisotopic=False)[0]
        assert mass == pytest.approx(781.808)

    def test_mass_with_internal_mod(self):
        seq = pt.ProFormaAnnotation.parse("PEP[+10]TIDE")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(791.349)

    def test_mass_with_nterm_mod(self):
        seq = pt.ProFormaAnnotation.parse("[+10]-PEPTIDE")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(791.349)

    def test_mass_with_cterm_mod(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE-[+10]")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(791.349)

    def test_mass_with_multiple_mods(self):
        seq = pt.ProFormaAnnotation.parse("[+5]-PEP[+10]TIDE-[-15]")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(781.349)

    def test_mass_with_labile_mod(self):
        seq = pt.ProFormaAnnotation.parse("{+10}PEPTIDE")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(791.349)

    def test_mass_with_unknown_mod(self):
        seq = pt.ProFormaAnnotation.parse("[+10]?PEPTIDE")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(791.349)

    def test_mass_with_interval_mod(self):
        seq = pt.ProFormaAnnotation.parse("PEP(TID)[+10]E")
        mass = seq._base_mass(monoisotopic=True)[0]
        assert mass == pytest.approx(791.349)

    def test_mass_ave_and_mono(self):
        sequence = "PEPTIDE"

        mass_mono = pt.mass(sequence, charge=1, ion_type="y", monoisotopic=True)
        mass_avg = pt.mass(sequence, charge=1, ion_type="y", monoisotopic=False)
        assert mass_mono == pytest.approx(PEPTIDE_Y1_MONO_MASS)
        assert mass_avg == pytest.approx(PEPTIDE_Y1_AVG_MASS)

    def test_mass_with_deltas(self):
        sequence = "PEPTIDE"

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            deltas=10.0,
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS + 10)

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            deltas=-10.0,
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS - 10)

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            deltas="C2",
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS + 24)

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            deltas={"C": 2},
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS + 24)

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            deltas={"C": -2},
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS - 24)

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            deltas={"C": 2, 10: 1},
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS + 34)

        # ensure raises value error for charged formula in deltas
        with pytest.raises(ValueError):
            mass = pt.mass(
                sequence,
                charge=1,
                ion_type="y",
                monoisotopic=True,
                deltas={"C:z+1": 2, 10: 1},
            )

    def test_mass_with_isotopes(self):
        sequence = "PEPTIDE"

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            isotopes=1,
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS + pt.C13_NEUTRON_MASS)

        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            isotopes=2,
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS + 2 * pt.C13_NEUTRON_MASS)

        # assert throws error
        with pytest.raises(ValueError):
            mass = pt.mass(
                sequence,
                charge=1,
                ion_type="y",
                isotopes=-2,
            )

        # Custom isotopic composition
        mass = pt.mass(
            sequence,
            charge=1,
            ion_type="y",
            isotopes={"13C": 2},
        )
        assert mass == pytest.approx(PEPTIDE_Y1_MONO_MASS + 2 * pt.C13_NEUTRON_MASS)

        # ensure raises error with unknown isotope
        with pytest.raises(KeyError):
            mass = pt.mass(
                sequence,
                charge=1,
                ion_type="y",
                isotopes={"15C": 2},
            )

        with pytest.raises(ValueError, match="negative element counts"):
            mass = pt.mass(
                sequence,
                charge=1,
                ion_type="y",
                isotopes={"13C": 200},
            )

        # will throw error since we check composition directly
        with pytest.raises(ValueError):
            mass = pt.mass(
                sequence,
                charge=1,
                ion_type="y",
                isotopes={"13C": 2000},
                calculate_with_composition=True,
            )


if __name__ == "__main__":
    pytest.main([__file__])

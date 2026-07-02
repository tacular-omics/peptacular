"""Regression tests for mass/composition correctness bugs."""

import peptacular as pt


def _comp_mass(seq: str) -> float:
    a = pt.parse(seq)
    return sum(e.get_mass(monoisotopic=True) * n for e, n in a.comp().items())


class TestGlobalIsotopeMass:
    """Regression: mass()/mz() must apply global isotope labels (<13C>, <15N>)."""

    def test_13c_label_shifts_mass(self):
        base = pt.mass("PEPTIDE")
        labeled = pt.mass("<13C>PEPTIDE")
        # 34 carbons * (13.00335 - 12.0) ~= 34.1 Da
        assert labeled - base > 33.0
        assert labeled - base < 35.0

    def test_15n_label_shifts_mass(self):
        base = pt.mass("PEPTIDE")
        labeled = pt.mass("<15N>PEPTIDE")
        assert labeled > base

    def test_mass_matches_composition_path(self):
        for seq in ["<13C>PEPTIDE", "<15N>PEPTIDE", "<13C><15N>PEPTIDE", "<13C>PEM[Oxidation]TIDEK"]:
            fast = pt.parse(seq).mass(charge=0)
            comp = pt.parse(seq).mass(charge=0, calculate_with_composition=True)
            assert abs(fast - comp) < 1e-6, f"{seq}: fast={fast} comp={comp}"
            assert abs(fast - _comp_mass(seq)) < 1e-6


class TestNegativeCompositionDelta:
    """Regression: comp() must keep atom-removing modifications (e.g. Amidated's O-1)."""

    def test_amidated_composition_removes_oxygen(self):
        base = pt.parse("PEPTIDEK").comp()
        amid = pt.parse("PEPTIDEK-[Amidated]").comp()
        # find oxygen ElementInfo key
        o_delta = sum(n for e, n in amid.items() if str(e) == "O") - sum(n for e, n in base.items() if str(e) == "O")
        assert o_delta == -1

    def test_comp_mass_matches_mass_for_amidated(self):
        for seq in ["PEPTIDEK-[Amidated]", "PEPTIDE-[Amidated]", "[Acetyl]-PEPTIDEK-[Amidated]"]:
            assert abs(_comp_mass(seq) - pt.parse(seq).mass(charge=0)) < 1e-6, seq


class TestCondenseAmbiguityMass:
    """Regression: condensing ambiguity to X-notation must preserve neutral mass."""

    def test_condense_preserves_mass(self):
        for seq in ["PEPT(?ID)E", "PE(?PT)IDE", "P(?EPTIDE)", "PEPT(?ID)[Phospho]E"]:
            a = pt.parse(seq)
            before = a.mass(charge=0)
            after = a.condense_ambiguity_to_xnotation().mass(charge=0)
            assert abs(before - after) < 1e-4, f"{seq}: {before} -> {after}"

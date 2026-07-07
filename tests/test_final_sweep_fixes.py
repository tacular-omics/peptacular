"""Regression tests for bugs found in the final pre-publish verification sweep."""

import math

import peptacular as pt


class TestRepeatedModificationComposition:
    """comp() must scale a repeated modification's composition by its count, matching mass()."""

    def test_single_vs_repeated_same_mod(self):
        base = pt.parse("PEPTFIDE").mass()
        for n in (1, 2, 3):
            s = "PEPTF" + "[UNIMOD:310]" * n + "IDE"
            a = pt.parse(s)
            m = a.mass()
            cm = sum(e.get_mass() * c for e, c in a.comp().items())
            assert abs(m - cm) < 1e-4, f"n={n}: mass={m} comp={cm}"

    def test_mass_scales_linearly_with_count(self):
        m1 = pt.mass("PEPTF[UNIMOD:310]IDE")
        m2 = pt.mass("PEPTF[UNIMOD:310][UNIMOD:310]IDE")
        single_mod_mass = m1 - pt.mass("PEPTFIDE")
        assert abs((m2 - m1) - single_mod_mass) < 1e-6


class TestFunctionalFragmentChargeDefault:
    """The functional fragment() API must match the OOP .fragment() smart charge default."""

    def test_matches_oop_for_positive_precursor(self):
        oop = pt.parse("PEPTIDE/3").fragment(ion_types=["b", "y"])
        func = pt.fragment("PEPTIDE/3", ion_types=["b", "y"])
        assert len(oop) == len(func)
        assert sorted(f.charge_state for f in oop) == sorted(f.charge_state for f in func)

    def test_matches_oop_for_negative_precursor(self):
        oop_charges = sorted(set(f.charge_state for f in pt.parse("PEPTIDE/-3").fragment(ion_types=["b"])))
        func_charges = sorted(set(f.charge_state for f in pt.fragment("PEPTIDE/-3", ion_types=["b"])))
        assert oop_charges == func_charges
        assert all(c < 0 for c in func_charges)  # regression: used to silently return positive charge=1

    def test_explicit_charges_still_respected(self):
        r = pt.fragment("PEPTIDE", ion_types=["b"], charges=[1, 2])
        assert sorted(set(f.charge_state for f in r)) == [1, 2]


class TestOverlappingSubsequenceMatches:
    """Subsequence search must find overlapping occurrences, not just non-overlapping ones."""

    def test_find_subsequence_indices_overlapping(self):
        assert pt.find_subsequence_indices("IIIII", "II") == [0, 1, 2, 3]

    def test_coverage_overlapping(self):
        assert pt.coverage("AAAAA", ["AAA"]) == [1, 1, 1, 1, 1]

    def test_non_overlapping_case_unaffected(self):
        assert pt.find_subsequence_indices("PEPTIDEPEPTIDE", "PEPTIDE") == [0, 7]


class TestCountResiduesDoesNotMutate:
    """count_residues/percent_residues must not mutate a caller-supplied ProFormaAnnotation."""

    def test_count_residues_no_mutation(self):
        ann = pt.parse("PEPCTIDE")
        ann.static_mods = {"[Carbamidomethyl]@C": 1}
        before = ann.serialize()
        pt.count_residues(ann)
        assert ann.serialize() == before

    def test_percent_residues_no_mutation(self):
        ann = pt.parse("PEPCTIDE")
        ann.static_mods = {"[Carbamidomethyl]@C": 1}
        before = ann.serialize()
        pt.percent_residues(ann)
        assert ann.serialize() == before


class TestIP2ConversionAdjacentBrackets:
    """convert_ip2_sequence must handle adjacent modification brackets in every position."""

    def test_nterm_modification(self):
        assert pt.convert_ip2_sequence("K.(-1)PEP(phospho)TIDE.K") == "[-1]-PEP[phospho]TIDE"

    def test_multiple_cterm_modifications(self):
        assert pt.convert_ip2_sequence("K.PEPTIDE(2)(3).K") == "PEPTIDE[2]-[3]"

    def test_complex_modifications(self):
        assert pt.convert_ip2_sequence("-.(1)PEP(phospho)TIDE(2)(3).-") == "[1]-PEP[phospho]TIDE[2]-[3]"

    def test_internal_adjacent_mods_no_longer_crash(self):
        # Regression: two adjacent parenthetical mods sandwiched between residues
        # previously produced an unparseable 'PEP[mod1]-[mod2]TIDE'.
        result = pt.convert_ip2_sequence("K.PEP(mod1)(mod2)TIDE.K")
        assert result == "PEP[mod1][mod2]TIDE"
        pt.parse(result)  # must be parseable

    def test_nterm_adjacent_mods_no_longer_crash(self):
        result = pt.convert_ip2_sequence("K.(mod1)(mod2)PEPTIDE.K")
        assert result == "[mod1][mod2]-PEPTIDE"
        pt.parse(result)


class TestRandomizerIntervalCoverage:
    """Randomly generated intervals must be able to reach the sequence's final residue."""

    def test_intervals_can_reach_last_residue(self):
        import random as _random

        _random.seed(42)
        found = False
        for _ in range(500):
            a = pt.ProFormaAnnotation.random()
            if a.has_intervals and any(iv.end == len(a) for iv in a.intervals):
                found = True
                break
        assert found, "no randomly generated interval reached the last residue in 500 tries"


class TestGeneratePartitionsOverlap:
    """generate_partitions must respect aa_overlap=0 (no silent window overlap) where possible."""

    def test_default_overlap_zero_minimizes_overlap(self):
        from peptacular.property.core import generate_partitions
        from peptacular.property.data import HydrophobicityScale

        import peptacular.property.core as core_mod

        seen = []
        orig = core_mod.calc_property

        def spy(sequence, *a, **kw):
            seen.append(sequence)
            return orig(sequence, *a, **kw)

        core_mod.calc_property = spy
        try:
            seq = "ACDEFGHIKLMNPQRST"  # 17 residues, not divisible by 5
            generate_partitions(seq, HydrophobicityScale.KYTE_DOOLITTLE, num_windows=5, aa_overlap=0)
        finally:
            core_mod.calc_property = orig

        # Reconstruct window positions and count overlapping adjacent pairs.
        positions = []
        cursor = 0
        for w in seen:
            start = seq.find(w, max(0, cursor - len(w)))
            positions.append((start, start + len(w)))
            cursor = start + len(w)
        overlapping_pairs = sum(1 for i in range(len(positions) - 1) if positions[i][1] > positions[i + 1][0])
        # Regression: previously EVERY adjacent pair silently overlapped despite
        # aa_overlap=0 (3 of 4 pairs). Now only the unavoidable final-window clamp
        # (pulled left to end exactly at the sequence) may overlap with its neighbor.
        assert overlapping_pairs <= 1, f"positions={positions} overlapping_pairs={overlapping_pairs}"

    def test_requested_overlap_is_honored_for_most_pairs(self):
        from peptacular.property.core import generate_partitions
        from peptacular.property.data import HydrophobicityScale

        # Should not raise, and should return the requested number of windows.
        result = generate_partitions("ACDEFGHIKL", HydrophobicityScale.KYTE_DOOLITTLE, num_windows=3, aa_overlap=2)
        assert len(result) == 3
        assert all(isinstance(v, float) and not math.isnan(v) for v in result)

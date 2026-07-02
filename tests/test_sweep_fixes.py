"""Regression tests for bugs found in the full-sweep review."""

import peptacular as pt


class TestModifyEnumeration:
    """pt.modify must enumerate positional isomers by default (was collapsing them)."""

    def test_enumerates_all_placements(self):
        # 1 phospho on the 5 S/T/Y sites of SPTYSGTK -> empty + one per site = 6
        r = pt.modify("SPTYSGTK", internal_variable={"STY": [79.966]}, max_variable_mods=1)
        assert len(r) == 6

    def test_unique_peptidoforms_collapses(self):
        r = pt.modify("SPTYSGTK", internal_variable={"STY": [79.966]}, max_variable_mods=1, unique_peptidoforms=True)
        assert len(r) == 2  # empty + one representative

    def test_batch_parity(self):
        rb = pt.modify(["SPTYSGTK", "SPTYSGTK"], internal_variable={"STY": [79.966]}, max_variable_mods=1)
        assert [len(x) for x in rb] == [6, 6]

    def test_string_mod_value_is_single_mod(self):
        # A bare string mod must be applied whole, not shredded into per-character mods.
        r = pt.modify("SPTIDE", internal_variable={"S": "Phospho"}, max_variable_mods=1)
        assert set(r) == {"SPTIDE", "S[Phospho]PTIDE"}


class TestCondenseModsInplaceContract:
    """condense_mods_to_intervals(inplace=False) must not mutate the original."""

    def test_original_unchanged(self):
        a = pt.parse("PE(P[Acetyl]T)[Phospho]IDE")
        before = a.serialize()
        result = a.condense_mods_to_intervals(inplace=False)
        assert a.serialize() == before  # original untouched
        assert result.serialize() != before  # a real change was returned

    def test_copy_intervals_are_independent(self):
        a = pt.parse("PE(P[Acetyl]T)[Phospho]IDE")
        c = a.copy()
        assert c.intervals[0] is not a.intervals[0]


class TestNonSpecificDigestion:
    """Non-specific digestion must include the full-length peptide and single residues."""

    def test_includes_full_length(self):
        a = pt.parse("PEPTIDEK")
        spans = list(a.nonspecific_spans(min_len=1, max_len=8))
        assert len(spans) == 36  # sum_{L=1..8}(8-L+1)
        assert any(sp.end - sp.start == 8 for sp in spans)

    def test_single_residue_yields_itself(self):
        a = pt.parse("K")
        peps = [a[sp].serialize() for sp in a.nonspecific_spans()]
        assert peps == ["K"]

    def test_consistent_with_enzymatic_full_length(self):
        a = pt.parse("PEPTIDEK")
        non = {a[sp].serialize() for sp in a.nonspecific_spans(min_len=1, max_len=100)}
        enz = {a[sp].serialize() for sp in a.simple_digest(cleave_on="K", missed_cleavages=10, min_len=1, max_len=100)}
        assert enz <= non  # every enzymatic product (incl. full length) is a non-specific product


class TestSemiEnzymaticMaxLen:
    """Semi-enzymatic digestion must not drop in-range children of over-length parents."""

    def test_in_range_semi_peptides_present(self):
        a = pt.parse("AAKAAA")
        p5 = {a[sp].serialize() for sp in a.simple_digest(cleave_on="K", missed_cleavages=1, semi=True, min_len=1, max_len=5)}
        # These are all length <= 5 and valid semi-tryptic peptides of the mc=1 parent AAKAAA
        for expected in ("AKAAA", "AAKAA", "AAKA", "KAAA"):
            assert expected in p5

    def test_all_within_length_bounds(self):
        a = pt.parse("AAKAAA")
        for sp in a.simple_digest(cleave_on="K", missed_cleavages=1, semi=True, min_len=2, max_len=4):
            assert 2 <= (sp.end - sp.start) <= 4

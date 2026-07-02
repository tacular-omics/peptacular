"""Tests for enzyme cleavage rules (cleave_on / restrict / cterminal / ambiguity codes)."""

import peptacular as pt
from peptacular.digestion.core import _convert_to_aa_set, generate_regex


def _peptides(sequence: str, **kwargs) -> list[str]:
    ann = pt.parse(sequence)
    kwargs.setdefault("min_len", 1)
    kwargs.setdefault("max_len", 100)
    return [ann[span].serialize() for span in ann.simple_digest(**kwargs)]


class TestCleaveRules:
    def test_cterminal_cleavage(self):
        # Cleave C-terminal to K or R (trypsin-like)
        assert _peptides("AKBKPCKDR", cleave_on="KR") == ["AK", "BK", "PCK", "DR"]

    def test_restrict_after_blocks_cleavage(self):
        # Do not cleave when the residue after the cut site is P
        assert _peptides("AKBKPCKDR", cleave_on="KR", restrict_after="P") == ["AK", "BKPCK", "DR"]

    def test_restrict_before_blocks_cleavage(self):
        # Do not cleave after K when it is preceded by P
        assert _peptides("PKAKG", cleave_on="K", restrict_before="P") == ["PKAK", "G"]

    def test_nterminal_cleavage(self):
        # cterminal=False cleaves N-terminal to the residue (Asp-N style)
        assert _peptides("DKAEKG", cleave_on="K", cterminal=False) == ["D", "KAE", "KG"]


class TestAmbiguityCodes:
    def test_b_expands_to_d_and_n(self):
        assert _peptides("ADANG", cleave_on="B") == ["AD", "AN", "G"]

    def test_convert_aa_set_expansions(self):
        assert _convert_to_aa_set("B") == {"D", "N"}
        assert _convert_to_aa_set("J") == {"I", "L"}
        assert _convert_to_aa_set("Z") == {"E", "Q"}
        assert _convert_to_aa_set("KR") == {"K", "R"}

    def test_convert_aa_set_subtraction(self):
        # X is all residues; X-KR removes K and R
        full = _convert_to_aa_set("X")
        assert _convert_to_aa_set("X-KR") == full - {"K", "R"}

    def test_convert_aa_set_none_is_empty(self):
        assert _convert_to_aa_set(None) == set()


class TestMissedCleavages:
    def test_count_monotonic_in_missed_cleavages(self):
        ann = pt.parse("AKBKCKDK")
        counts = [len(list(ann.simple_digest(cleave_on="K", missed_cleavages=mc, min_len=1, max_len=100))) for mc in range(3)]
        assert counts == sorted(counts)
        assert counts[0] < counts[-1]

    def test_missed_cleavage_peptides_are_superset(self):
        ann = pt.parse("AKBKCK")
        p0 = set(_peptides("AKBKCK", cleave_on="K", missed_cleavages=0))
        p1 = set(_peptides("AKBKCK", cleave_on="K", missed_cleavages=1))
        assert p0.issubset(p1)


class TestGenerateRegex:
    def test_none_and_empty_are_nonspecific(self):
        # Empty pattern -> non-specific (regression guard for the earlier off-by-one fix)
        for cleave_on in (None, ""):
            assert generate_regex(cleave_on=cleave_on).pattern == ""

    def test_cterminal_uses_lookbehind(self):
        rgx = generate_regex(cleave_on="KR", cterminal=True)
        assert "(?<=" in rgx.pattern

    def test_nterminal_uses_lookahead(self):
        rgx = generate_regex(cleave_on="KR", cterminal=False)
        assert "(?=" in rgx.pattern

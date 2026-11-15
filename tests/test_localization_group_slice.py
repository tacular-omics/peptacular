"""A slice or digest must not separate a localization group's label from its primary tag (#85)."""

import pytest

import peptacular as pt
from peptacular.diagnostics import PeptacularError


def test_digest_splitting_localization_group_raises():
    with pytest.raises(PeptacularError, match="#g1"):
        pt.digest("EM[#g1(0.10)]EVTKS[Phospho#g1(0.90)]ESPEK", enzyme="trypsin", missed_cleavages=0)


def test_digest_splitting_localization_group_primary_in_first_half_raises():
    with pytest.raises(PeptacularError, match="#g1"):
        pt.digest("EM[Phospho#g1(0.90)]EVTKS[#g1(0.10)]ESPEK", enzyme="trypsin", missed_cleavages=0)


def test_slice_splitting_localization_group_raises():
    with pytest.raises(PeptacularError, match="#g1"):
        pt.parse("EM[#g1(0.10)]EVTKS[Phospho#g1(0.90)]ESPEK")[0:6]


def test_group_wholly_inside_one_product_digests_and_keeps_score():
    products = pt.digest("PEPS[Phospho#g1(0.90)]T[#g1(0.10)]IDEK", enzyme="trypsin", missed_cleavages=0)
    serialized = [seq for seq, _ in products]
    assert any("Phospho#g1(0.90)" in seq and "#g1(0.10)" in seq for seq in serialized)


class TestSplitKeepsLocalizationGroups:
    """Splitting is not slicing: it keeps every position, so the guard above must not fire.

    `split()` returns one piece per residue and keeps all of them, so every member of a
    localization group is still present across the pieces. Nothing is dropped and no reference is
    orphaned, which is the loss #85 is about. Before this change the per-residue slices tripped
    the guard, so `pt.split`, `pt.count_residues`, `pt.percent_residues`, `pt.permutations`,
    `pt.combinations`, `pt.combinations_with_replacement`, `pt.product` and `pt.sort` all raised
    `PeptacularError` on any peptide whose group spans more than one residue.
    """

    SEQ = "PEPT[Phospho#g1]IDEK[#g1]"
    SCORED = "EM[#g1(0.10)]EVTKS[Phospho#g1(0.90)]ESPEK"

    def test_split_keeps_the_primary_tag_and_the_reference(self):
        assert pt.split(self.SEQ) == ["P", "E", "P", "T[Phospho#g1]", "I", "D", "E", "K[#g1]"]

    def test_every_piece_parses_on_its_own(self):
        assert all(pt.parse(piece).serialize() == piece for piece in pt.split(self.SEQ))

    def test_split_round_trips_through_join_with_scores(self):
        assert pt.join(pt.split(self.SCORED)) == self.SCORED

    def test_count_residues_counts_a_group_that_spans_two_residues(self):
        assert pt.count_residues(self.SEQ) == {"P": 2, "E": 2, "T[Phospho#g1]": 1, "I": 1, "D": 1, "K[#g1]": 1}

    @pytest.mark.parametrize(
        ("name", "call"),
        [
            ("percent_residues", lambda s: pt.percent_residues(s)),
            ("permutations", lambda s: pt.permutations(s, size=2)),
            ("combinations", lambda s: pt.combinations(s, size=2)),
            ("combinations_with_replacement", lambda s: pt.combinations_with_replacement(s, size=2)),
            ("product", lambda s: pt.product(s, repeat=1)),
            ("sort", lambda s: pt.sort(s)),
        ],
    )
    def test_the_other_functions_that_split_no_longer_raise(self, name: str, call):
        assert call(self.SEQ)

    def test_a_partial_slice_of_the_same_peptide_still_refuses(self):
        # The boundary: splitting keeps the whole peptide, slicing [0:4] keeps the primary tag and
        # throws the reference away, so only the second is a loss and only the second raises.
        assert pt.split(self.SEQ)[3] == "T[Phospho#g1]"
        with pytest.raises(PeptacularError, match="#g1"):
            pt.parse(self.SEQ)[0:4]

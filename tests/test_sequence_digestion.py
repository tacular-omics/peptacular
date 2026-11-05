"""Coverage for the functional digestion API (peptacular.sequence.digestion).

Each function dispatches on scalar vs. batch input; these tests exercise both
paths and confirm they agree with the underlying ProFormaAnnotation methods.
"""

import re

import pytest

import peptacular as pt
from peptacular.sequence.digestion import (
    cleavage_sites,
    digest,
    left_semi_digest,
    nonspecific_digest,
    right_semi_digest,
    semi_digest,
    simple_cleavage_sites,
    simple_digest,
)

SEQ = "TIDERTIDEKTIDE"
PARENT = "PEPTIDEKTIDERTIDE"


class TestLeftSemiDigest:
    def test_scalar(self):
        result = left_semi_digest(SEQ, min_len=2, max_len=5)
        assert ("TIDER", pt.spans.Span(0, 5, 0)) in result
        assert all(2 <= len(s) <= 5 for s, _ in result)

    def test_batch_matches_scalar(self):
        scalar = left_semi_digest(SEQ, min_len=2, max_len=5)
        batch = left_semi_digest([SEQ, SEQ], min_len=2, max_len=5)
        assert batch == [scalar, scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert left_semi_digest(a, min_len=2, max_len=5) == left_semi_digest(SEQ, min_len=2, max_len=5)

    def test_default_bounds(self):
        assert len(left_semi_digest(SEQ)) > 0


class TestRightSemiDigest:
    def test_scalar(self):
        result = right_semi_digest(SEQ, min_len=2, max_len=5)
        assert all(2 <= len(s) <= 5 for s, _ in result)

    def test_batch_matches_scalar(self):
        scalar = right_semi_digest(SEQ, min_len=2, max_len=5)
        batch = right_semi_digest([SEQ], min_len=2, max_len=5)
        assert batch == [scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert right_semi_digest(a, min_len=2, max_len=5) == right_semi_digest(SEQ, min_len=2, max_len=5)


class TestSemiDigest:
    def test_scalar_is_union_of_left_and_right(self):
        left = {s for s, _ in left_semi_digest(SEQ, min_len=2, max_len=5)}
        right = {s for s, _ in right_semi_digest(SEQ, min_len=2, max_len=5)}
        both = {s for s, _ in semi_digest(SEQ, min_len=2, max_len=5)}
        assert left <= both and right <= both

    def test_batch(self):
        assert len(semi_digest([SEQ, SEQ], min_len=2, max_len=5)) == 2

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert semi_digest(a, min_len=2, max_len=5) == semi_digest(SEQ, min_len=2, max_len=5)


class TestNonspecificDigest:
    def test_scalar(self):
        result = nonspecific_digest(SEQ, min_len=2, max_len=3)
        assert all(2 <= len(s) <= 3 for s, _ in result)

    def test_batch_matches_scalar(self):
        scalar = nonspecific_digest(SEQ, min_len=2, max_len=3)
        batch = nonspecific_digest([SEQ, SEQ], min_len=2, max_len=3)
        assert batch == [scalar, scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert nonspecific_digest(a, min_len=2, max_len=3) == nonspecific_digest(SEQ, min_len=2, max_len=3)


class TestCleavageSites:
    def test_scalar_regex(self):
        assert cleavage_sites(SEQ, re.compile("(?<=[KR])")) == [5, 10]

    def test_batch(self):
        assert cleavage_sites([SEQ, SEQ], re.compile("(?<=[KR])")) == [[5, 10], [5, 10]]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert cleavage_sites(a, re.compile("(?<=[KR])")) == cleavage_sites(SEQ, re.compile("(?<=[KR])"))


class TestSimpleCleavageSites:
    def test_scalar(self):
        assert simple_cleavage_sites(SEQ, cleave_on="KR") == [5, 10]

    def test_restrict_and_cterminal_kwargs(self):
        sites = simple_cleavage_sites(SEQ, cleave_on="KR", restrict_before="", restrict_after="", cterminal=True)
        assert sites == [5, 10]

    def test_batch(self):
        assert simple_cleavage_sites([SEQ, SEQ], cleave_on="KR") == [[5, 10], [5, 10]]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert simple_cleavage_sites(a, cleave_on="KR") == simple_cleavage_sites(SEQ, cleave_on="KR")


class TestDigest:
    def test_scalar_regex(self):
        result = digest(SEQ, re.compile("(?<=[KR])"), missed_cleavages=1, min_len=1, max_len=100)
        peptides = {s for s, _ in result}
        assert "TIDER" in peptides
        assert "TIDERTIDEK" in peptides

    def test_batch_matches_scalar(self):
        scalar = digest(SEQ, re.compile("(?<=[KR])"), missed_cleavages=0, min_len=1, max_len=100)
        batch = digest([SEQ, SEQ], re.compile("(?<=[KR])"), missed_cleavages=0, min_len=1, max_len=100)
        assert batch == [scalar, scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert digest(a, re.compile("(?<=[KR])"), min_len=1, max_len=100) == digest(SEQ, re.compile("(?<=[KR])"), min_len=1, max_len=100)

    def test_batch_with_explicit_parallel_kwargs(self):
        result = digest([SEQ, SEQ], re.compile("(?<=[KR])"), missed_cleavages=0, min_len=1, max_len=100, n_workers=2, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_semi_flag(self):
        semi_result = digest(SEQ, re.compile("(?<=[KR])"), semi=True, min_len=1, max_len=100)
        non_semi_result = digest(SEQ, re.compile("(?<=[KR])"), semi=False, min_len=1, max_len=100)
        assert len(semi_result) >= len(non_semi_result)


class TestSimpleDigest:
    def test_scalar(self):
        result = simple_digest(SEQ, cleave_on="KR", missed_cleavages=1, min_len=1, max_len=100)
        peptides = {s for s, _ in result}
        assert "TIDER" in peptides
        assert "TIDERTIDEK" in peptides

    def test_batch_matches_scalar(self):
        scalar = simple_digest(SEQ, cleave_on="KR", missed_cleavages=0, min_len=1, max_len=100)
        batch = simple_digest([SEQ, SEQ], cleave_on="KR", missed_cleavages=0, min_len=1, max_len=100)
        assert batch == [scalar, scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert simple_digest(a, cleave_on="KR", min_len=1, max_len=100) == simple_digest(SEQ, cleave_on="KR", min_len=1, max_len=100)

    def test_restrict_and_cterminal_kwargs(self):
        result = simple_digest(SEQ, cleave_on="KR", restrict_before="", restrict_after="", cterminal=True, min_len=1, max_len=100)
        assert len(result) > 0


class TestDigestProductAnnotations:
    """A digestion product is a different molecule from the sequence it was cut out of.

    The parent's charge and its labile and unlocalised modifications do not belong to any
    one product, so copying them onto every product gives every product the wrong mass.
    These tests assert masses rather than serialized strings: a product can serialize
    plausibly and still weigh the wrong amount, and the reverse.
    """

    WATER = pt.chem_mass("H2O")

    def _product_masses(self, sequence: str) -> list[float]:
        return [pt.parse(peptide).neutral_mass() for peptide, _ in digest(sequence, "trypsin")]

    def test_charge_is_not_copied_onto_products(self):
        assert self._product_masses(f"{PARENT}/2") == pytest.approx(self._product_masses(PARENT))

    def test_adduct_charge_is_not_copied_onto_products(self):
        assert self._product_masses(f"{PARENT}/[Na:z+1]") == pytest.approx(self._product_masses(PARENT))

    def test_labile_mod_is_not_copied_onto_products(self):
        assert self._product_masses("{Glycan:Hex}" + PARENT) == pytest.approx(self._product_masses(PARENT))

    def test_unlocalised_mod_is_not_copied_onto_products(self):
        assert self._product_masses("[Phospho]?" + PARENT) == pytest.approx(self._product_masses(PARENT))

    @pytest.mark.parametrize(
        ("sequence", "residue_level"),
        [
            (PARENT, PARENT),
            (f"{PARENT}/3", PARENT),
            (f"{PARENT}/[Na:z+1]", PARENT),
            ("{Glycan:Hex}" + PARENT, PARENT),
            ("[Phospho]?" + PARENT, PARENT),
            ("<13C>" + PARENT, "<13C>" + PARENT),
            ("[Acetyl]-" + PARENT + "-[Amidated]", "[Acetyl]-" + PARENT + "-[Amidated]"),
            ("PEM[Oxidation]TIDEKTIDERTIDE", "PEM[Oxidation]TIDEKTIDERTIDE"),
        ],
    )
    def test_product_masses_balance_against_the_parent(self, sequence, residue_level):
        """Hydrolysing n-1 peptide bonds adds n-1 waters and nothing else.

        ``residue_level`` is the parent written without the annotations a product cannot
        inherit, so the balance also states which of them are dropped and which are kept.
        """
        masses = self._product_masses(sequence)
        expected = pt.parse(residue_level).neutral_mass() + (len(masses) - 1) * self.WATER
        assert sum(masses) == pytest.approx(expected)

    def test_isotope_label_is_kept(self):
        """A global label applies to every residue, so every product keeps it."""
        labelled = self._product_masses("<13C>" + PARENT)
        unlabelled = self._product_masses(PARENT)
        assert all(a > b for a, b in zip(labelled, unlabelled, strict=True))

    def test_static_rule_is_kept_and_applies_per_residue(self):
        """A static rule applies to whichever products carry the residue it targets."""
        products = dict(digest("<[Carbamidomethyl]@C>PEPCTIDEKTIDER", "trypsin"))
        masses = {pt.parse(peptide).stripped_sequence: pt.parse(peptide).neutral_mass() for peptide in products}
        assert masses["PEPCTIDEK"] == pytest.approx(pt.parse("PEPCTIDEK").neutral_mass() + pt.chem_mass("C2H3NO"))
        assert masses["TIDER"] == pytest.approx(pt.parse("TIDER").neutral_mass())

    def test_terminal_mods_stay_on_their_own_product(self):
        masses = self._product_masses("[Acetyl]-" + PARENT + "-[Amidated]")
        bare = self._product_masses(PARENT)
        assert masses[0] == pytest.approx(bare[0] + pt.chem_mass("C2H2O"))
        assert masses[-1] == pytest.approx(bare[-1] - pt.chem_mass("O") + pt.chem_mass("NH"))
        assert masses[1:-1] == pytest.approx(bare[1:-1])

    def test_digesting_does_not_mutate_the_caller_annotation(self):
        annotation = pt.parse("{Glycan:Hex}[Phospho]?" + PARENT + "/2")
        before = annotation.serialize()
        digest(annotation, "trypsin")
        assert annotation.serialize() == before

    @pytest.mark.parametrize(
        ("function", "kwargs"),
        [
            (digest, {"enzyme": "trypsin"}),
            (simple_digest, {"cleave_on": "KR"}),
            (semi_digest, {"min_len": 3}),
            (nonspecific_digest, {"min_len": 3}),
            (left_semi_digest, {"min_len": 3}),
            (right_semi_digest, {"min_len": 3}),
        ],
    )
    def test_every_digest_drops_the_parent_analyte_annotations(self, function, kwargs):
        charged = "{Glycan:Hex}[Phospho]?PEPTIDEKTIDER/2"
        for peptide, span in function(charged, **kwargs):
            product = pt.parse(peptide)
            assert product.charge is None
            assert not product.has_labile_mods
            assert not product.has_unknown_mods
            assert product.neutral_mass() == pytest.approx(pt.parse("PEPTIDEKTIDER")[span].neutral_mass())

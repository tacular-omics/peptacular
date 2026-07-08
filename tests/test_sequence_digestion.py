"""Coverage for the functional digestion API (peptacular.sequence.digestion).

Each function dispatches on scalar vs. batch input; these tests exercise both
paths and confirm they agree with the underlying ProFormaAnnotation methods.
"""

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
        assert cleavage_sites(SEQ, "(?<=[KR])") == [5, 10]

    def test_batch(self):
        assert cleavage_sites([SEQ, SEQ], "(?<=[KR])") == [[5, 10], [5, 10]]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert cleavage_sites(a, "(?<=[KR])") == cleavage_sites(SEQ, "(?<=[KR])")


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
        result = digest(SEQ, "(?<=[KR])", missed_cleavages=1, min_len=1, max_len=100)
        peptides = {s for s, _ in result}
        assert "TIDER" in peptides
        assert "TIDERTIDEK" in peptides

    def test_batch_matches_scalar(self):
        scalar = digest(SEQ, "(?<=[KR])", missed_cleavages=0, min_len=1, max_len=100)
        batch = digest([SEQ, SEQ], "(?<=[KR])", missed_cleavages=0, min_len=1, max_len=100)
        assert batch == [scalar, scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert digest(a, "(?<=[KR])", min_len=1, max_len=100) == digest(SEQ, "(?<=[KR])", min_len=1, max_len=100)

    def test_batch_with_explicit_parallel_kwargs(self):
        result = digest(
            [SEQ, SEQ], "(?<=[KR])", missed_cleavages=0, min_len=1, max_len=100, n_workers=2, chunksize=1, method="sequential"
        )
        assert len(result) == 2

    def test_semi_flag(self):
        semi_result = digest(SEQ, "(?<=[KR])", semi=True, min_len=1, max_len=100)
        non_semi_result = digest(SEQ, "(?<=[KR])", semi=False, min_len=1, max_len=100)
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

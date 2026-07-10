"""Coverage for the functional combinatoric API (peptacular.sequence.combinatoric):
permutations, product, combinations, combinations_with_replacement.
"""

import peptacular as pt
from peptacular.sequence.combinatoric import combinations, combinations_with_replacement, permutations, product

SEQ = "PET"


class TestPermutations:
    def test_scalar(self):
        assert permutations(SEQ) == ["PET", "PTE", "EPT", "ETP", "TPE", "TEP"]

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert permutations(a) == permutations(SEQ)

    def test_batch(self):
        result = permutations([SEQ, SEQ])
        assert result == [permutations(SEQ), permutations(SEQ)]

    def test_batch_with_parallel_kwargs(self):
        result = permutations([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_size_kwarg(self):
        result = permutations(SEQ, size=2)
        assert all(len(pt.strip_mods(r)) == 2 for r in result)


class TestProduct:
    def test_scalar(self):
        result = product(SEQ, 2)
        assert result[0] == "PP"
        assert len(result) == 9

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert product(a, 2) == product(SEQ, 2)

    def test_batch(self):
        result = product([SEQ, SEQ], 2)
        assert result == [product(SEQ, 2), product(SEQ, 2)]

    def test_batch_with_parallel_kwargs(self):
        result = product([SEQ, SEQ], 2, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_repeat_none_defaults_to_length(self):
        result = product(SEQ, None)
        assert len(result) == 3**3


class TestCombinations:
    def test_scalar(self):
        assert combinations(SEQ, 2) == ["PE", "PT", "ET"]

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert combinations(a, 2) == combinations(SEQ, 2)

    def test_batch(self):
        result = combinations([SEQ, SEQ], 2)
        assert result == [combinations(SEQ, 2), combinations(SEQ, 2)]

    def test_batch_with_parallel_kwargs(self):
        result = combinations([SEQ, SEQ], 2, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_size_none_defaults_to_length(self):
        assert combinations(SEQ, None) == ["PET"]


class TestCombinationsWithReplacement:
    def test_scalar(self):
        assert combinations_with_replacement(SEQ, 2) == ["PP", "PE", "PT", "EE", "ET", "TT"]

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert combinations_with_replacement(a, 2) == combinations_with_replacement(SEQ, 2)

    def test_batch(self):
        result = combinations_with_replacement([SEQ, SEQ], 2)
        assert result == [combinations_with_replacement(SEQ, 2), combinations_with_replacement(SEQ, 2)]

    def test_batch_with_parallel_kwargs(self):
        result = combinations_with_replacement([SEQ, SEQ], 2, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_size_none_defaults_to_length(self):
        result = combinations_with_replacement(SEQ, None)
        assert all(len(r) == len(SEQ) for r in result)
        assert result[0] == "PPP"

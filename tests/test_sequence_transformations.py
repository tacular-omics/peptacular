"""Coverage for the functional sequence-transformation API (peptacular.sequence.transformations):
reverse, shuffle, shift, span_to_sequence, split, sort, join.
"""

import peptacular as pt
from peptacular.sequence.transformations import join, reverse, shift, shuffle, sort, span_to_sequence, split
from peptacular.spans import Span

SEQ = "PEPTIDE"


class TestReverse:
    def test_scalar(self):
        assert reverse(SEQ) == "EDITPEP"

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert reverse(a) == reverse(SEQ)

    def test_keep_nterm(self):
        assert reverse(SEQ, keep_nterm=2) == "PEEDITP"

    def test_batch(self):
        assert reverse([SEQ, SEQ]) == ["EDITPEP", "EDITPEP"]

    def test_batch_with_parallel_kwargs(self):
        result = reverse([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert result == ["EDITPEP", "EDITPEP"]


class TestShuffle:
    def test_scalar(self):
        result = shuffle(SEQ, seed=0)
        assert sorted(result) == sorted(SEQ)

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert shuffle(a, seed=0) == shuffle(SEQ, seed=0)

    def test_keep_nterm(self):
        result = shuffle(SEQ, seed=0, keep_nterm=2)
        assert result.startswith("PE")

    def test_batch(self):
        result = shuffle([SEQ, SEQ], seed=0)
        assert result[0] == result[1]

    def test_batch_with_parallel_kwargs(self):
        result = shuffle([SEQ, SEQ], seed=0, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestShift:
    def test_scalar(self):
        assert shift(SEQ, 2) == "PTIDEPE"

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert shift(a, 2) == shift(SEQ, 2)

    def test_keep_nterm(self):
        assert shift(SEQ, 2, keep_nterm=2) == "PEIDEPT"

    def test_batch(self):
        assert shift([SEQ, SEQ], 2) == ["PTIDEPE", "PTIDEPE"]

    def test_batch_with_parallel_kwargs(self):
        result = shift([SEQ, SEQ], 2, n_workers=1, chunksize=1, method="sequential")
        assert result == ["PTIDEPE", "PTIDEPE"]


class TestSpanToSequence:
    def test_scalar_tuple(self):
        assert span_to_sequence(SEQ, (0, 4, 0)) == "PEPT"

    def test_scalar_span_object(self):
        assert span_to_sequence(SEQ, Span(0, 4, 0)) == "PEPT"

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert span_to_sequence(a, (0, 4, 0)) == span_to_sequence(SEQ, (0, 4, 0))

    def test_batch(self):
        result = span_to_sequence([SEQ, SEQ], (0, 4, 0))
        assert result == ["PEPT", "PEPT"]

    def test_batch_with_parallel_kwargs(self):
        result = span_to_sequence([SEQ, SEQ], (0, 4, 0), n_workers=1, chunksize=1, method="sequential")
        assert result == ["PEPT", "PEPT"]


class TestSplit:
    def test_scalar(self):
        assert split(SEQ) == ["P", "E", "P", "T", "I", "D", "E"]

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert split(a) == split(SEQ)

    def test_batch(self):
        result = split([SEQ, SEQ])
        assert result == [split(SEQ), split(SEQ)]

    def test_batch_with_parallel_kwargs(self):
        result = split([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestSort:
    def test_scalar(self):
        assert sort(SEQ) == "DEEIPPT"

    def test_scalar_annotation(self):
        a = pt.parse(SEQ)
        assert sort(a) == sort(SEQ)

    def test_reverse_kwarg(self):
        assert sort(SEQ, reverse=True) == "".join(sorted(SEQ, reverse=True))

    def test_key_kwarg(self):
        result = sort(SEQ, key=lambda c: -ord(c[0]))
        assert result == "".join(sorted(SEQ, key=lambda c: -ord(c)))

    def test_batch(self):
        assert sort([SEQ, SEQ]) == ["DEEIPPT", "DEEIPPT"]

    def test_batch_with_parallel_kwargs(self):
        result = sort([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert result == ["DEEIPPT", "DEEIPPT"]


class TestJoin:
    def test_scalar_str_list(self):
        assert join(["PEPTIDE", "MODIFIED"]) == "PEPTIDEMODIFIED"

    def test_scalar_annotation_list(self):
        annots = [pt.parse("PEPTIDE"), pt.parse("MODIFIED")]
        assert join(annots) == "PEPTIDEMODIFIED"

    def test_batch(self):
        result = join([["PEPTIDE", "MODIFIED"], ["PEPTIDE", "MODIFIED"]])
        assert result == ["PEPTIDEMODIFIED", "PEPTIDEMODIFIED"]

    def test_batch_with_parallel_kwargs(self):
        result = join(
            [["PEPTIDE", "MODIFIED"], ["PEPTIDE", "MODIFIED"]], n_workers=1, chunksize=1, method="sequential"
        )
        assert result == ["PEPTIDEMODIFIED", "PEPTIDEMODIFIED"]

"""Coverage for the functional fragmentation API (peptacular.sequence.fragmentation):
fragment(), frag(), and fast_fragment() scalar-vs-batch dispatch.
"""

import peptacular as pt
from peptacular.sequence.fragmentation import fast_fragment, frag, fragment

SEQ = "PEPTIDE"


class TestFragment:
    def test_scalar(self):
        result = fragment(SEQ, ion_types=["b", "y"], charges=[1])
        assert len(result) == 14  # full 1..n ladder: 7 b-ions + 7 y-ions

    def test_batch_matches_scalar(self):
        scalar = fragment(SEQ, ion_types=["b", "y"], charges=[1])
        batch = fragment([SEQ, SEQ], ion_types=["b", "y"], charges=[1])
        assert len(batch) == 2
        assert [f.mass for f in batch[0]] == [f.mass for f in scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        scalar = fragment(SEQ, ion_types=["b"], charges=[1])
        from_annot = fragment(a, ion_types=["b"], charges=[1])
        assert [f.mass for f in from_annot] == [f.mass for f in scalar]

    def test_charge_default_matches_precursor(self):
        # Regression: fragment() must use the OOP smart charge default, not charges=(1,).
        oop = pt.parse("PEPTIDE/3").fragment(ion_types=["b"])
        func = fragment("PEPTIDE/3", ion_types=["b"])
        assert sorted(f.charge_state for f in oop) == sorted(f.charge_state for f in func)

    def test_batch_with_parallel_kwargs(self):
        result = fragment([SEQ, SEQ], ion_types=["b"], charges=[1], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestFrag:
    def test_scalar(self):
        f = frag(SEQ, ion_type="y", charge=1, position=3)
        assert f.mass > 0

    def test_batch_matches_scalar(self):
        scalar = frag(SEQ, ion_type="y", charge=1, position=3)
        batch = frag([SEQ, SEQ], ion_type="y", charge=1, position=3)
        assert len(batch) == 2
        assert batch[0].mass == batch[1].mass == scalar.mass

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        scalar = frag(SEQ, ion_type="y", charge=1, position=3)
        from_annot = frag(a, ion_type="y", charge=1, position=3)
        assert from_annot.mass == scalar.mass

    def test_batch_with_parallel_kwargs(self):
        result = frag([SEQ, SEQ], ion_type="y", charge=1, position=3, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestFastFragment:
    def test_scalar(self):
        result = fast_fragment(SEQ, ion_types=["b", "y"], charges=[1])
        assert isinstance(result, dict)
        assert len(result) == 2  # (B,1) and (Y,1) keys

    def test_batch_matches_scalar(self):
        scalar = fast_fragment(SEQ, ion_types=["b", "y"], charges=[1])
        batch = fast_fragment([SEQ, SEQ], ion_types=["b", "y"], charges=[1])
        assert len(batch) == 2
        assert batch[0] == batch[1] == scalar

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        scalar = fast_fragment(SEQ, ion_types=["b"], charges=[1])
        from_annot = fast_fragment(a, ion_types=["b"], charges=[1])
        assert from_annot == scalar

    def test_batch_with_parallel_kwargs(self):
        result = fast_fragment([SEQ, SEQ], ion_types=["b"], charges=[1], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

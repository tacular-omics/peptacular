"""Coverage for the functional mass/mz/composition API (peptacular.sequence.mass_funcs)."""

from tacular import IonType

import peptacular as pt
from peptacular.sequence.mass_funcs import comp, mass, mz

SEQ = "PEPTIDE"


class TestMass:
    def test_scalar(self):
        assert mass(SEQ) > 0

    def test_scalar_annotation(self):
        assert mass(pt.parse(SEQ)) == mass(SEQ)

    def test_batch_matches_scalar(self):
        scalar = mass(SEQ)
        batch = mass([SEQ, SEQ])
        assert batch == [scalar, scalar]

    def test_batch_with_parallel_kwargs(self):
        result = mass([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_ion_type_and_charge(self):
        assert mass(SEQ, ion_type=IonType.Y, charge=2) > 0

    def test_monoisotopic_false(self):
        assert mass(SEQ, monoisotopic=False) != mass(SEQ, monoisotopic=True)

    def test_calculate_with_composition(self):
        assert abs(mass(SEQ, calculate_with_composition=True) - mass(SEQ)) < 1e-6


class TestMz:
    def test_scalar(self):
        assert mz(SEQ, charge=2) > 0

    def test_scalar_annotation(self):
        assert mz(pt.parse(SEQ), charge=2) == mz(SEQ, charge=2)

    def test_batch_matches_scalar(self):
        scalar = mz(SEQ, charge=2)
        batch = mz([SEQ, SEQ], charge=2)
        assert batch == [scalar, scalar]

    def test_batch_with_parallel_kwargs(self):
        result = mz([SEQ, SEQ], charge=2, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_ion_type(self):
        assert mz(SEQ, ion_type=IonType.B, charge=1) > 0


class TestComp:
    def test_scalar(self):
        result = comp(SEQ)
        assert sum(result.values()) > 0

    def test_scalar_annotation(self):
        assert comp(pt.parse(SEQ)) == comp(SEQ)

    def test_batch_matches_scalar(self):
        scalar = comp(SEQ)
        batch = comp([SEQ, SEQ])
        assert batch == [scalar, scalar]

    def test_batch_with_parallel_kwargs(self):
        result = comp([SEQ, SEQ], n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2

    def test_ion_type_and_charge(self):
        result = comp(SEQ, ion_type=IonType.Y, charge=1)
        assert sum(result.values()) > 0

    def test_does_not_mutate_input_annotation(self):
        a = pt.parse(SEQ)
        before = a.serialize()
        comp(a)
        assert a.serialize() == before

"""Coverage for the functional isotopic-distribution API (peptacular.sequence.isotope)."""

import pytest

import peptacular as pt
import peptacular.sequence
from peptacular.sequence.isotope import estimate_isotopic_distribution, isotopic_distribution

SEQ = "PEPTIDE"


class TestIsotopicDistribution:
    def test_scalar(self):
        result = isotopic_distribution(SEQ, max_isotopes=3)
        assert len(result) == 3
        assert result[0].mass > 0

    def test_batch_matches_scalar(self):
        scalar = isotopic_distribution(SEQ, max_isotopes=3)
        batch = isotopic_distribution([SEQ, SEQ], max_isotopes=3)
        assert batch == [scalar, scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert isotopic_distribution(a, max_isotopes=3) == isotopic_distribution(SEQ, max_isotopes=3)

    def test_charge_and_ion_type_kwargs(self):
        result = isotopic_distribution(SEQ, ion_type="y", charge=2, max_isotopes=2)
        assert len(result) == 2

    def test_batch_with_parallel_kwargs(self):
        result = isotopic_distribution([SEQ, SEQ], max_isotopes=2, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestEstimateIsotopicDistribution:
    def test_scalar(self):
        result = estimate_isotopic_distribution(SEQ, max_isotopes=3)
        assert len(result) == 3

    def test_batch_matches_scalar(self):
        scalar = estimate_isotopic_distribution(SEQ, max_isotopes=3)
        batch = estimate_isotopic_distribution([SEQ, SEQ], max_isotopes=3)
        assert batch == [scalar, scalar]

    def test_annotation_input(self):
        a = pt.parse(SEQ)
        assert estimate_isotopic_distribution(a, max_isotopes=3) == estimate_isotopic_distribution(SEQ, max_isotopes=3)

    def test_mono_matches_exact_distribution(self):
        # The estimate's monoisotopic peak should land on the true neutral mass
        # (regression coverage for the averagine-anchoring fix).
        exact = isotopic_distribution(SEQ, max_isotopes=1, min_abundance_threshold=0.0)
        est = estimate_isotopic_distribution(SEQ, max_isotopes=1, min_abundance_threshold=0.0)
        assert abs(exact[0].mass - est[0].mass) < 1e-4

    def test_batch_with_parallel_kwargs(self):
        result = estimate_isotopic_distribution([SEQ, SEQ], max_isotopes=2, n_workers=1, chunksize=1, method="sequential")
        assert len(result) == 2


class TestEstimateIsotopicDistributionIsReachable:
    """The two functions named ``estimate_isotopic_distribution`` are reachable and distinct.

    ``peptacular.sequence`` holds the sequence-based estimate, documented in
    ``docs/migration.rst`` and mirrored by ``ProFormaAnnotation.estimate_isotopic_distribution``.
    ``pt.estimate_isotopic_distribution`` is the averagine estimate from a neutral mass,
    which is what it has been since 4.x.
    """

    def test_imports_from_the_public_sequence_namespace(self):
        from peptacular.sequence import estimate_isotopic_distribution as public

        assert "estimate_isotopic_distribution" in peptacular.sequence.__all__
        assert public is estimate_isotopic_distribution
        assert public(SEQ, max_isotopes=2) == estimate_isotopic_distribution(SEQ, max_isotopes=2)

    def test_pt_name_takes_a_neutral_mass(self):
        mass = pt.parse(SEQ).neutral_mass()
        from_mass = pt.estimate_isotopic_distribution(mass, max_isotopes=2)
        from_sequence = estimate_isotopic_distribution(SEQ, max_isotopes=2)
        assert [peak.mass for peak in from_mass] == pytest.approx([peak.mass for peak in from_sequence])

"""
Tests for isotopic distribution charge handling.
"""

import peptacular as pt
from peptacular.constants import ELECTRON_MASS
from peptacular.isotope import _convolve_distributions, estimate_isotopic_distribution, isotopic_distribution


class TestAveragineAnchoring:
    """The averagine estimate supplies only the shape; its monoisotopic peak must
    be anchored to the requested mass (regression: it used to drift by 10-34 Da)."""

    def test_monoisotopic_peak_lands_on_requested_mass(self):
        # distribution_resolution=None avoids peak-mass rounding, so anchoring is exact.
        for target in (500.0, 800.0, 1500.0, 3000.0):
            dist = estimate_isotopic_distribution(target, max_isotopes=5, min_abundance_threshold=0.0, distribution_resolution=None)
            assert abs(dist[0].mass - target) < 1e-6, f"{target}: got {dist[0].mass}"

    def test_default_resolution_lands_within_rounding(self):
        # With distribution_resolution=5 the anchor is exact up to the peak rounding.
        dist = estimate_isotopic_distribution(800.0, max_isotopes=5, min_abundance_threshold=0.0)
        assert abs(dist[0].mass - 800.0) < 1e-4

    def test_neutron_count_mode_unaffected(self):
        # In neutron-count mode masses are offsets (0, 1, 2, ...) and must not be shifted.
        dist = estimate_isotopic_distribution(800.0, max_isotopes=3, min_abundance_threshold=0.0, use_neutron_count=True)
        assert dist[0].mass == 0.0

    def test_estimate_mono_matches_exact_mono(self):
        for seq in ("PEPTIDE", "PEPTIDEKR"):
            a = pt.parse(seq)
            exact = a.isotopic_distribution(max_isotopes=1, min_abundance_threshold=0.0, distribution_resolution=None)[0].mass
            est = a.estimate_isotopic_distribution(max_isotopes=1, min_abundance_threshold=0.0, distribution_resolution=None)[0].mass
            assert abs(exact - est) < 1e-4, f"{seq}: exact={exact} est={est}"


class TestConvolveThreshold:
    """Regression test for the abundance-threshold pruning in the convolution."""

    def test_low_abundance_pairing_does_not_drop_later_peaks(self):
        """A below-threshold pairing must not abort the whole inner loop.

        ``dist2`` is not sorted by descending abundance, so the old ``break`` would
        exit on the first (low-abundance) entry and discard the following
        high-abundance entry. Only the below-threshold pairing should be skipped.
        """
        dist1 = {0.0: (1.0, 0)}
        # First entry is below-threshold when multiplied; second is well above.
        dist2 = {1.0: (0.1, 0), 2.0: (0.9, 1)}
        result = _convolve_distributions(dist1, dist2, None, 0.5, 5)
        # The high-abundance pairing (2.0) must survive; the low one (1.0) is pruned.
        assert 2.0 in result
        assert 1.0 not in result


class TestIsotopeChargeState:
    """Regression tests for the ``charge_state`` electron-mass correction.

    Previously the per-charge electron-mass correction was only applied when the
    formula contained non-integer counts (``delta_mass != 0``), so for ordinary
    integer formulas ``charge_state`` was silently ignored and every charge state
    returned identical masses.
    """

    def test_charge_state_shifts_mass_by_electron(self):
        """A positive charge removes one electron mass per charge from every peak."""
        formula = {"C": 12, "H": 6, "N": 3}
        neutral = isotopic_distribution(formula)
        charged1 = isotopic_distribution(formula, charge_state=1)
        charged2 = isotopic_distribution(formula, charge_state=2)

        assert abs((neutral[0].mass - charged1[0].mass) - ELECTRON_MASS) < 1e-9
        assert abs((neutral[0].mass - charged2[0].mass) - 2 * ELECTRON_MASS) < 1e-9

    def test_charge_state_none_and_zero_are_uncorrected(self):
        """``charge_state`` of ``None`` or ``0`` leaves masses unchanged."""
        formula = {"C": 12, "H": 6, "N": 3}
        base = isotopic_distribution(formula, charge_state=None)
        zero = isotopic_distribution(formula, charge_state=0)
        assert [round(x.mass, 8) for x in base] == [round(x.mass, 8) for x in zero]

    def test_charge_states_differ(self):
        """Different charge states must not produce identical masses."""
        formula = {"C": 12, "H": 6, "N": 3}
        m1 = [x.mass for x in isotopic_distribution(formula, charge_state=1)]
        m2 = [x.mass for x in isotopic_distribution(formula, charge_state=2)]
        assert m1 != m2

    def test_annotation_isotope_mono_matches_mass(self):
        """The mono-isotopic peak must match the electron-corrected charged mass."""
        annot = pt.parse("PEPTIDE")
        for charge in (1, 2):
            expected = annot.mass(charge=charge)
            mono = annot.isotopic_distribution(charge=charge, distribution_resolution=None)[0].mass
            assert abs(expected - mono) < 1e-5


def test_annotation_isotopes_use_total_intrinsic_and_external_charge():
    from peptacular import parse
    from peptacular.isotope import isotopic_distribution

    annotation = parse("PEP[Formula:CH2:z+1]TIDE/2")
    fragment = annotation.frag(calculate_composition=True)
    assert fragment.charge_state == 3
    expected = isotopic_distribution(fragment.composition, charge_state=3)
    assert annotation.isotopic_distribution() == expected

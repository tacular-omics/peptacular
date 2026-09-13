"""Tests for aggregated BRAIN isotope distributions."""

from collections.abc import Mapping

import pytest
from tacular import ELEMENT_LOOKUP

import peptacular as pt
from peptacular.constants import ELECTRON_MASS
from peptacular.isotope import (
    _brain_coefficients,
    brain_isotopic_distribution,
    estimate_isotopic_distribution,
    isotopic_distribution,
)


def _reference_distribution(composition: Mapping[str, int], length: int) -> tuple[list[float], list[float]]:
    """Directly convolve probabilities and exact mass moments."""

    distribution: dict[int, tuple[float, float]] = {0: (1.0, 0.0)}
    for element, count in composition.items():
        info = ELEMENT_LOOKUP[element]
        offsets = ELEMENT_LOOKUP.get_neutron_offsets_and_abundances(info)
        masses = ELEMENT_LOOKUP.get_masses_and_abundances(info)
        entries = [(int(offset), float(abundance), float(mass)) for (offset, abundance), (mass, _) in zip(offsets, masses, strict=True) if abundance > 0.0]
        minimum = min(offset for offset, _, _ in entries)
        atom = [(offset - minimum, abundance, mass) for offset, abundance, mass in entries]
        for _ in range(count):
            updated: dict[int, tuple[float, float]] = {}
            for prior_offset, (prior_probability, prior_moment) in distribution.items():
                for atom_offset, atom_probability, atom_mass in atom:
                    offset = prior_offset + atom_offset
                    if offset >= length:
                        continue
                    probability = prior_probability * atom_probability
                    moment = prior_moment * atom_probability + probability * atom_mass
                    old_probability, old_moment = updated.get(offset, (0.0, 0.0))
                    updated[offset] = old_probability + probability, old_moment + moment
            distribution = updated

    maximum = max(probability for probability, _ in distribution.values())
    probabilities = [distribution.get(index, (0.0, 0.0))[0] / maximum for index in range(length)]
    centers = [moment / probability if probability else 0.0 for probability, moment in (distribution.get(index, (0.0, 0.0)) for index in range(length))]
    return probabilities, centers


@pytest.mark.parametrize(
    "composition",
    [
        {"C": 6, "H": 12, "O": 6},
        {"C": 12, "H": 22, "N": 2, "O": 11},
        {"C": 8, "H": 18, "N": 1, "O": 6, "P": 1, "S": 1},
    ],
)
def test_brain_matches_direct_convolution(composition):
    observed = brain_isotopic_distribution(composition, max_isotopes=16, min_abundance_threshold=0.0)
    expected_probability, expected_center = _reference_distribution(composition, 16)
    observed_by_offset = {peak.neutron_count: peak for peak in observed}
    for offset, probability in enumerate(expected_probability):
        if probability == 0.0:
            assert offset not in observed_by_offset
            continue
        peak = observed_by_offset[offset]
        assert peak.abundance == pytest.approx(probability, rel=1e-11, abs=1e-14)
        assert peak.mass == pytest.approx(expected_center[offset], rel=1e-12, abs=1e-8)


def test_aggregates_fine_structure_into_nominal_peaks():
    distribution = isotopic_distribution({"C": 12, "H": 6, "N": 3}, max_isotopes=3, min_abundance_threshold=0.0)
    assert [peak.neutron_count for peak in distribution] == [0, 1, 2]
    assert len(distribution) == 3


def test_sulfur_36_uses_four_neutron_offset():
    distribution = isotopic_distribution({"S": 1}, max_isotopes=5, min_abundance_threshold=0.0)
    assert [peak.neutron_count for peak in distribution] == [0, 1, 2, 4]
    assert distribution[-1].mass == pytest.approx(ELEMENT_LOOKUP["36S"].mass)


def test_elements_with_isotope_gaps_remain_sparse():
    distribution = isotopic_distribution({"Cl": 2}, max_isotopes=5, min_abundance_threshold=0.0)
    assert [peak.neutron_count for peak in distribution] == [0, 2, 4]


def test_fixed_isotope_labels_do_not_convolve():
    distribution = isotopic_distribution({"13C": 2}, max_isotopes=10, min_abundance_threshold=0.0)
    assert distribution == [pt.IsotopicData(ELEMENT_LOOKUP["13C"].mass * 2, 0, 1.0)]


def test_adaptive_high_mass_envelope_extends_past_old_limit():
    distribution = estimate_isotopic_distribution(50_000.0, min_abundance_threshold=0.01)
    assert len(distribution) > 32
    assert max(distribution, key=lambda peak: peak.abundance).neutron_count > 25
    assert distribution[-1].abundance >= 0.01


def test_adaptive_envelope_keeps_weak_leading_peaks():
    distribution = estimate_isotopic_distribution(50_000.0, min_abundance_threshold=0.01)
    assert distribution[0].neutron_count == 0
    assert distribution[0].abundance < 0.01


def test_formula_cache_is_bounded_and_public_results_are_distinct():
    first = isotopic_distribution({"C": 100}, max_isotopes=10)
    second = isotopic_distribution({"C": 100}, max_isotopes=10)
    assert first == second
    assert first is not second
    assert _brain_coefficients.cache_info().maxsize == 4096


class TestAveragineAnchoring:
    def test_monoisotopic_peak_lands_on_requested_mass(self):
        for target in (500.0, 800.0, 1500.0, 3000.0):
            distribution = estimate_isotopic_distribution(target, max_isotopes=5, min_abundance_threshold=0.0)
            assert distribution[0].mass == pytest.approx(target, abs=1e-8)

    def test_estimate_mono_matches_exact_mono(self):
        for sequence in ("PEPTIDE", "PEPTIDEKR"):
            annotation = pt.parse(sequence)
            exact = annotation.isotopic_distribution(max_isotopes=1, min_abundance_threshold=0.0)[0].mass
            estimated = annotation.estimate_isotopic_distribution(max_isotopes=1, min_abundance_threshold=0.0)[0].mass
            assert estimated == pytest.approx(exact, abs=1e-4)


class TestIsotopeChargeState:
    def test_charge_state_shifts_mass_by_electron(self):
        formula = {"C": 12, "H": 6, "N": 3}
        neutral = isotopic_distribution(formula)
        charged1 = isotopic_distribution(formula, charge_state=1)
        charged2 = isotopic_distribution(formula, charge_state=2)
        assert neutral[0].mass - charged1[0].mass == pytest.approx(ELECTRON_MASS)
        assert neutral[0].mass - charged2[0].mass == pytest.approx(2 * ELECTRON_MASS)

    def test_annotation_isotope_mono_matches_mass(self):
        annotation = pt.parse("PEPTIDE")
        for charge in (1, 2):
            expected = annotation.mass(charge=charge)
            observed = annotation.isotopic_distribution(charge=charge, max_isotopes=1)[0].mass
            assert observed == pytest.approx(expected, abs=1e-5)


def test_annotation_isotopes_use_total_intrinsic_and_external_charge():
    annotation = pt.parse("PEP[Formula:CH2:z+1]TIDE/2")
    fragment = annotation.frag(calculate_composition=True)
    assert fragment.charge_state == 3
    expected = isotopic_distribution(fragment.composition, charge_state=3)
    assert annotation.isotopic_distribution() == expected


@pytest.mark.parametrize("formula", [{"C": -1}, {"C": float("nan")}, {"C": True}])
def test_invalid_compositions_are_rejected(formula):
    with pytest.raises(ValueError):
        isotopic_distribution(formula)

"""Tests for aggregated BRAIN isotope distributions."""

from collections.abc import Mapping

import pytest
from tacular import ELEMENT_LOOKUP, ElementInfo

import peptacular as pt
from peptacular.constants import ELECTRON_MASS
from peptacular.isotope import (
    _brain_coefficients,
    brain_isotopic_distribution,
    estimate_isotopic_distribution,
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
    distribution = brain_isotopic_distribution({"C": 12, "H": 6, "N": 3}, max_isotopes=3, min_abundance_threshold=0.0)
    assert [peak.neutron_count for peak in distribution] == [0, 1, 2]
    assert len(distribution) == 3


def test_sulfur_36_uses_four_neutron_offset():
    distribution = brain_isotopic_distribution({"S": 1}, max_isotopes=5, min_abundance_threshold=0.0)
    assert [peak.neutron_count for peak in distribution] == [0, 1, 2, 4]
    assert distribution[-1].mass == pytest.approx(ELEMENT_LOOKUP["36S"].mass)


def test_elements_with_isotope_gaps_remain_sparse():
    distribution = brain_isotopic_distribution({"Cl": 2}, max_isotopes=5, min_abundance_threshold=0.0)
    assert [peak.neutron_count for peak in distribution] == [0, 2, 4]


def test_fixed_isotope_labels_do_not_convolve():
    distribution = brain_isotopic_distribution({"13C": 2}, max_isotopes=10, min_abundance_threshold=0.0)
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
    first = brain_isotopic_distribution({"C": 100}, max_isotopes=10)
    second = brain_isotopic_distribution({"C": 100}, max_isotopes=10)
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
        neutral = brain_isotopic_distribution(formula)
        charged1 = brain_isotopic_distribution(formula, charge=1)
        charged2 = brain_isotopic_distribution(formula, charge=2)
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
    fragment = annotation.frag(calculate_with_composition=True)
    assert fragment.charge_state == 3
    # two external protons each add the H binding term (PROTON_MASS, not H - e); the CH2 carrier does not
    binding = 2 * pt.HYDROGEN_BINDING_MASS
    expected = [
        pt.IsotopicData(mass=p.mass + binding, neutron_count=p.neutron_count, abundance=p.abundance)
        for p in brain_isotopic_distribution(fragment.composition, charge=3)
    ]
    assert annotation.isotopic_distribution() == expected
    assert annotation.isotopic_distribution()[0].mass == pytest.approx(fragment.mass, abs=1e-11)


@pytest.mark.parametrize("formula", [{"C": -1}, {"C": float("nan")}, {"C": True}])
def test_invalid_compositions_are_rejected(formula):
    with pytest.raises(ValueError):
        brain_isotopic_distribution(formula)


@pytest.mark.parametrize("mass", [-1, float("inf"), float("nan")])
@pytest.mark.parametrize("estimate", [pt.estimate_averagine_comp, estimate_isotopic_distribution])
def test_averagine_rejects_invalid_target_mass(estimate, mass):
    with pytest.raises(ValueError, match="finite and non-negative"):
        estimate(mass)


@pytest.mark.parametrize("limit", [0, -1, True, 2.5])
def test_isotope_window_requires_a_positive_integer(limit):
    with pytest.raises(ValueError, match="max_isotopes"):
        brain_isotopic_distribution({"C": 6}, max_isotopes=limit)


@pytest.mark.parametrize("threshold", [-0.1, 1.1, True, float("nan"), float("inf")])
def test_relative_abundance_threshold_is_validated(threshold):
    with pytest.raises(ValueError, match="min_abundance_threshold"):
        brain_isotopic_distribution({"C": 6}, min_abundance_threshold=threshold)


@pytest.mark.parametrize("charge", [True, 1.5, "2"])
def test_isotope_charge_requires_an_integer(charge):
    with pytest.raises(ValueError, match="charge must"):
        brain_isotopic_distribution({"C": 6}, charge=charge)


def test_zero_threshold_returns_complete_small_envelope_with_isotope_gaps():
    complete = brain_isotopic_distribution({"Cl": 2}, min_abundance_threshold=0)
    bounded = brain_isotopic_distribution({"Cl": 2}, max_isotopes=5, min_abundance_threshold=0)
    assert complete == bounded
    assert [peak.neutron_count for peak in complete] == [0, 2, 4]


def test_zero_threshold_large_envelope_requires_an_explicit_bound():
    with pytest.raises(ValueError, match="max_isotopes is required"):
        brain_isotopic_distribution({"C": 5000}, min_abundance_threshold=0)
    bounded = brain_isotopic_distribution({"C": 5000}, max_isotopes=8, min_abundance_threshold=0)
    assert len(bounded) == 8
    assert bounded[0].mass == pytest.approx(60_000)


@pytest.mark.parametrize("precision", [None, 2])
def test_merging_isotope_envelopes_preserves_total_abundance(precision):
    first = [pt.IsotopicData(101.004, 1, 0.25), pt.IsotopicData(100.001, 0, 1.0)]
    second = [pt.IsotopicData(100.001, 0, 0.5), pt.IsotopicData(101.003, 1, 0.125)]
    merged = pt.merge_isotopic_distributions(first, second, merge_precision=precision)
    assert sum(peak.abundance for peak in merged) == pytest.approx(1.875)
    assert [peak.mass for peak in merged] == sorted(peak.mass for peak in merged)
    assert merged[0].abundance == pytest.approx(1.5)
    if precision is None:
        assert len(merged) == 3
    else:
        assert merged == [pt.IsotopicData(100.0, 0, 1.5), pt.IsotopicData(101.0, 1, 0.375)]
    assert first[0] == pt.IsotopicData(101.004, 1, 0.25)


def test_merging_empty_envelopes_returns_no_peaks():
    assert pt.merge_isotopic_distributions() == []
    assert pt.merge_isotopic_distributions([], []) == []


def test_formula_isotopic_distribution_alias_removed():
    # 5.0: peptacular.isotope.isotopic_distribution duplicated brain_isotopic_distribution.
    import peptacular.isotope as iso

    assert not hasattr(iso, "isotopic_distribution")
    assert "isotopic_distribution" not in iso.__all__


def _composition_mass(composition: Mapping[ElementInfo, int]) -> float:
    return sum(element.get_mass(monoisotopic=True) * count for element, count in composition.items())


@pytest.mark.parametrize("target", [500.0, 1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 10000.0])
def test_averagine_comp_mass_is_within_one_dalton_of_target(target):
    """Rounding each element on its own drifted up to 21 Da from the target (issue #87)."""

    composition = pt.averagine_comp(target)
    assert abs(_composition_mass(composition) - target) < 1.0


@pytest.mark.parametrize("target", [500.0, 1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 10000.0])
def test_averagine_comp_stays_near_fractional_estimate(target):
    """The integer composition is still an averagine composition: non-negative and close to the estimate."""

    fractional = pt.estimate_averagine_comp(target)
    composition = pt.averagine_comp(target)
    for element, ideal in fractional.items():
        count = composition.get(element, 0)
        assert isinstance(count, int)
        assert count >= 0
        assert abs(count - ideal) <= 2


_ENVELOPE_PEPTIDES = [
    "SESPEKVAL",
    "PEPTIDEK",
    "ACDEFGHIK",
    "LLLLLLLLLK",
    "GGGGGGGGGGGGGGR",
    "MPEPTIDERPEPTIDEKLLSVAK",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAK",
    "YGGFLRSGGRGIK",
    "DDLSLGEKHLLKWDRVA",
    "SAMPLEPEPTIDESEQUENCE",
]
# Largest k=1..5 relative abundance error over _ENVELOPE_PEPTIDES on main before this change (issue #87).
_MAIN_MAX_ENVELOPE_ERROR_PERCENT = 101.78
# The same error with the #87 search in place, measured at 93.3805% before the #107 cache was
# added and rounded up here. An improvement to the search must lower it, so the test below
# bounds the error instead of pinning it, like the #87 baseline above.
_SEARCH_MAX_ENVELOPE_ERROR_PERCENT = 93.39
_ENVELOPE_PEAKS = 5


def _true_envelope(peptide: str) -> tuple[float, list[float]]:
    """Return the monoisotopic mass and the k=0..5 abundances of the peptide's true envelope.

    The BRAIN envelope of the real composition is anchored on the peak nearest the
    monoisotopic mass, so a selenium-containing peptide does not shift the reference.
    """

    mass = float(pt.mass(peptide, monoisotopic=True))
    distribution = brain_isotopic_distribution(pt.comp(peptide), max_isotopes=_ENVELOPE_PEAKS + 8, min_abundance_threshold=0.0)
    anchor = min(range(len(distribution)), key=lambda index: abs(distribution[index].mass - mass))
    assert abs(distribution[anchor].mass - mass) < 0.01
    return mass, [distribution[anchor + k].abundance for k in range(_ENVELOPE_PEAKS + 1)]


def _envelope_error_percent(peptide: str) -> float:
    """Largest relative abundance error of estimate_isotopic_distribution over k=1..5, in percent."""

    mass, truth = _true_envelope(peptide)
    estimate = estimate_isotopic_distribution(mass, max_isotopes=_ENVELOPE_PEAKS + 8, min_abundance_threshold=0.0)
    return 100.0 * max(abs(estimate[k].abundance - truth[k]) / truth[k] for k in range(1, _ENVELOPE_PEAKS + 1))


def test_averagine_envelope_worst_case_is_no_worse_than_main(capsys):
    """The largest k=1..5 relative abundance error is no worse than the rounding baseline (issue #87)."""

    errors = {peptide: _envelope_error_percent(peptide) for peptide in _ENVELOPE_PEPTIDES}
    worst = max(errors.values())
    for peptide, error in errors.items():
        print(f"{peptide:28s} worst k=1..5 error {error:8.2f}%")
    print(f"largest k=1..5 error: {worst:.2f}% (main: {_MAIN_MAX_ENVELOPE_ERROR_PERCENT:.2f}%)")
    assert worst <= _MAIN_MAX_ENVELOPE_ERROR_PERCENT


@pytest.mark.parametrize("target", [1000.0, 2000.0, 5000.0])
def test_averagine_comp_minimises_envelope_error_among_mass_feasible_compositions(target):
    """No other integer composition within two of the estimate and one dalton of the target has a closer envelope."""

    from itertools import product
    from math import ceil, floor

    from peptacular import isotope

    fractional = pt.estimate_averagine_comp(target)
    elements = list(fractional)
    reference = isotope._envelope_ratios(fractional)
    windows = [range(max(0, ceil(fractional[e] - 2)), floor(fractional[e] + 2) + 1) for e in elements]
    feasible = []
    for counts in product(*windows):
        composition = {e: n for e, n in zip(elements, counts, strict=True) if n}
        if composition and abs(_composition_mass(composition) - target) < 1.0:
            feasible.append(composition)

    chosen = dict(pt.averagine_comp(target))
    chosen_error = isotope._envelope_error(chosen, reference)
    assert all(isotope._envelope_error(composition, reference) >= chosen_error - 1e-12 for composition in feasible)


# Compositions returned by averagine_comp on main before the cache was added (issue #107).
_AVERAGINE_COMP_BEFORE_CACHE = {
    500.0: {"C": 21, "H": 37, "N": 7, "O": 7},
    1000.0: {"C": 42, "H": 71, "N": 12, "O": 16},
    1500.0: {"C": 66, "H": 105, "N": 19, "O": 19, "S": 1},
    2000.0: {"C": 88, "H": 142, "N": 24, "O": 27, "S": 1},
    3000.0: {"C": 132, "H": 209, "N": 37, "O": 41, "S": 1},
    5000.0: {"C": 220, "H": 351, "N": 61, "O": 68, "S": 2},
    10000.0: {"C": 441, "H": 697, "N": 125, "O": 133, "S": 4},
}


@pytest.mark.parametrize("target", list(_AVERAGINE_COMP_BEFORE_CACHE))
def test_averagine_comp_values_match_the_uncached_search(target):
    """Caching leaves the returned composition unchanged, on a cold and on a warm cache (issue #107)."""

    from peptacular import isotope

    expected = _AVERAGINE_COMP_BEFORE_CACHE[target]
    isotope._averagine_counts.cache_clear()
    cold = pt.averagine_comp(target)
    warm = pt.averagine_comp(target)
    assert {str(element): count for element, count in cold.items()} == expected
    assert warm == cold
    assert isotope._averagine_counts.cache_info().hits >= 1


def test_averagine_envelope_worst_case_is_no_worse_than_before_the_cache():
    """The worst k=1..5 error is no worse than before the cache, and the warm cache gives the same number (issue #107)."""

    from peptacular import isotope

    isotope._averagine_counts.cache_clear()
    cold = max(_envelope_error_percent(peptide) for peptide in _ENVELOPE_PEPTIDES)
    warm = max(_envelope_error_percent(peptide) for peptide in _ENVELOPE_PEPTIDES)
    assert warm == pytest.approx(cold, abs=1e-12)
    assert cold <= _SEARCH_MAX_ENVELOPE_ERROR_PERCENT


def test_averagine_comp_returns_a_fresh_counter_each_call():
    """Mutating a returned composition cannot change what the next caller gets (issue #107)."""

    first = pt.averagine_comp(1000.0)
    expected = dict(first)
    carbon = next(element for element in first if element.symbol == "C")

    first[carbon] += 100
    first[ELEMENT_LOOKUP["Se"]] = 3
    first.clear()

    second = pt.averagine_comp(1000.0)
    assert dict(second) == expected
    assert second is not first
    second[carbon] = 0
    assert dict(pt.averagine_comp(1000.0)) == expected


def test_estimate_isotopic_distribution_is_unaffected_by_mutating_averagine_comp():
    """An edited averagine_comp result does not leak into the isotope envelope (issue #107)."""

    before = estimate_isotopic_distribution(1000.0, max_isotopes=5)
    pt.averagine_comp(1000.0).clear()
    assert estimate_isotopic_distribution(1000.0, max_isotopes=5) == before


def test_averagine_comp_searches_each_mass_once():
    """Repeated calls with one mass search once, and estimate_isotopic_distribution shares that search (issue #107)."""

    from peptacular import isotope

    isotope._averagine_counts.cache_clear()
    for _ in range(10):
        pt.averagine_comp(1000.0)
    info = isotope._averagine_counts.cache_info()
    assert (info.misses, info.hits) == (1, 9)

    for _ in range(5):
        estimate_isotopic_distribution(1000.0, max_isotopes=5)
    info = isotope._averagine_counts.cache_info()
    assert (info.misses, info.hits) == (1, 14)

    pt.averagine_comp(1000.5)
    assert isotope._averagine_counts.cache_info().misses == 2

    # An int key reaches the same entry as the float. Value equality would hold with no cache
    # at all, so what names the shared entry is the miss count, which must not move.
    before = isotope._averagine_counts.cache_info()
    assert pt.averagine_comp(1000) == pt.averagine_comp(1000.0)
    after = isotope._averagine_counts.cache_info()
    assert (after.misses, after.hits) == (before.misses, before.hits + 2)


def test_averagine_comp_cache_is_bounded():
    """The cache keeps at most _AVERAGINE_CACHE_SIZE masses, like the other bounded caches in the module (issue #107)."""

    from peptacular import isotope

    isotope._averagine_counts.cache_clear()
    assert isotope._AVERAGINE_CACHE_SIZE == 1024
    assert isotope._averagine_counts.cache_info().maxsize == isotope._AVERAGINE_CACHE_SIZE
    for index in range(20):
        pt.averagine_comp(1000.0 + index)
    assert isotope._averagine_counts.cache_info().currsize == 20


def test_averagine_comp_rejects_invalid_mass_every_time():
    """A refused mass is not cached as a result (issue #107)."""

    for _ in range(2):
        with pytest.raises(pt.PeptacularError):
            pt.averagine_comp(-1.0)


# --- small masses through the averagine search (#96 follow-up) -------------------------
# A composition too small to show six isotope peaks gave _envelope_ratios fewer ratios than
# the fractional reference, and the strict zip comparing them raised ValueError. 74 of the
# 128 masses on the quarter-dalton grid below reached that comparison. v5.1.0 returned a
# composition for every one of them.
_SMALL_MASS_GRID = [round(step * 0.25, 4) for step in range(1, 129)]


def test_water_sized_mass_returns_the_composition_it_did_in_5_1_0():
    comp = pt.averagine_comp(18.0106)
    assert {element.symbol: count for element, count in comp.items()} == {"H": 2, "O": 1}


def test_estimate_isotopic_distribution_works_for_a_water_sized_mass():
    distribution = pt.estimate_isotopic_distribution(18.0106)
    assert len(distribution) == 3


@pytest.mark.parametrize("mass", _SMALL_MASS_GRID)
def test_every_small_mass_on_the_quarter_dalton_grid_returns_a_composition(mass: float):
    comp = pt.averagine_comp(mass)
    assert all(count > 0 for count in comp.values())
    assert all(isinstance(count, int) for count in comp.values())


def test_the_envelope_term_still_discriminates_between_feasible_candidates():
    # Non-vacuity for the padding: the error is still a number that orders candidates, and a
    # composition that cannot produce the envelope's peaks is still worse than one that can.
    from peptacular import isotope

    reference = isotope._envelope_ratios(isotope.estimate_averagine_comp(1000.0))
    assert len(reference) == isotope._AVERAGINE_ENVELOPE_PEAKS - 1
    good = dict(isotope._averagine_counts(1000.0))
    worse = dict(good)
    worse[next(iter(worse))] = 1
    assert isotope._envelope_error(good, reference) < isotope._envelope_error(worse, reference)

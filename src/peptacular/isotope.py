"""Aggregated isotope envelopes calculated with the BRAIN recurrence.

The implementation uses Newton-Girard power-series identities to calculate
one peak per nominal neutron offset. It also calculates the exact
probability-weighted center mass of each aggregated peak. Isotope fine
structure is intentionally outside the scope of this module.

References
----------
Dittwald et al. (2014), BRAIN 2.0, doi:10.1007/s13361-013-0796-5.
"""

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from functools import lru_cache
from math import ceil, floor, isfinite, sqrt
from typing import Final, cast

from tacular import ELEMENT_LOOKUP, FRAGMENT_ION_LOOKUP, ElementInfo, FragmentIonInfo, IonType, IonTypeLiteral

from . import constants

CARBON = ELEMENT_LOOKUP["C"]
HYDROGEN = ELEMENT_LOOKUP["H"]
NITROGEN = ELEMENT_LOOKUP["N"]
OXYGEN = ELEMENT_LOOKUP["O"]
SULFUR = ELEMENT_LOOKUP["S"]

AVERAGINE_RATIOS: Final[dict[ElementInfo, float]] = {
    CARBON: 0.044179,
    HYDROGEN: 0.069749,
    NITROGEN: 0.012344,
    OXYGEN: 0.013352,
    SULFUR: 0.0004,
}
DEFAULT_MIN_RELATIVE_ABUNDANCE: Final[float] = 0.001
MAX_ADAPTIVE_ISOTOPES: Final[int] = 4096


@dataclass(frozen=True, slots=True)
class IsotopicData:
    """One aggregated nominal isotope peak.

    ``mass`` is the exact probability-weighted center mass,
    ``neutron_count`` is the nominal offset from the all-light composition,
    and ``abundance`` is relative to the most abundant returned peak.
    """

    mass: float
    neutron_count: int
    abundance: float


type ElementPattern = tuple[tuple[int, float, float], ...]
type CompositionSignature = tuple[tuple[str, int], ...]


def estimate_averagine_comp(neutral_mass: float, ion_type: str | IonType | IonTypeLiteral = "p") -> Mapping[ElementInfo, float]:
    """Estimate an elemental composition from molecular mass.

    The fragment ion composition is treated as a fixed component. Its mass is
    removed before the averagine ratios are applied, which avoids counting the
    terminal composition twice.
    """

    mass = float(neutral_mass)
    if not isfinite(mass) or mass < 0.0:
        raise ValueError(f"neutral_mass must be finite and non-negative, got {neutral_mass!r}")

    ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]
    fixed_mass = sum(element.get_mass(monoisotopic=True) * count for element, count in ion_info.composition.items())
    scalable_mass = max(0.0, mass - fixed_mass)
    composition: dict[ElementInfo, float] = {element: ratio * scalable_mass for element, ratio in AVERAGINE_RATIOS.items()}

    for element, count in ion_info.composition.items():
        composition[element] = composition.get(element, 0.0) + count

    return composition


def averagine_comp(neutral_mass: float) -> Counter[ElementInfo]:
    """Return the nearest integer peptide averagine composition."""

    composition = Counter({element: int(floor(count + 0.5)) for element, count in estimate_averagine_comp(neutral_mass).items()})
    return Counter({element: count for element, count in composition.items() if count != 0})


@lru_cache(maxsize=256)
def _element_pattern(element: str) -> ElementPattern:
    """Return normalized nominal offsets, probabilities, and exact masses."""

    element_info = ELEMENT_LOOKUP[element]
    if element_info.mass_number is not None:
        return ((0, 1.0, element_info.get_mass(monoisotopic=True)),)

    offsets = ELEMENT_LOOKUP.get_neutron_offsets_and_abundances(element_info)
    masses = ELEMENT_LOOKUP.get_masses_and_abundances(element_info)
    entries = [
        (int(offset), float(abundance), float(mass))
        for (offset, abundance), (mass, mass_abundance) in zip(offsets, masses, strict=True)
        if abundance > 0.0 and mass_abundance > 0.0
    ]
    if not entries:
        raise ValueError(f"no positive natural isotope abundances are available for {element!r}")

    total = sum(abundance for _, abundance, _ in entries)
    minimum_offset = min(offset for offset, _, _ in entries)
    return tuple((offset - minimum_offset, abundance / total, mass) for offset, abundance, mass in entries)


def _canonical_composition(
    chemical_formula: Mapping[str | ElementInfo, int | float],
) -> tuple[CompositionSignature, float]:
    """Canonicalize integer counts and retain the mass effect of rounding."""

    counts: Counter[str] = Counter()
    mass_correction = 0.0
    for raw_element, raw_count in chemical_formula.items():
        element = str(raw_element)
        if isinstance(raw_count, bool) or not isinstance(raw_count, (int, float)):
            raise ValueError(f"element counts must be numeric and non-negative, got {raw_count!r} for {element}")
        count = float(raw_count)
        if not isfinite(count) or count < 0.0:
            raise ValueError(f"element counts must be finite and non-negative, got {raw_count!r} for {element}")
        integer_count = int(floor(count + 0.5))
        if integer_count:
            counts[element] += integer_count
        mass_correction += (count - integer_count) * ELEMENT_LOOKUP[element].get_mass(monoisotopic=True)
    return tuple(sorted(counts.items())), mass_correction


@lru_cache(maxsize=512)
def _element_recurrence_coefficients(element: str, length: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Compile one element polynomial for abundance and center-mass recurrences."""

    pattern = _element_pattern(element)
    base_probability = pattern[0][1]
    if pattern[0][0] != 0 or base_probability <= 0.0:
        raise ValueError(f"the lightest isotope of {element!r} must have positive abundance")

    abundance_ratio = [0.0] * length
    mass_ratio = [0.0] * length
    abundance_ratio[0] = 1.0
    for offset, probability, mass in pattern:
        if offset < length:
            ratio = probability / base_probability
            if offset > 0:
                abundance_ratio[offset] += ratio
            mass_ratio[offset] += ratio * mass

    log_coefficients = [0.0] * length
    for n in range(1, length):
        correction = sum(k * log_coefficients[k] * abundance_ratio[n - k] for k in range(1, n))
        log_coefficients[n] = abundance_ratio[n] - correction / n

    mass_quotient = [0.0] * length
    for n in range(length):
        mass_quotient[n] = mass_ratio[n] - sum(abundance_ratio[k] * mass_quotient[n - k] for k in range(1, n + 1))

    return tuple(log_coefficients), tuple(mass_quotient)


@lru_cache(maxsize=4096)
def _brain_coefficients(composition: CompositionSignature, length: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Calculate relative probabilities and exact center masses."""

    aggregate_log = [0.0] * length
    aggregate_mass = [0.0] * length
    for element, count in composition:
        log_coefficients, mass_quotient = _element_recurrence_coefficients(element, length)
        for index in range(length):
            aggregate_log[index] += count * log_coefficients[index]
            aggregate_mass[index] += count * mass_quotient[index]

    distribution = [0.0] * length
    distribution[0] = 1.0
    for n in range(1, length):
        scale = max(abs(value) for value in distribution[:n])
        if scale > 1e100 or 0.0 < scale < 1e-100:
            for index in range(n):
                distribution[index] /= scale
        value = sum(k * aggregate_log[k] * distribution[n - k] for k in range(1, n + 1)) / n
        if not isfinite(value):
            raise ArithmeticError("BRAIN isotope recurrence did not produce a finite distribution")
        distribution[n] = max(0.0, value)

    center_masses: list[float] = []
    for n, probability in enumerate(distribution):
        if probability == 0.0:
            center_masses.append(0.0)
            continue
        weighted_mass = sum(aggregate_mass[k] * distribution[n - k] for k in range(n + 1))
        center_masses.append(weighted_mass / probability)

    return tuple(distribution), tuple(center_masses)


def _distribution_moments(composition: CompositionSignature) -> tuple[float, float, int, int]:
    mean = 0.0
    variance = 0.0
    maximum_offset = 0
    theoretical_maximum = 0
    for element, count in composition:
        pattern = _element_pattern(element)
        element_mean = sum(offset * probability for offset, probability, _ in pattern)
        element_second = sum(offset * offset * probability for offset, probability, _ in pattern)
        mean += count * element_mean
        variance += count * max(0.0, element_second - element_mean * element_mean)
        element_maximum = max(offset for offset, _, _ in pattern)
        maximum_offset = max(maximum_offset, element_maximum)
        theoretical_maximum += count * element_maximum
    return mean, variance, maximum_offset, theoretical_maximum


def _adaptive_length(
    composition: CompositionSignature,
    min_abundance_threshold: float,
    max_isotopes: int | None,
) -> int:
    if max_isotopes is not None:
        if isinstance(max_isotopes, bool) or not isinstance(max_isotopes, int) or max_isotopes < 1:
            raise ValueError(f"max_isotopes must be positive or None, got {max_isotopes!r}")
        theoretical_maximum = _distribution_moments(composition)[3]
        return min(max_isotopes, theoretical_maximum + 1)

    mean, variance, maximum_offset, theoretical_maximum = _distribution_moments(composition)
    if min_abundance_threshold == 0.0:
        length = theoretical_maximum + 1
        if length > MAX_ADAPTIVE_ISOTOPES:
            raise ValueError("max_isotopes is required when requesting a zero abundance threshold for a large composition")
        return max(1, length)

    length = max(8, ceil(mean + 8.0 * sqrt(variance) + 2 * maximum_offset + 4))
    length = max(1, min(length, theoretical_maximum + 1))
    trailing_window = max(2, maximum_offset + 1)
    while True:
        probabilities, _ = _brain_coefficients(composition, length)
        maximum = max(probabilities)
        relative = [probability / maximum for probability in probabilities]
        window = relative[-min(trailing_window, length) :]
        apex = max(enumerate(relative), key=lambda item: item[1])[0]
        if apex < length - 1 and all(value < min_abundance_threshold for value in window):
            return length
        if length >= theoretical_maximum + 1:
            return length
        if length >= MAX_ADAPTIVE_ISOTOPES:
            raise ArithmeticError(f"adaptive isotope envelope exceeded the {MAX_ADAPTIVE_ISOTOPES}-peak safety limit")
        length = min(MAX_ADAPTIVE_ISOTOPES, theoretical_maximum + 1, length * 2)


def brain_isotopic_distribution(
    chemical_formula: Mapping[str | ElementInfo, int | float],
    max_isotopes: int | None = None,
    min_abundance_threshold: float = DEFAULT_MIN_RELATIVE_ABUNDANCE,
    charge_state: int | None = None,
) -> list[IsotopicData]:
    """Calculate an aggregated nominal isotope distribution with BRAIN.

    The returned envelope begins at the all-light composition and continues
    through the last peak at or above ``min_abundance_threshold`` relative to
    the apex. Leading peaks are retained even when they are below the threshold.
    ``max_isotopes`` limits the nominal offset window when supplied.
    """

    threshold = float(min_abundance_threshold)
    if isinstance(min_abundance_threshold, bool) or not isfinite(threshold) or not 0.0 <= threshold <= 1.0:
        raise ValueError(f"min_abundance_threshold must be in [0, 1], got {min_abundance_threshold!r}")
    if charge_state is not None and (isinstance(charge_state, bool) or not isinstance(charge_state, int)):
        raise ValueError(f"charge_state must be an integer or None, got {charge_state!r}")

    composition, mass_correction = _canonical_composition(chemical_formula)
    length = _adaptive_length(composition, threshold, max_isotopes)
    probabilities, center_masses = _brain_coefficients(composition, length)
    maximum = max(probabilities)
    relative = [probability / maximum for probability in probabilities]

    if threshold == 0.0:
        end = length
    else:
        significant = [index for index, abundance in enumerate(relative) if abundance >= threshold]
        end = significant[-1] + 1 if significant else relative.index(max(relative)) + 1

    particle_mass_offset = 0.0 if charge_state in (None, 0) else -charge_state * constants.ELECTRON_MASS
    return [
        IsotopicData(
            mass=center_masses[index] + mass_correction + particle_mass_offset,
            neutron_count=index,
            abundance=relative[index],
        )
        for index in range(end)
        if probabilities[index] > 0.0
    ]


def isotopic_distribution(
    chemical_formula: Mapping[str | ElementInfo, int | float],
    max_isotopes: int | None = None,
    min_abundance_threshold: float = DEFAULT_MIN_RELATIVE_ABUNDANCE,
    charge_state: int | None = None,
) -> list[IsotopicData]:
    """Calculate an aggregated nominal isotope distribution with BRAIN."""

    return brain_isotopic_distribution(
        chemical_formula,
        max_isotopes=max_isotopes,
        min_abundance_threshold=min_abundance_threshold,
        charge_state=charge_state,
    )


def estimate_isotopic_distribution(
    neutral_mass: float,
    max_isotopes: int | None = None,
    min_abundance_threshold: float = DEFAULT_MIN_RELATIVE_ABUNDANCE,
) -> list[IsotopicData]:
    """Estimate an aggregated peptide isotope envelope with averagine."""

    mass = float(neutral_mass)
    if not isfinite(mass) or mass < 0.0:
        raise ValueError(f"neutral_mass must be finite and non-negative, got {neutral_mass!r}")
    composition = averagine_comp(mass)
    distribution = isotopic_distribution(
        cast(Mapping[str | ElementInfo, int | float], composition),
        max_isotopes=max_isotopes,
        min_abundance_threshold=min_abundance_threshold,
    )
    if not distribution:
        return []
    shift = mass - distribution[0].mass
    return [IsotopicData(mass=peak.mass + shift, neutron_count=peak.neutron_count, abundance=peak.abundance) for peak in distribution]


def merge_isotopic_distributions(*distributions: list[IsotopicData], merge_precision: int | None = None) -> list[IsotopicData]:
    """Merge isotope distributions by summing abundances at each mass."""

    merged: dict[float, tuple[float, int]] = {}
    for distribution in distributions:
        for isotope in distribution:
            mass = round(isotope.mass, merge_precision) if merge_precision is not None else isotope.mass
            if mass in merged:
                merged[mass] = (merged[mass][0] + isotope.abundance, merged[mass][1])
            else:
                merged[mass] = (isotope.abundance, isotope.neutron_count)
    return [IsotopicData(mass=mass, neutron_count=neutron_count, abundance=abundance) for mass, (abundance, neutron_count) in sorted(merged.items())]


__all__ = [
    "AVERAGINE_RATIOS",
    "IsotopicData",
    "averagine_comp",
    "brain_isotopic_distribution",
    "estimate_averagine_comp",
    "estimate_isotopic_distribution",
    "isotopic_distribution",
    "merge_isotopic_distributions",
]

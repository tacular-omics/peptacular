from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final, cast

from tacular import ELEMENT_LOOKUP, FRAGMENT_ION_LOOKUP, ElementInfo, FragmentIonInfo, IonType, IonTypeLiteral

from . import constants

CARBON = ELEMENT_LOOKUP["C"]
HYDROGEN = ELEMENT_LOOKUP["H"]
NITROGEN = ELEMENT_LOOKUP["N"]
OXYGEN = ELEMENT_LOOKUP["O"]
SULFUR = ELEMENT_LOOKUP["S"]

AVERAGINE_RATIOS: Final[dict[ElementInfo, float]] = {CARBON: 0.044179, HYDROGEN: 0.069749, NITROGEN: 0.012344, OXYGEN: 0.013352, SULFUR: 0.0004}


def _chem_mass(
    formula: Mapping[str, int | float],
    monoisotopic: bool = True,
) -> float:
    m: float = 0.0
    for element, count in formula.items():
        m += ELEMENT_LOOKUP[element].get_mass(monoisotopic=monoisotopic) * count
    return m


@dataclass(frozen=True, slots=True)
class IsotopicData:
    mass: float
    neutron_count: int
    abundance: float


def estimate_averagine_comp(neutral_mass: float, ion_type: str | IonType | IonTypeLiteral = "p") -> Mapping[ElementInfo, float]:
    """
    Estimate elemental composition from molecular mass using the averagine model.

    .. code-block:: python

        # Example usage
        >>> round(estimate_averagine_comp(1000)['C'], 2)
        44.18

    """

    composition: dict[ElementInfo, float] = {element: ratio * neutral_mass for element, ratio in AVERAGINE_RATIOS.items()}
    ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]

    for elem_info, count in ion_info.composition.items():
        if elem_info in composition:
            composition[elem_info] += count
        else:
            composition[elem_info] = count

    return composition


def averagine_comp(neutral_mass: float) -> Counter[ElementInfo]:
    comp: Mapping[ElementInfo, int | float] = estimate_averagine_comp(neutral_mass)
    composition: Counter[ElementInfo] = Counter()
    for element, count in comp.items():
        composition[element] = int(round(count))

    # pop zeros
    for elem in list(composition.keys()):
        if composition[elem] == 0:
            del composition[elem]

    return composition


def isotopic_distribution(
    chemical_formula: Mapping[str | ElementInfo, int | float],
    max_isotopes: int | None = 10,
    min_abundance_threshold: float = 0.001,  # based on the most abundant peak
    distribution_resolution: int | None = 5,
    use_neutron_count: bool = False,
    conv_min_abundance_threshold: float = 1e-14,
    charge_state: int | None = None,
) -> list[IsotopicData]:
    """
    Calculate the isotopic distribution for a given chemical formula.

    :param chemical_formula: Chemical formula with element counts. Non-integer counts will be rounded.
    :type chemical_formula: Dict[str, Union[int, float]]
    :param max_isotopes: Maximum number of isotopes to track during convolution. Limits memory usage
                         but may affect accuracy. Use `min_abundance_threshold` to control final output.
    :type max_isotopes: Optional[int]
    :param min_abundance_threshold: Minimum relative abundance threshold for returned isotopes.
    :type min_abundance_threshold: float
    :param distribution_resolution: Decimal places for rounding masses during convolution. Simulates
                                    mass spectrometer resolution.
    :type distribution_resolution: Optional[int]
    :param use_neutron_count: If True, use neutron offsets instead of absolute masses.
    :type use_neutron_count: bool
    :param conv_min_abundance_threshold: Minimum abundance threshold during convolution. Works with
                                         `max_isotopes` to limit intermediate isotopes tracked.
    :type conv_min_abundance_threshold: float

    :return: Isotopic distribution normalized so the most abundant peak = 1.0.
    :rtype: List[IsotopicData]

    .. code-block:: python

        # Example usage
        >>> formula = {'C': 12, 'H': 6, 'N': 3}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(192.06, 1.0), (193.05, 0.01), (193.06, 0.13)]

        # Example usage with Isotopes Specified, '13C' is assumed to be
        # static and will not be used in convolution, though mass is corrected
        >>> formula = {'C': 12, 'H': 6, 'N': 3, '13C': 1}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(205.06, 1.0), (206.06, 0.01), (206.06, 0.13)]

        >>> formula = {'C': 12, 'H': 6, 'N': 3, 'Li': 0}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(192.06, 1.0), (193.05, 0.01), (193.06, 0.13)]

        >>> formula = {'C': 100, 'H': 100, 'N': 100, 'O': 100}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(4300.58, 0.92), (4301.58, 1.0), (4302.59, 0.54)]

        # Example usage
        >>> formula = {'C': 100, 'H': 150, 'N': 15, 'O': 20}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(1881.12, 0.92), (1882.12, 1.0), (1883.12, 0.54)]

        # Example usage with neutron count
        >>> formula = {'C': 100, 'H': 150, 'N': 15, 'O': 20}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5, True)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(0.0, 0.86), (1.0, 1.0), (2.0, 0.61)]

        # Example usage with negative values
        >>> formula = {'C': 100, 'Li': 10, 'N': 15, 'O': 20}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5, True)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(-1.0, 0.54), (0.0, 1.0), (1.0, 0.81)]

        # Example usage with floating values
        >>> formula = {'C': 100, 'H': 150.5, 'N': 15, 'O': 20}
        >>> result = isotopic_distribution(formula, 3, 0.0, 5, False)
        >>> [(round(r.mass, 2), round(r.abundance, 2)) for r in result]
        [(1881.62, 0.92), (1882.63, 1.0), (1883.63, 0.54)]

        # Example usage with floating values
        >>> formula = {'C': -100, 'H': 150.5, 'N': 15, 'O': 20}
        >>> isotopic_distribution(formula, 3, 0.0, 5, False)
        Traceback (most recent call last):
        ValueError: Negative values are not allowed in the chemical formula.

    """

    composition: dict[str, float | int] = {}
    for key, value in chemical_formula.items():
        if isinstance(key, str):
            composition[key] = value
        else:
            composition[str(key)] = value

    # Just use these to correct mass
    if charge_state == 0 or charge_state is None:
        particle_mass_offset = 0.0
    else:
        particle_mass_offset = -charge_state * constants.ELECTRON_MASS

    # check if any values are negative:
    if any(v < 0 for v in composition.values()):
        raise ValueError("Negative values are not allowed in the chemical formula.")

    # Remove 0 values from the chemical formula
    for elem in list(composition.keys()):
        if composition[elem] == 0:
            del composition[elem]

    delta_mass = 0.0
    if not all(isinstance(v, int) for v in composition.values()):
        prior_mass = _chem_mass(composition)
        composition = {k: int(round(v)) for k, v in composition.items()}
        post_mass = _chem_mass(composition)
        delta_mass = prior_mass - post_mass

    total_distribution = {0.0: (1.0, 0)}  # Start with a base distribution (mass/offset: (abundance, neutron_count))
    for element, count in composition.items():
        elemental_distribution = _calculate_elemental_distribution(
            element,
            int(count),
            use_neutron_count,
            conv_min_abundance_threshold,
            max_isotopes,
        )
        total_distribution = _convolve_distributions(
            total_distribution,
            elemental_distribution,
            max_isotopes,
            conv_min_abundance_threshold,
            distribution_resolution,
        )

    # Normalize abundances before filtering to ensure consistent abundance calculation
    abundances = [abundance for abundance, _ in total_distribution.values()]
    if not abundances:
        return []
    max_abundance = max(abundances)
    normalized_distribution = [
        (mass, abundance / max_abundance, neutron_count)
        for mass, (abundance, neutron_count) in sorted(total_distribution.items())
        if abundance / max_abundance >= min_abundance_threshold
    ]

    if not use_neutron_count and (delta_mass != 0.0 or particle_mass_offset != 0.0):
        normalized_distribution = [
            (mass + delta_mass + particle_mass_offset, abundance, neutron_count) for mass, abundance, neutron_count in normalized_distribution
        ]

    return [IsotopicData(mass=mass, neutron_count=neutron_count, abundance=abundance) for mass, abundance, neutron_count in normalized_distribution]


def merge_isotopic_distributions(*distributions: list[IsotopicData], merge_precision: int | None = None) -> list[IsotopicData]:
    """
    Merge multiple isotopic distributions by summing abundances at each mass.

    :param distributions: Isotopic distributions to merge.
    :type distributions: List[IsotopicData]
    :param merge_precision: Decimal places for rounding masses during merge.
    :type merge_precision: Optional[int]

    :return: Merged isotopic distribution.
    :rtype: List[IsotopicData]

    .. code-block:: python

        # Example usage
        >>> d1 = [IsotopicData(1.0, 0, 0.5), IsotopicData(2.0, 1, 0.5)]
        >>> d2 = [IsotopicData(1.0, 0, 0.5), IsotopicData(2.0, 1, 0.5)]
        >>> result = merge_isotopic_distributions(d1, d2)
        >>> [(r.mass, r.abundance) for r in result]
        [(1.0, 1.0), (2.0, 1.0)]

    """
    merged_distribution: dict[float, tuple[float, int]] = {}
    for distribution in distributions:
        for isotope in distribution:
            mass = isotope.mass
            if merge_precision is not None:
                mass = round(mass, merge_precision)
            # Merge the distributions, keep neutron count from first occurrence
            if mass in merged_distribution:
                merged_distribution[mass] = (
                    merged_distribution[mass][0] + isotope.abundance,
                    merged_distribution[mass][1],  # Keep existing neutron count
                )
            else:
                merged_distribution[mass] = (isotope.abundance, isotope.neutron_count)

    return [
        IsotopicData(mass=mass, neutron_count=neutron_count, abundance=abundance)
        for mass, (abundance, neutron_count) in sorted(merged_distribution.items(), key=lambda x: x[0])
    ]


def estimate_isotopic_distribution(
    neutral_mass: float,
    max_isotopes: int | None = 10,
    min_abundance_threshold: float = 0.001,
    distribution_resolution: int | None = 5,
    use_neutron_count: bool = False,
    conv_min_abundance_threshold: float = 1e-14,
) -> list[IsotopicData]:
    """
    Estimate isotopic distribution from molecular mass using the averagine model.

    :param neutral_mass: Neutral mass of the molecule.
    :type neutral_mass: float
    :param max_isotopes: Maximum number of isotopes to track during convolution.
    :type max_isotopes: int
    :param min_abundance_threshold: Minimum relative abundance threshold for returned isotopes.
    :type min_abundance_threshold: float
    :param distribution_resolution: Decimal places for rounding masses during convolution.
    :type distribution_resolution: Optional[int]
    :param use_neutron_count: If True, use neutron offsets instead of absolute masses.
    :type use_neutron_count: bool
    :param conv_min_abundance_threshold: Minimum abundance threshold during convolution.
    :type conv_min_abundance_threshold: float

    :return: Predicted isotopic distribution normalized so the most abundant peak = 1.0.
    :rtype: List[IsotopicData]

    .. code-block:: python

        # Example usage: the monoisotopic peak is anchored to the requested mass
        >>> result = estimate_isotopic_distribution(800, 3, 0.0, 5)
        >>> [(round(r.mass, 3), round(r.abundance, 3)) for r in result]
        [(800.0, 1.0), (801.003, 0.379), (802.007, 0.07)]

        # Example usage with neutron count
        >>> result = estimate_isotopic_distribution(800, 3, 0.0, 5, True)
        >>> [(round(r.mass, 3), round(r.abundance, 3)) for r in result]
        [(0.0, 1.0), (1.0, 0.426), (2.0, 0.113)]

    """
    # Calculate the total number of each atom in the molecule based on its molecular mass
    total_atoms: Counter[ElementInfo] = averagine_comp(neutral_mass)

    distribution = isotopic_distribution(
        cast(Mapping[str | ElementInfo, int | float], total_atoms),
        max_isotopes,
        min_abundance_threshold,
        distribution_resolution,
        use_neutron_count,
        conv_min_abundance_threshold,
    )

    # Averagine only supplies the rough *shape* (peak spacings and abundances). Anchor the
    # envelope so the monoisotopic peak lands on the requested mass, otherwise the absolute
    # masses drift by the averagine composition's own monoisotopic-mass offset. When
    # ``use_neutron_count`` is set the masses are neutron offsets, so no anchoring is applied.
    if not use_neutron_count and distribution:
        averagine_mono_mass = sum(element.get_mass(monoisotopic=True) * count for element, count in total_atoms.items())
        shift = neutral_mass - averagine_mono_mass
        distribution = [
            IsotopicData(mass=d.mass + shift, neutron_count=d.neutron_count, abundance=d.abundance) for d in distribution
        ]

    return distribution


def _convolve_distributions(
    dist1: dict[float, tuple[float, int]],
    dist2: dict[float, tuple[float, int]],
    max_isotopes: int | None,
    min_abundance_threshold: float,
    distribution_resolution: int | None,
) -> dict[float, tuple[float, int]]:
    """
    Convolve two distributions to calculate the distribution of their sum.

    :param dist1: First distribution mapping mass/offset to (abundance, neutron_count).
    :type dist1: Dict[float, Tuple[float, int]]
    :param dist2: Second distribution mapping mass/offset to (abundance, neutron_count).
    :type dist2: Dict[float, Tuple[float, int]]
    :param max_isotopes: Maximum number of isotopes to retain based on abundance.
    :type max_isotopes: Optional[int]
    :param min_abundance_threshold: Minimum abundance threshold for keeping isotopes.
    :type min_abundance_threshold: float
    :param distribution_resolution: Decimal places for rounding masses.
    :type distribution_resolution: Optional[int]

    :return: Convolved isotopic distribution as mass-to-(abundance, neutron_count) mapping.
    :rtype: Dict[float, Tuple[float, int]]

    .. code-block:: python

        # Example usage
        >>> d1 = {1.0: (0.5, 0), 2.0: (0.5, 1)}
        >>> d2 = {1.0: (0.5, 0), 2.0: (0.5, 1)}
        >>> result = _convolve_distributions(d1, d2, 3, 0.0, 5)
        >>> sorted(result.items())
        [(2.0, (0.25, 0)), (3.0, (0.5, 1)), (4.0, (0.25, 2))]

    """

    result: dict[float, tuple[float, int]] = {}
    for mass1, (abundance1, neutron1) in dist1.items():
        for mass2, (abundance2, neutron2) in dist2.items():
            new_abundance = abundance1 * abundance2
            if new_abundance < min_abundance_threshold:
                # ``dist2`` is not guaranteed to be sorted by descending abundance, so
                # skip only this low-abundance pairing rather than aborting the whole
                # inner loop (which could drop later high-abundance pairings).
                continue

            new_mass = mass1 + mass2
            if distribution_resolution is not None:
                new_mass = round(new_mass, distribution_resolution)

            new_neutron_count = neutron1 + neutron2
            if new_mass in result:
                result[new_mass] = (
                    result[new_mass][0] + new_abundance,
                    new_neutron_count,
                )
            else:
                result[new_mass] = (new_abundance, new_neutron_count)

    # Apply max isotopes limit and sort by abundance
    sorted_result = sorted(result.items(), key=lambda x: x[1][0], reverse=True)
    if max_isotopes is not None:
        # Retain only the top `max_isotopes` isotopes based on abundance
        return dict(sorted_result[:max_isotopes])

    return result


def _calculate_elemental_distribution_slow(
    element: str,
    count: int,
    use_neutron_count: bool,
    min_abundance_threshold: float = 10e-9,
) -> dict[float, tuple[float, int]]:
    """
    Calculate the isotopic distribution for an element.

    :param element: The element.
    :type element: str
    :param count: The number of atoms.
    :type count: int
    :param use_neutron_count: Whether to use neutron offsets instead of masses.
    :type use_neutron_count: bool

    :return: The isotopic distribution mapping mass/offset to (abundance, neutron_count).
    :rtype: Dict[float, Tuple[float, int]]

    .. code-block:: python

        # Example usage
        >>> _calculate_elemental_distribution_slow('C', 2, False)
        {24.0: (0.9787144899999999, 0), 25.00335483507: (0.02117102, 1), 26.00670967014: (0.00011448999999999998, 2)}

        # Example using neutron count
        >>> _calculate_elemental_distribution_slow('C', 2, True)
        {0.0: (0.9787144899999999, 0), 1.0: (0.02117102, 1), 2.0: (0.00011448999999999998, 2)}

    """

    if use_neutron_count is True:
        isotopes = ELEMENT_LOOKUP.get_neutron_offsets_and_abundances(element)
    else:
        isotopes = ELEMENT_LOOKUP.get_masses_and_abundances(element)

    # Start with a distribution for an element not present (mass=0, abundance=1, neutron_count=0)
    distribution = {0.0: (1.0, 0)}
    for _ in range(count):
        # Update the distribution by convolving it with the isotopes' distribution each time
        # Convert isotopes list to dict with neutron count tracking
        isotope_distribution: dict[float, tuple[float, int]] = {}
        for i, (mass_or_offset, abundance) in enumerate(isotopes):
            # When use_neutron_count=True, mass_or_offset is already the neutron offset (int)
            # When use_neutron_count=False, mass_or_offset is mass, and we track neutron by index
            neutron_count = int(mass_or_offset) if use_neutron_count else i
            isotope_distribution[mass_or_offset] = (abundance, neutron_count)
        distribution = _convolve_distributions(distribution, isotope_distribution, None, min_abundance_threshold, None)
    return distribution


def _calculate_elemental_distribution(
    element: str | ElementInfo,
    count: int,
    use_neutron_count: bool,
    min_abundance_threshold: float = 1e-14,
    max_isotopes: int | None = None,
) -> dict[float, tuple[float, int]]:
    """Calculate elemental isotopic distribution using binary exponentiation for efficiency."""

    if not isinstance(element, ElementInfo):
        elem_info = ELEMENT_LOOKUP[element]
    else:
        elem_info = element

    if elem_info.mass_number is not None:
        # Monoisotopic elements have only one isotope
        if use_neutron_count:
            return {0.0: (1.0, 0)}
        else:
            mass = elem_info.get_mass(monoisotopic=True)
            return {mass: (1.0, 0)}

    if use_neutron_count:
        isotopes = ELEMENT_LOOKUP.get_neutron_offsets_and_abundances(element)
    else:
        isotopes = ELEMENT_LOOKUP.get_masses_and_abundances(element)

    # Build base distribution
    isotope_distribution = {}
    for i, (mass_or_offset, abundance) in enumerate(isotopes):
        neutron_count = int(mass_or_offset) if use_neutron_count else i
        isotope_distribution[mass_or_offset] = (abundance, neutron_count)

    # Fast exponentiation by squaring: compute isotope_distribution^count
    result = {0.0: (1.0, 0)}
    base = isotope_distribution

    while count > 0:
        if count % 2 == 1:
            result = _convolve_distributions(result, base, max_isotopes, min_abundance_threshold, None)
        base = _convolve_distributions(base, base, max_isotopes, min_abundance_threshold, None)
        count //= 2

    return result


class IsotopeLookup:
    def __init__(
        self,
        mass_step: int = 50,
        max_isotopes: int = 25,
        min_abundance_threshold: float = 0.005,
        use_neutron_count: bool = True,
        is_abundance_sum: bool = True,
    ):
        self.mass_step = mass_step
        self.max_isotopes = max_isotopes
        self.min_abundance_threshold = min_abundance_threshold
        self.use_neutron_count = use_neutron_count
        self.is_abundance_sum = is_abundance_sum
        self.lookup: dict[int, list[IsotopicData]] = {}

    def _get_mass_bin(self, mass: float) -> int:
        """Calculate the mass bin for a given mass."""
        return round(mass / self.mass_step) * self.mass_step

    def _generate_pattern(self, mass_bin: int) -> list[IsotopicData]:
        """Generate isotope pattern for a specific mass bin."""
        iso_pattern: list[IsotopicData] = estimate_isotopic_distribution(
            neutral_mass=mass_bin,
            max_isotopes=self.max_isotopes,
            min_abundance_threshold=self.min_abundance_threshold,
            use_neutron_count=self.use_neutron_count,
        )
        return iso_pattern

    def get_isotope_pattern(self, mass: float) -> list[IsotopicData]:
        """
        Get the isotope pattern for a given mass.

        Generates and caches the pattern if not already present.
        """
        mass_bin = self._get_mass_bin(mass)

        # Check if pattern exists, if not generate it
        if mass_bin not in self.lookup:
            self.lookup[mass_bin] = self._generate_pattern(mass_bin)

        return self.lookup[mass_bin]

    def clear_cache(self):
        """Clear the cached isotope patterns."""
        self.lookup.clear()

    def get_cache_size(self) -> int:
        """Return the number of cached isotope patterns."""
        return len(self.lookup)

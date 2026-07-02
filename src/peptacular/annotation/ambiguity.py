from __future__ import annotations

import itertools
from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

from tacular import IonType

from ..annotation.parser import Interval
from ..constants import ModType

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def condense_ambiguity_to_xnotation(annotation: ProFormaAnnotation, inplace: bool = False) -> ProFormaAnnotation:
    """
    Condense ambiguity intervals in the annotation to X[+Mod]
    """
    if not inplace:
        return condense_ambiguity_to_xnotation(annotation.copy(), inplace=True)

    elems = annotation.split()
    for elem in elems:
        if elem.has_intervals:
            # drop unknown and labile and charge / adducts
            elem_annot = elem.filter_mods([ModType.INTERNAL, ModType.STATIC, ModType.ISOTOPE, ModType.INTERVAL])
            # Use the neutral mass with no charge; charge=1 would fold an extra
            # proton into the mass baked into the condensed X[+...] residue.
            mass = elem_annot.mass(ion_type=IonType.NEUTRAL, charge=0)

            elem.sequence = "X"
            elem.clear_internal_mods()
            elem.clear_intervals()
            elem.set_internal_mods_at_index(0, mass)

    new_annot = annotation.__class__.join(elems)
    annotation.update(new_annot)
    return annotation


def annotate_ambiguity(
    annotation: ProFormaAnnotation,
    forward_coverage: list[int],
    reverse_coverage: list[int],
    mass_shift: Any | None,
    *,
    add_mods_to_intervals: bool,
    sort_mods: bool,
    inplace: bool,
) -> ProFormaAnnotation:
    """
    Generate ambiguity intervals based on the coverage of the sequence.

    Args:
        annotation: The annotation to modify
        forward_coverage: Forward coverage data
        reverse_coverage: Reverse coverage data
        mass_shift: Optional mass shift to apply
        inplace: Whether to modify annotation in place

    Returns:
        Annotation with ambiguity intervals added

    Raises:
        ValueError: If annotation already contains intervals or coverage length doesn't match
    """
    if not inplace:
        return annotate_ambiguity(
            annotation.copy(),
            forward_coverage,
            reverse_coverage,
            mass_shift,
            add_mods_to_intervals=add_mods_to_intervals,
            sort_mods=sort_mods,
            inplace=True,
        )

    if annotation.has_intervals:
        raise ValueError("Annotation should not contain intervals")

    _validate_coverage_lengths(forward_coverage, reverse_coverage, len(annotation))

    forward_intervals = _construct_ambiguity_intervals(forward_coverage, reverse=False)
    reverse_intervals = _construct_ambiguity_intervals(reverse_coverage, reverse=True)
    ambiguity_intervals = _combine_ambiguity_intervals(forward_intervals, reverse_intervals)

    intervals = [Interval(start, end + 1, True, None, validate=annotation._validate) for start, end in ambiguity_intervals]

    annotation.extend_intervals(intervals)

    if mass_shift is not None:
        _apply_mass_shift(annotation, forward_coverage, reverse_coverage, mass_shift)

    # add_mods_to_intervals
    if add_mods_to_intervals:
        annotation.condense_mods_to_intervals(inplace=True)

    if sort_mods:
        annotation.sort_mods(inplace=True)

    return annotation


def _construct_ambiguity_intervals(counts: list[int], reverse: bool) -> list[tuple[int, int]]:
    """
    Construct intervals for sequences of zeros in the counts list. When reverse is false, start from the left hand side
    and move to the right. When reverse is true, start from the right hand side and move to the left. Intervals start
    at 0 and end on any positive value. Both are inclusive. Returned intervals should be in forwards format, that is
    have a starting value less than the ending value.

    :param counts: List of integers (typically counts)
    :param reverse: If True, reverse the list before processing
    :return: List of intervals [start, end] indicating runs of zeros

    .. code-block:: python

        # [0, 1, 1, 1, 0, 0, 0]
        # [1, 1, 0, 0, 1, 1, 1] # ambiguity
        >>> _construct_ambiguity_intervals([0, 1, 1, 1, 0, 0, 0], reverse=False)
        [(0, 1), (4, 6)]

        # [0, 0, 1, 1, 1, 1, 0]
        # [1, 1, 0, 0, 0, 1, 1] # ambiguity
        >>> _construct_ambiguity_intervals([0, 0, 1, 1, 1, 0, 0], reverse=True)
        [(0, 1), (4, 6)]

        >>> _construct_ambiguity_intervals([0, 1, 1, 1, 0, 0, 1], reverse=False)
        [(0, 1), (4, 6)]
    """

    if reverse:
        ambiguity_intervals = _construct_ambiguity_intervals(counts[::-1], reverse=False)
        ambiguity_intervals = [(len(counts) - 1 - end, len(counts) - 1 - start) for start, end in ambiguity_intervals]
        # sort the intervals
        ambiguity_intervals.sort(key=lambda x: x[0])
        return ambiguity_intervals

    ambiguity_intervals: list[tuple[int, int]] = []
    current_interval: tuple[int, int] | None = None
    for i, cnt in enumerate(counts):
        if cnt == 0:
            if current_interval is not None:
                current_interval = (current_interval[0], i)
            else:
                current_interval = (i, i)

        else:
            if current_interval is not None:
                current_interval = (current_interval[0], i)
                ambiguity_intervals.append(current_interval)
                current_interval = None
            else:
                continue

    if current_interval is not None:
        current_interval = (current_interval[0], len(counts) - 1)
        ambiguity_intervals.append(current_interval)

    return ambiguity_intervals


def _combine_ambiguity_intervals(
    *interval_lists: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """
    Merge multiple lists of ambiguity intervals into a single list of common ambiguity intervals.

    This function identifies positions that are ambiguous across all provided interval lists.
    For a position to be considered ambiguous in the result, it must be contained in at least
    one interval from each input list. The function then constructs optimized intervals
    covering these common ambiguous positions.

    Intervals are represented as tuples (start, end) where:
    - start is inclusive
    - end is exclusive

    Intervals with identical start and end values (zero-length intervals) are removed.

    :param interval_lists: Variable number of lists containing ambiguity intervals
    :type interval_lists: List[Tuple[int, int]]

    :return: A list of merged intervals representing positions that are ambiguous across all input lists
    :rtype: List[Tuple[int, int]]

    .. code-block:: python

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0,1)])
        [(0, 1)]

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0,1), (4,5)])
        [(0, 1), (4, 5)]

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0, 4), (5, 6)])
        [(0, 1), (5, 6)]

        >>> _combine_ambiguity_intervals([(2, 5)], [(3, 6)])
        [(3, 5)]

        >>> _combine_ambiguity_intervals([(0, 1)], [(4, 6)])
        []
    """

    # First, collect all unique intervals from input lists
    all_intervals: set[tuple[int, int]] = set()
    for interval_list in interval_lists:
        for interval in interval_list:
            all_intervals.add(interval)

    # Remove intervals where start == end (these are not ambiguous)
    filtered_intervals = {(start, end) for start, end in all_intervals if start != end}

    # Find all possible indices that are covered by any interval
    all_indices: set[int] = set()
    for start, end in filtered_intervals:
        for i in range(start, end):
            all_indices.add(i)

    # For each index, check if it's contained in at least one interval from each input list
    common_indices: set[int] = set()
    for idx in all_indices:
        is_common = True
        for interval_list in interval_lists:
            if not any(start <= idx < end for start, end in interval_list):
                is_common = False
                break
        if is_common:
            common_indices.add(idx)

    # If no common indices found, return empty list
    if not common_indices:
        return []

    # Construct new intervals from the common indices
    result: list[tuple[int, int]] = []
    if common_indices:
        sorted_indices = sorted(common_indices)
        start = sorted_indices[0]
        for i in range(1, len(sorted_indices)):
            if sorted_indices[i] > sorted_indices[i - 1] + 1:
                # Gap found, close the current interval and start a new one
                result.append((start, sorted_indices[i - 1] + 1))
                start = sorted_indices[i]
        # Add the last interval
        result.append((start, sorted_indices[-1] + 1))

    return result


def _get_mass_shift_interval(forward_coverage: list[int], reverse_coverage: list[int]) -> tuple[int, int] | None:
    """
        Determine the interval where a mass shift should be placed based on fragment ion coverage.

    This function examines the forward and reverse ion coverage to identify the region where
    a mass shift (such as a modification) should be positioned. It returns the start and end
    indices (inclusive) of this region, or None if no suitable region is found.

    The mass shift interval is determined by:
    1. Finding the highest position with forward ion coverage
    2. Finding the lowest position with reverse ion coverage
    3. The mass shift belongs between these two positions

    :param forward_coverage: Binary list indicating forward ion coverage (1) or no coverage (0)
    :type forward_coverage: List[int]
    :param reverse_coverage: Binary list indicating reverse ion coverage (1) or no coverage (0)
    :type reverse_coverage: List[int]

    :return: A tuple containing the start and end indices (inclusive) for the mass shift,
             or None if there is no valid interval
    :rtype: Optional[Tuple[int, int]]

    .. code-block:: python

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,0,1,1,1])
        (3, 3)

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,1,1,1,1])
        (3, 3)

        >>> _get_mass_shift_interval([1,1,0,0,0,0,0], [0,0,0,0,1,1,1])
        (2, 3)

        >>> _get_mass_shift_interval([0,0,0,0,0,0,0], [0,0,0,0,1,1,1])
        (0, 3)

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,0,0,0,0])
        (3, 6)

        >>> _get_mass_shift_interval([1,1,1,1,1,0,0], [0,0,0,0,1,1,1]) # None

    """

    highest_forward_fragment = [i for i, cnt in enumerate(forward_coverage) if cnt > 0]
    if len(highest_forward_fragment) == 0:
        highest_forward_fragment = -1
    else:
        highest_forward_fragment = max(highest_forward_fragment)

    highest_reverse_fragment = [i for i, cnt in enumerate(reverse_coverage) if cnt > 0]
    if len(highest_reverse_fragment) == 0:
        highest_reverse_fragment = len(reverse_coverage)
    else:
        highest_reverse_fragment = min(highest_reverse_fragment)

    if highest_forward_fragment >= highest_reverse_fragment:
        return None

    if highest_forward_fragment == highest_reverse_fragment - 1:
        return (highest_forward_fragment + 1, highest_forward_fragment + 1)

    # if there is a gap between the two, return the gap
    return (highest_forward_fragment + 1, highest_reverse_fragment - 1)


def _validate_coverage_lengths(forward_coverage: list[int], reverse_coverage: list[int], seq_len: int) -> None:
    """Validate that coverage lengths match sequence length"""
    if len(forward_coverage) != len(reverse_coverage) != seq_len:
        raise ValueError(f"Coverage length does not match sequence length: {len(forward_coverage)} != {len(reverse_coverage)} != {seq_len}")


def _apply_mass_shift(
    annotation: ProFormaAnnotation,
    forward_coverage: list[int],
    reverse_coverage: list[int],
    mass_shift: Any,
) -> None:
    """Apply mass shift to annotation based on coverage"""

    mass_shift_interval = _get_mass_shift_interval(forward_coverage, reverse_coverage)

    if mass_shift_interval is None:
        annotation.append_labile_mod(mass_shift)
    elif mass_shift_interval[0] == mass_shift_interval[1]:
        # Add modification to specific position
        annotation.append_internal_mod_at_index(mass_shift_interval[0], mass_shift)
    else:
        # Add interval modification

        # check if interval already exists
        found_int = annotation.get_interval(mass_shift_interval[0], mass_shift_interval[1] + 1)

        if found_int is not None:
            # add mod to existing interval
            found_int.append_mod(mass_shift)
            return

        mod_interval = Interval(
            mass_shift_interval[0],
            mass_shift_interval[1] + 1,
            False,
            None,
            validate=annotation._validate,
        )
        annotation.append_interval(mod_interval)


def group_by_ambiguity(annotations: Iterable[ProFormaAnnotation], precision: int = 4) -> list[tuple[ProFormaAnnotation, ...]]:
    annotation_masses: list[tuple[ProFormaAnnotation, set[int]]] = []

    if precision < 0 or precision > 10:
        raise ValueError(f"Precision must be an integer between 0 and 10, got {precision}")

    mult = 10**precision

    for annotation in annotations:
        elements = annotation.split()
        masses = [int(round(element.mass(ion_type=IonType.NEUTRAL), precision)) * mult for element in elements]
        forward_masses = [mass for mass in itertools.accumulate(masses)]
        reverse_masses = [mass for mass in itertools.accumulate(reversed(masses))]
        unique_masses = set(forward_masses) | set(reverse_masses)
        # print(annotation.serialize(), sorted(list(map(int, unique_masses))))
        annotation_masses.append((annotation, unique_masses))

    # sort annotations by length of unique masses (descending)
    annotation_masses.sort(key=lambda x: len(x[1]), reverse=True)

    groups: list[list[ProFormaAnnotation]] = []
    added_indexes: set[int] = set()

    for i, (annotation, unique_masses) in enumerate(annotation_masses):
        if i not in added_indexes:
            # Check if this annotation is a subset of any already processed annotation
            # If so, skip it (it can't start its own group)
            is_subset_of_earlier = False
            for j in range(i):
                if j not in added_indexes:  # Only check group leaders
                    _, earlier_masses = annotation_masses[j]
                    if unique_masses.issubset(earlier_masses):
                        is_subset_of_earlier = True
                        break

            if not is_subset_of_earlier:
                # Create new group with current annotation
                current_group = [annotation]
                added_indexes.add(i)

                # Check all remaining annotations for subsets
                for j, (other_annotation, other_unique_masses) in enumerate(annotation_masses[i + 1 :], start=i + 1):
                    if j not in added_indexes and other_unique_masses.issubset(unique_masses):
                        current_group.append(other_annotation)
                        added_indexes.add(j)

                groups.append(current_group)

    return [tuple(group) for group in groups]


def unique_fragments(annotations: Iterable[ProFormaAnnotation], precision: int = 4) -> list[int]:
    annotation_masses: list[tuple[ProFormaAnnotation, set[int]]] = []

    if precision < 0 or precision > 10:
        raise ValueError(f"Precision must be an integer between 0 and 10, got {precision}")

    mult = 10**precision

    for annotation in annotations:
        elements = annotation.split()
        masses = [int(round(element.mass(ion_type=IonType.NEUTRAL), precision)) * mult for element in elements]
        forward_masses = [mass for mass in itertools.accumulate(masses)]
        reverse_masses = [mass for mass in itertools.accumulate(reversed(masses))]
        unique_masses = set(forward_masses) | set(reverse_masses)
        # print(annotation.serialize(), sorted(list(map(int, unique_masses))))
        annotation_masses.append((annotation, unique_masses))

    # Second pass: find globally unique masses for each annotation
    results: list[tuple[ProFormaAnnotation, int]] = []

    for i, (annotation, masses) in enumerate(annotation_masses):
        # Get all other masses from other annotations
        all_other_masses: set[int] = set()
        for j, (_, other_masses) in enumerate(annotation_masses):
            if i != j:
                all_other_masses.update(other_masses)

        # Find masses that are unique to this annotation (not in any other)
        globally_unique_masses = masses - all_other_masses
        results.append((annotation, len(globally_unique_masses)))

    return [count for _, count in results]

from __future__ import annotations

import itertools
from collections.abc import Generator
from typing import TYPE_CHECKING, Any

from ..constants import ModType

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation
    from .mod import Mods


def generate_permutations(annotation: ProFormaAnnotation, size: int | None = None) -> Generator[ProFormaAnnotation]:
    """
    Generate all permutations of the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to permute
        size: Number of elements in each permutation (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each permutation

    Note:
        Non-internal modifications (terminal, labile, etc.) are preserved
        and applied to each permutation.
    """
    if size is None:
        size = len(annotation)

    # Extract mods that should be preserved across permutations
    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    # Generate permutations of split amino acids
    split_aas = annotation_copy.split()
    from .annotation import ProFormaAnnotation

    for permutation in itertools.permutations(split_aas, size):
        # Create new annotation from permuted sequence
        combined_sequence = "".join(aa.serialize() for aa in permutation)
        result = ProFormaAnnotation.parse_single(combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


def generate_product(annotation: ProFormaAnnotation, repeat: int | None = None) -> Generator[ProFormaAnnotation]:
    """
    Generate the Cartesian product of the annotation sequence with itself.

    Args:
        annotation: The ProFormaAnnotation to generate products from
        repeat: Number of repetitions (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each product combination

    Note:
        Non-internal modifications are preserved and applied to each product.
    """
    if repeat is None:
        repeat = len(annotation)

    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    split_aas = annotation_copy.split()

    from .annotation import ProFormaAnnotation

    for product in itertools.product(split_aas, repeat=repeat):
        combined_sequence = "".join(aa.serialize() for aa in product)
        result = ProFormaAnnotation.parse_single(combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


def generate_combinations(annotation: ProFormaAnnotation, r: int | None = None) -> Generator[ProFormaAnnotation]:
    """
    Generate all combinations of the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to generate combinations from
        r: Length of each combination (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each combination

    Note:
        Non-internal modifications are preserved and applied to each combination.
    """
    if r is None:
        r = len(annotation)

    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    split_aas = annotation_copy.split()
    from .annotation import ProFormaAnnotation

    for combination in itertools.combinations(split_aas, r=r):
        combined_sequence = "".join(aa.serialize() for aa in combination)
        result = ProFormaAnnotation.parse_single(combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


def generate_combinations_with_replacement(annotation: ProFormaAnnotation, r: int | None = None) -> Generator[ProFormaAnnotation]:
    """
    Generate all combinations of the annotation sequence with replacement.

    Args:
        annotation: The ProFormaAnnotation to generate combinations from
        r: Length of each combination (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each combination with replacement

    Note:
        Non-internal modifications are preserved and applied to each combination.
    """
    if r is None:
        r = len(annotation)

    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    split_aas = annotation_copy.split()
    from .annotation import ProFormaAnnotation

    for combination in itertools.combinations_with_replacement(split_aas, r=r):
        combined_sequence = "".join(aa.serialize() for aa in combination)
        result = ProFormaAnnotation.parse_single(combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


# Helper functions


def _extract_non_internal_mods(
    annotation: ProFormaAnnotation,
) -> dict[ModType, Mods[Any] | int | None]:
    """Extract non-internal modifications from annotation"""

    non_internal_mod_types = [
        ModType.NTERM,
        ModType.CTERM,
        ModType.LABILE,
        ModType.ISOTOPE,
        ModType.STATIC,
        ModType.UNKNOWN,
        ModType.CHARGE,
    ]

    d: dict[ModType, Any] = {}
    d[ModType.NTERM] = annotation.nterm_mods
    d[ModType.CTERM] = annotation.cterm_mods
    d[ModType.LABILE] = annotation.labile_mods
    d[ModType.ISOTOPE] = annotation.isotope_mods
    d[ModType.STATIC] = annotation.static_mods
    d[ModType.UNKNOWN] = annotation.unknown_mods
    d[ModType.CHARGE] = annotation.charge

    # remove mods
    annotation.clear_mods(mods=non_internal_mod_types)

    return d


def _apply_extracted_mods(annotation: ProFormaAnnotation, mods_dict: dict[ModType, Mods[Any] | int | None]) -> None:
    """Apply extracted modifications back to annotation"""
    for mod_type, mod_value in mods_dict.items():
        match mod_type:
            case ModType.NTERM:
                annotation.nterm_mods = mod_value
            case ModType.CTERM:
                annotation.cterm_mods = mod_value
            case ModType.LABILE:
                annotation.labile_mods = mod_value
            case ModType.ISOTOPE:
                annotation.isotope_mods = mod_value
            case ModType.STATIC:
                annotation.static_mods = mod_value
            case ModType.UNKNOWN:
                annotation.unknown_mods = mod_value
            case ModType.CHARGE:
                annotation.charge = mod_value
            case _:
                raise ValueError(f"Unknown ModType: {mod_type}")

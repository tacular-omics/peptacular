from __future__ import annotations

import itertools
import re
import warnings
from collections.abc import Generator, Iterable, Mapping
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def get_mod_index_from_aa(peptide: str, mod_aa: str) -> set[int]:
    # Given a peptide sequence and a modification, return the indices of the modification
    sites: set[int] = set()
    for aa in mod_aa:
        sites.update(i for i, peptide_aa in enumerate(peptide) if peptide_aa == aa)
    return sites


def get_mod_index_from_regex(peptide: str, mod: str) -> set[int]:
    # given a regex, ensure that it is just a single site
    # then return a set of all indexes

    sites: set[int] = set()
    match = re.search(mod, peptide)
    if match:
        sites.add(match.start())
        if match.end() - match.start() != 0:
            warnings.warn("Using start of regex", UserWarning, stacklevel=2)
    return sites


def get_mod_index(peptide: str, mod: str | None, is_regex: bool) -> set[int]:
    if mod is None:
        # represents all sites
        return set(range(len(peptide)))

    if is_regex:
        return get_mod_index_from_regex(peptide, mod)
    return get_mod_index_from_aa(peptide, mod)


def get_sites(
    peptide: str, mods: Mapping[str | None, Iterable[Any]], is_regex: bool
) -> dict[int, list[Any]]:
    sites: dict[int, list[Any]] = {}
    for mod_pattern, mod_values in mods.items():
        mod_sites = get_mod_index(peptide, mod_pattern, is_regex)
        for site in mod_sites:
            if site not in sites:
                sites[site] = []
            if isinstance(mod_values, str):
                sites[site].extend(mod_values)
            else:
                sites[site].extend(mod_values)
    return sites


def ensure_single_static_mod(mods: dict[int, list[Any]]) -> None:
    # warn if multiple static mods are applied to the same site
    for site, mod in mods.items():
        if len(mod) > 1:
            warnings.warn(
                f"Multiple static modifications applied to site {site}: {mod}",
                UserWarning,
                stacklevel=2,
            )


def apply_mods(
    annotation: ProFormaAnnotation,
    nterm: Mapping[str, Iterable[Any]] | None = None,
    cterm: Mapping[str, Iterable[Any]] | None = None,
    internal: Mapping[str, Iterable[Any]] | None = None,
    inplace: bool = True,
    is_regex: bool = False,
) -> ProFormaAnnotation:
    if inplace is False:
        return apply_mods(
            annotation.copy(),
            nterm=nterm,
            cterm=cterm,
            internal=internal,
            inplace=True,
        )

    # get sites
    nterm_sites = get_sites(annotation.sequence, nterm, is_regex) if nterm else {}
    cterm_sites = get_sites(annotation.sequence, cterm, is_regex) if cterm else {}
    internal_sites = (
        get_sites(annotation.sequence, internal, is_regex) if internal else {}
    )

    # warn if more than 1 static mod per site
    # ensure_single_static_mod(nterm_sites)
    # ensure_single_static_mod(cterm_sites)
    # ensure_single_static_mod(internal_sites)

    if len(nterm_sites) > 0:
        annotation.extend_nterm_mods(nterm_sites.get(0, []))
    if len(cterm_sites) > 0:
        annotation.extend_cterm_mods(cterm_sites.get(len(annotation) - 1, []))
    if len(internal_sites) > 0:
        for site, mods in internal_sites.items():
            annotation.extend_internal_mods_at_index(site, mods)

    return annotation


def apply_static_mods_infront(
    annotation: ProFormaAnnotation,
    internal_static: Mapping[str, Iterable[Any]] | None = None,
    is_regex: bool = False,
) -> ProFormaAnnotation:
    if is_regex:
        raise NotImplementedError("Regex not supported for infront static mods")

    if internal_static is None:
        return annotation

    for mod_aa, mod_values in internal_static.items():
        for mod in mod_values:
            annotation.add_static_mod_by_residue(residue=mod_aa, mod=mod)

    return annotation


def modify(
    annotation: ProFormaAnnotation,
    nterm_static: Mapping[str | None, Iterable[Any]] | None = None,
    cterm_static: Mapping[str | None, Iterable[Any]] | None = None,
    internal_static: Mapping[str | None, Iterable[Any]] | None = None,
    labile_static: Mapping[str | None, Iterable[Any]] | None = None,
    nterm_variable: Mapping[str | None, Iterable[Any]] | None = None,
    cterm_variable: Mapping[str | None, Iterable[Any]] | None = None,
    internal_variable: Mapping[str | None, Iterable[Any]] | None = None,
    labile_variable: Mapping[str | None, Iterable[Any]] | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    use_static_notation: bool = False,
    unique_peptidoforms: bool = False,
    inplace: bool = False,
) -> Generator[ProFormaAnnotation]:
    """
    Generate all possible combinations of modifications for a peptide.

    Args:
        annotation: The base ProFormaAnnotation to modify
        nterm_static: Static modifications for N-terminus
        cterm_static: Static modifications for C-terminus
        internal_static: Static modifications for internal residues
        labile_static: Static labile modifications
        nterm_variable: Variable modifications for N-terminus
        cterm_variable: Variable modifications for C-terminus
        internal_variable: Variable modifications for internal residues
        labile_variable: Variable labile modifications
        max_variable_mods: Maximum number of variable modifications to apply
        use_regex: Whether to use regex for modification site matching
        unique_peptidoforms: If True, only yield unique modification patterns (by mod composition, not position)

    Yields:
        ProFormaAnnotation: Each possible modification combination
    """

    # Get all modification sites
    nterm_variable_sites = (
        get_sites(annotation.stripped_sequence, nterm_variable, use_regex)
        if nterm_variable
        else {}
    )
    cterm_variable_sites = (
        get_sites(annotation.stripped_sequence, cterm_variable, use_regex)
        if cterm_variable
        else {}
    )
    internal_variable_sites = (
        get_sites(annotation.stripped_sequence, internal_variable, use_regex)
        if internal_variable
        else {}
    )

    # Get terminal modifications
    nterm_variable_mods: list[Any] = nterm_variable_sites.get(0, [])
    cterm_variable_mods: list[Any] = cterm_variable_sites.get(len(annotation) - 1, [])

    # Create list of all possible variable modification sites with their mods
    variable_site_mod_pairs: list[tuple[str, int, Any]] = []

    # Add N-term variable mods
    for mod in nterm_variable_mods:
        variable_site_mod_pairs.append(("nterm", 0, mod))

    # Add C-term variable mods
    for mod in cterm_variable_mods:
        variable_site_mod_pairs.append(("cterm", len(annotation.sequence) - 1, mod))

    # Add internal variable mods (including terminal positions since they're treated separately)
    for site, mods in internal_variable_sites.items():
        for mod in mods:
            variable_site_mod_pairs.append(("internal", site, mod))

    # Add labile variable mods
    if labile_variable:
        labile_variable_sites = get_sites(
            annotation.sequence, labile_variable, use_regex
        )
        for site, mods in labile_variable_sites.items():
            for mod in mods:
                variable_site_mod_pairs.append(("labile", site, mod))

    if use_static_notation:
        static_modified_annotation = apply_static_mods_infront(
            annotation,
            internal_static=internal_static,
            is_regex=use_regex,
        )
    else:
        static_modified_annotation = apply_mods(
            annotation,
            nterm=nterm_static,
            cterm=cterm_static,
            internal=internal_static,
            inplace=inplace,
            is_regex=use_regex,
        )

    if labile_static:
        labile_static_sites = get_sites(annotation.sequence, labile_static, use_regex)
        for _, mods in labile_static_sites.items():
            for mod in mods:
                static_modified_annotation.append_labile_mod(mod)

    # Track seen modification combinations if first_proteoform_only
    seen_mod_combinations: set[tuple[Any, ...]] | None = None
    if unique_peptidoforms:
        seen_mod_combinations = set()

    # Generate all combinations of variable modifications up to max_variable_mods
    for num_var_mods in range(max_variable_mods + 1):
        for var_mod_combination in itertools.combinations(
            variable_site_mod_pairs, num_var_mods
        ):
            # Check for site conflicts
            site_conflicts = {}
            has_conflict = False

            for mod_type, site, _ in var_mod_combination:
                key = (mod_type, site)
                if key in site_conflicts:
                    has_conflict = True
                    break
                site_conflicts[key] = True

            if has_conflict:
                continue

            # OPTIMIZATION: Check uniqueness early if first_proteoform_only
            if unique_peptidoforms:
                # Create a hashable signature of the modification combination
                # Sort by (mod_type, mod_value) to make position-independent
                mod_signature = tuple(
                    sorted(
                        (mod_type, _get_mod_value(mod))
                        for mod_type, _, mod in var_mod_combination
                    )
                )

                if seen_mod_combinations is None:
                    raise RuntimeError("seen_mod_combinations should be initialized")

                if mod_signature in seen_mod_combinations:
                    continue
                seen_mod_combinations.add(mod_signature)

            # Start with base annotation and apply static mods
            modified_annotation = static_modified_annotation.copy()

            # Apply variable modifications directly
            for mod_type, site, mod in var_mod_combination:
                if mod_type == "nterm":
                    modified_annotation.append_nterm_mod(mod)

                elif mod_type == "cterm":
                    modified_annotation.append_cterm_mod(mod)

                elif mod_type == "internal":
                    modified_annotation.append_internal_mod_at_index(site, mod)

                elif mod_type == "labile":
                    modified_annotation.append_labile_mod(mod)

            yield modified_annotation


def _get_mod_value(mod: Any) -> str:
    """Helper function to get a hashable mod value for uniqueness checking."""
    if isinstance(mod, str):
        return mod
    return str(mod)

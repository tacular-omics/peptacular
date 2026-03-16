from collections.abc import Sequence
from typing import Any, overload

from ..annotation import (
    ProFormaAnnotation,
)
from ..constants import parallelMethod, parallelMethodLiteral
from .parallel import parallel_apply_internal
from .util import get_annotation_input


def _parse_chimeric_single(s: str, validate: bool = False) -> list[ProFormaAnnotation]:
    return list(ProFormaAnnotation.parse_chimeric(s, validate=validate))


@overload
def parse_chimeric(
    s: str,
    validate: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[ProFormaAnnotation]: ...


@overload
def parse_chimeric(
    s: Sequence[str],
    validate: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[list[ProFormaAnnotation]]: ...


def parse_chimeric(
    s: str | Sequence[str],
    validate: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[ProFormaAnnotation] | list[list[ProFormaAnnotation]]:
    """Parse a chimeric ProForma string or list of strings into lists of ProFormaAnnotation objects."""
    if isinstance(s, Sequence) and not isinstance(s, str):
        return parallel_apply_internal(
            _parse_chimeric_single,
            s,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            validate=validate,
            reuse_pool=reuse_pool,
        )
    else:
        return _parse_chimeric_single(s, validate=validate)


def _serialize_chimeric_single(
    sequence: Sequence[ProFormaAnnotation | str],
) -> str:
    annots = [get_annotation_input(seq, copy=True) for seq in sequence]

    # first annot will have comound anme and global mods (iso and static)

    # ensure all annots share the same compound name and global mods
    compound_names = {annot.compound_name for annot in annots}
    if len(compound_names) > 1:
        raise ValueError("All annotations in a chimeric sequence must share the same compound name.")

    static_mods = {annot.static_mods for annot in annots}

    if len(static_mods) > 1:
        raise ValueError("All annotations in a chimeric sequence must share the same static modifications.")

    isotope_mods = {annot.isotope_mods for annot in annots}

    if len(isotope_mods) > 1:
        raise ValueError("All annotations in a chimeric sequence must share the same isotopic modifications.")

    for i, annot in enumerate(annots):
        if i == 0:
            continue
        annot.compound_name = None
        annot.static_mods = None
        annot.isotope_mods = None

    return "+".join(annot.serialize() for annot in annots)


@overload
def serialize_chimeric(
    sequence: Sequence[ProFormaAnnotation | str],
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> str: ...


@overload
def serialize_chimeric(
    sequence: Sequence[Sequence[ProFormaAnnotation | str]],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[str]: ...


def serialize_chimeric(
    sequence: Sequence[ProFormaAnnotation | str] | Sequence[Sequence[ProFormaAnnotation | str]],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> str | list[str]:
    """Serialize a chimeric peptide sequence or list of sequences to ProForma string format."""
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and all(isinstance(seq, Sequence) and not isinstance(seq, str) for seq in sequence):
        return parallel_apply_internal(
            _serialize_chimeric_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _serialize_chimeric_single(sequence)  # type: ignore[arg-type]


def _parse_single(s: str, validate: bool = False) -> ProFormaAnnotation:
    return ProFormaAnnotation.parse(s, validate=validate)


@overload
def parse(
    s: str,
    validate: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> ProFormaAnnotation: ...


@overload
def parse(
    s: Sequence[str],
    validate: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[ProFormaAnnotation]: ...


def parse(
    s: str | Sequence[str],
    validate: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> ProFormaAnnotation | list[ProFormaAnnotation]:
    """Parse a ProForma string or list of strings into ProFormaAnnotation object(s)."""
    if isinstance(s, Sequence) and not isinstance(s, str):
        return parallel_apply_internal(
            _parse_single,
            s,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            validate=validate,
            reuse_pool=reuse_pool,
        )
    else:
        return _parse_single(s, validate=validate)


def _serialize_single(
    sequence: str | ProFormaAnnotation,
) -> str:
    return get_annotation_input(sequence, copy=False).serialize()


@overload
def serialize(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> str: ...


@overload
def serialize(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[str]: ...


def serialize(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> str | list[str]:
    """Serialize a peptide sequence or list of sequences to ProForma string format."""
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _serialize_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _serialize_single(sequence)


def _sequence_length_single(sequence: str | ProFormaAnnotation) -> int:
    return len(get_annotation_input(sequence, copy=False))


@overload
def sequence_length(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> int: ...


@overload
def sequence_length(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[int]: ...


def sequence_length(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> int | list[int]:
    """Compute the length of the peptide sequence based on the unmodified sequence."""
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _sequence_length_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _sequence_length_single(sequence)


def _is_ambiguous_single(sequence: str | ProFormaAnnotation) -> bool:
    return get_annotation_input(sequence, copy=False).has_sequence_ambiguity


@overload
def is_ambiguous(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> bool: ...


@overload
def is_ambiguous(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[bool]: ...


def is_ambiguous(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> bool | list[bool]:
    """Check if the sequence contains ambiguous amino acids."""
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _is_ambiguous_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _is_ambiguous_single(sequence)


def _is_modified_single(sequence: str | ProFormaAnnotation) -> bool:
    return get_annotation_input(sequence, copy=False).has_mods()


@overload
def is_modified(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> bool: ...


@overload
def is_modified(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[bool]: ...


def is_modified(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> bool | list[bool]:
    """Check if the sequence contains any modifications."""
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _is_modified_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _is_modified_single(sequence)


def _count_residues_single(sequence: str | ProFormaAnnotation, include_mods: bool = True) -> dict[str, int]:
    return get_annotation_input(sequence, copy=False).condense_static_mods(inplace=True).count_residues(include_mods=include_mods)


@overload
def count_residues(
    sequence: str | ProFormaAnnotation,
    include_mods: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> dict[str, int]: ...


@overload
def count_residues(
    sequence: Sequence[str | ProFormaAnnotation],
    include_mods: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[dict[str, int]]: ...


def count_residues(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    include_mods: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> dict[str, int] | list[dict[str, int]]:
    """
    Counts the occurrences of each amino acid in the input sequence.

    include_mods: If True, modified residues are counted as distinct entities. If False,
                  only unmodified residues are counted.

    .. code-block:: python

        # Single sequence
        >>> count_residues('PEPTIDE')
        {'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1}

        # Single sequence
        >>> count_residues('PEP[Oxidation]TIDE-[+30]')
        {'P': 1, 'E': 1, 'P[Oxidation]': 1, 'T': 1, 'I': 1, 'D': 1, 'E-[+30]': 1}

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _count_residues_single,
            sequence,
            include_mods=include_mods,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _count_residues_single(sequence, include_mods=include_mods)


def _percent_residues_single(
    sequence: str | ProFormaAnnotation,
    include_mods: bool = True,
) -> dict[str, float]:
    return get_annotation_input(sequence, copy=False).condense_static_mods(inplace=True).percent_residues(include_mods=include_mods)


@overload
def percent_residues(
    sequence: str | ProFormaAnnotation,
    include_mods: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> dict[str, float]: ...


@overload
def percent_residues(
    sequence: Sequence[str | ProFormaAnnotation],
    include_mods: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[dict[str, float]]: ...


def percent_residues(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    include_mods: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> dict[str, float] | list[dict[str, float]]:
    """
    Calculates the percentage of each amino acid in the input sequence.

    include_mods: If True, modified residues are counted as distinct entities. If False,
                  only unmodified residues are counted.

    .. code-block:: python

        # Single sequence
        >>> d = percent_residues('PEPTIDE')
        >>> dict(map(lambda item: (item[0], round(item[1], 2)), d.items()))
        {'P': 28.57, 'E': 28.57, 'T': 14.29, 'I': 14.29, 'D': 14.29}

        # Single sequence with modification
        >>> d = percent_residues('PEP[Oxidation]TIDE-[+30]')
        >>> dict(map(lambda item: (item[0], round(item[1], 2)), d.items()))
        {'P': 14.29, 'E': 14.29, 'P[Oxidation]': 14.29, 'T': 14.29, 'I': 14.29, 'D': 14.29, 'E-[+30]': 14.29}

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _percent_residues_single,
            sequence,
            include_mods=include_mods,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _percent_residues_single(sequence, include_mods=include_mods)


def annotate_ambiguity(
    sequence: str | ProFormaAnnotation,
    forward_coverage: list[int],
    reverse_coverage: list[int],
    mass_shift: Any | None = None,
    add_mods_to_intervals: bool = False,
    sort_mods: bool = True,
    condense_to_xnotation: bool = False,
) -> str:
    """
    This function identifies regions in the sequence where there is insufficient fragment ion
    coverage and marks them as ambiguous using ProForma notation with parentheses.
    If a mass shift is provided, it will be added to the appropriate location.

    forward_coverage: Binary list indicating which positions have forward ion coverage (1) or not (0).
    reverse_coverage: Binary list indicating which positions have reverse ion coverage (1) or not (0).
    mass_shift: An optional mass shift to be added to the sequence at the appropriate position.
    add_mods_to_intervals: Whether to add modifications to interval annotations.
    sort_mods: Whether to sort modifications.
    condense_to_xnotation: Whether to condense ambiguity to X notation.

    .. code-block:: python

        # Add ambiguity intervals based on fragment ion coverage
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,1,0,0,0], [0,0,0,0,0,1,0])
        '(?PE)PTI(?DE)'

        # With a phosphorylation mass shift (note the '+' sign)
        >>> annotate_ambiguity('PEPTIDE', [1,1,1,0,0,0,0], [0,0,0,0,1,1,1], 79.966)
        'PEPT[+79.966]IDE'

        # Handling existing modifications
        >>> annotate_ambiguity('P[+10]EPTIDE', [1,1,1,0,0,0,0], [0,0,0,0,0,1,1])
        'P[+10]EP(?TI)DE'

        # When mass shift can't be localized to a specific residue
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,0,0,0,0], [0,0,0,0,0,1,0], 120)
        '(?PE)P(?TI)[+120](?DE)'

        # When mass shift is completely unlocalized, it becomes a labile modification
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,1,1,0,0], [0,0,1,1,1,1,0], 120)
        '{+120}(?PE)PTI(?DE)'

        # Complex example with multiple intervals
        >>> for_ions = list(map(int, '00011101001000000000000000000000000000'))
        >>> rev_ions = list(map(int, '00000000000110000000101111111111010100'))
        >>> annotate_ambiguity('SSGSIASSYVQWYQQRPGSAPTTVIYEDDERPSGVPDR', for_ions, rev_ions, 120)
        '(?SSGS)IA(?SS)(?YVQ)W[+120](?YQQRPGSA)(?PT)TVIYEDDER(?PS)(?GV)(?PDR)'
    """
    annot = get_annotation_input(sequence=sequence, copy=True).annotate_ambiguity(
        forward_coverage=forward_coverage,
        reverse_coverage=reverse_coverage,
        mass_shift=mass_shift,
        add_mods_to_intervals=add_mods_to_intervals,
        sort_mods=sort_mods,
        inplace=True,
    )

    if condense_to_xnotation:
        annot.condense_ambiguity_to_xnotation(inplace=True)

    return annot.serialize()


def _validate_single(
    sequence: str | ProFormaAnnotation,
) -> bool:
    try:
        get_annotation_input(sequence, copy=False).validate_annotation()
    except ValueError:
        return False
    return True


@overload
def validate(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> bool: ...


@overload
def validate(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[bool]: ...


def validate(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> bool | list[bool]:
    """
    Checks if the input sequence is a valid ProForma sequence.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _validate_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _validate_single(sequence)


def _generate_random_single(
    _item: Any = None,  # Dummy parameter for parallel processing
    min_length: int = 6,
    max_length: int = 20,
    mod_probability: float = 0.05,
    include_internal_mods: bool = True,
    include_nterm_mods: bool = True,
    include_cterm_mods: bool = True,
    include_labile_mods: bool = True,
    include_unknown_mods: bool = True,
    include_isotopic_mods: bool = True,
    include_static_mods: bool = True,
    include_intervals: bool = True,
    include_charge: bool = True,
    require_composition: bool = True,
) -> ProFormaAnnotation:
    return ProFormaAnnotation.random(
        min_length=min_length,
        max_length=max_length,
        mod_probability=mod_probability,
        include_internal_mods=include_internal_mods,
        include_nterm_mods=include_nterm_mods,
        include_cterm_mods=include_cterm_mods,
        include_labile_mods=include_labile_mods,
        include_unknown_mods=include_unknown_mods,
        include_isotopic_mods=include_isotopic_mods,
        include_static_mods=include_static_mods,
        include_intervals=include_intervals,
        include_charge=include_charge,
        require_composition=require_composition,
    )


@overload
def generate_random(
    count: None = None,
    min_length: int = 6,
    max_length: int = 20,
    mod_probability: float = 0.05,
    include_internal_mods: bool = True,
    include_nterm_mods: bool = True,
    include_cterm_mods: bool = True,
    include_labile_mods: bool = True,
    include_unknown_mods: bool = True,
    include_isotopic_mods: bool = True,
    include_static_mods: bool = True,
    include_intervals: bool = True,
    include_charge: bool = True,
    require_composition: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> ProFormaAnnotation: ...


@overload
def generate_random(
    count: int,
    min_length: int = 6,
    max_length: int = 20,
    mod_probability: float = 0.05,
    include_internal_mods: bool = True,
    include_nterm_mods: bool = True,
    include_cterm_mods: bool = True,
    include_labile_mods: bool = True,
    include_unknown_mods: bool = True,
    include_isotopic_mods: bool = True,
    include_static_mods: bool = True,
    include_intervals: bool = True,
    include_charge: bool = True,
    require_composition: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[ProFormaAnnotation]: ...


def generate_random(
    count: int | None = None,
    min_length: int = 6,
    max_length: int = 20,
    mod_probability: float = 0.05,
    include_internal_mods: bool = True,
    include_nterm_mods: bool = True,
    include_cterm_mods: bool = True,
    include_labile_mods: bool = True,
    include_unknown_mods: bool = True,
    include_isotopic_mods: bool = True,
    include_static_mods: bool = True,
    include_intervals: bool = True,
    include_charge: bool = True,
    require_composition: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> ProFormaAnnotation | list[ProFormaAnnotation]:
    """Generate random ProForma annotation(s) with configurable features.

    Args:
        count: Number of random sequences to generate. If None, generates a single sequence.
        min_length: Minimum sequence length
        max_length: Maximum sequence length
        mod_probability: Probability of adding modifications (0.0 to 1.0)
        include_internal_mods: Whether to generate internal modifications
        include_nterm_mods: Whether to generate N-terminal modifications
        include_cterm_mods: Whether to generate C-terminal modifications
        include_labile_mods: Whether to generate labile modifications
        include_unknown_mods: Whether to generate unknown position modifications
        include_isotopic_mods: Whether to generate isotopic modifications
        include_static_mods: Whether to generate static modifications
        include_intervals: Whether to generate intervals
        include_charge: Whether to generate charge state or adduct
        require_composition: If True, only modifications with composition are allowed (no mass-only)
        n_workers: Number of parallel workers (only used when count > 1)
        chunksize: Size of chunks for parallel processing
        method: Parallel processing method ('process', 'thread', or 'sequential')

    Returns:
        A single ProFormaAnnotation if count is None, otherwise a list of ProFormaAnnotations

    .. code-block:: python

        # Generate a single random sequence
        >>> seq = generate_random()
        >>> isinstance(seq, ProFormaAnnotation)
        True

        # Generate multiple random sequences
        >>> seqs = generate_random(count=10)
        >>> len(seqs)
        10

        # Generate without modifications
        >>> seq = generate_random(mod_probability=0.0)

        # Generate with only internal modifications
        >>> seq = generate_random(
        ...     include_nterm_mods=False,
        ...     include_cterm_mods=False,
        ...     include_labile_mods=False
        ... )
    """
    if count is None:
        return _generate_random_single(
            min_length=min_length,
            max_length=max_length,
            mod_probability=mod_probability,
            include_internal_mods=include_internal_mods,
            include_nterm_mods=include_nterm_mods,
            include_cterm_mods=include_cterm_mods,
            include_labile_mods=include_labile_mods,
            include_unknown_mods=include_unknown_mods,
            include_isotopic_mods=include_isotopic_mods,
            include_static_mods=include_static_mods,
            include_intervals=include_intervals,
            include_charge=include_charge,
            require_composition=require_composition,
        )
    else:
        # Generate count number of items
        items = [None] * count
        return parallel_apply_internal(
            _generate_random_single,
            items,
            min_length=min_length,
            max_length=max_length,
            mod_probability=mod_probability,
            include_internal_mods=include_internal_mods,
            include_nterm_mods=include_nterm_mods,
            include_cterm_mods=include_cterm_mods,
            include_labile_mods=include_labile_mods,
            include_unknown_mods=include_unknown_mods,
            include_isotopic_mods=include_isotopic_mods,
            include_static_mods=include_static_mods,
            include_intervals=include_intervals,
            include_charge=include_charge,
            require_composition=require_composition,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )

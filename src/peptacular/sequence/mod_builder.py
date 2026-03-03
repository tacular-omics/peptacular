from collections.abc import Iterable, Mapping, Sequence
from typing import Any, cast, overload

from ..annotation import ProFormaAnnotation
from ..constants import ModType, ModTypeLiteral, parallelMethod, parallelMethodLiteral
from .parallel import parallel_apply_internal
from .util import get_annotation_input

MOD_BUILDER_INPUT_TYPE = Mapping[str, Iterable[Any]]


def _modify_single(
    sequence: ProFormaAnnotation | str,
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
) -> list[str]:
    """Build modifications for a single sequence"""
    annotation = get_annotation_input(sequence, copy=True)
    return [
        annot.serialize()
        for annot in annotation.modify(
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
            inplace=False,
        )
    ]


@overload
def modify(
    sequence: ProFormaAnnotation | str,
    *,  # force keyword arguments
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def modify(
    sequence: Sequence[ProFormaAnnotation | str],
    *,  # force keyword arguments
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def modify(
    sequence: ProFormaAnnotation | str | Sequence[ProFormaAnnotation | str],
    *,  # force keyword arguments
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Build modified sequences by applying static and variable modifications to a sequence or list of sequences.

    Modifications are specified as dictionaries where keys represent the residue or terminus type and values are iterables of modifications.

    .. code-block:: python

        # Single sequence
        >>> results = modify('PEPTIDE', internal_variable={'P': [79.966]}, max_variable_mods=1)
        >>> len(results) > 0
        True

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> results = modify(sequences, internal_variable={'P': [79.966]}, max_variable_mods=1)
        >>> len(results)
        3


    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _modify_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
        )
    else:
        return _modify_single(
            sequence=sequence,
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
        )


def get_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | ModTypeLiteral | None = None,
) -> dict[ModType | ModTypeLiteral, Any]:
    """
    Parses a sequence with modifications and returns a dictionary where keys represent the position/type of the modifications.
    """

    return get_annotation_input(sequence, copy=True).get_mods(mods)


def set_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any] | None,
) -> str:
    return get_annotation_input(sequence, copy=True).set_mods(mods, inplace=True).serialize()


def append_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
) -> str:
    return get_annotation_input(sequence, copy=True).append_mods(mods, inplace=True).serialize()


def extend_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
) -> str:
    return get_annotation_input(sequence, copy=True).extend_mods(mods, inplace=True).serialize()


def condense_static_mods(
    sequence: str | ProFormaAnnotation,
) -> str:
    """
    Condenses static modifications into internal modifications.

    .. code-block:: python

        # Condenses static modifications to specified internal modifications
        >>> condense_static_mods('<13C><[100]@P>PEPTIDE')
        '<13C>P[100]EP[100]TIDE'

        # If residue is already modified, the static modification will be appended
        >>> condense_static_mods('<13C><[100]@P>P[10]EPTIDE')
        '<13C>P[10][100]EP[100]TIDE'

        # Example for unmodified sequences
        >>> condense_static_mods('PEPTIDE')
        'PEPTIDE'

        # Example for N-Term static modifications
        >>> condense_static_mods('<[Oxidation]@N-term>PEPTIDE')
        '[Oxidation]-PEPTIDE'

        # Example for C-Term static modifications
        >>> condense_static_mods('<[Oxidation]@C-term>PEPTIDE')
        'PEPTIDE-[Oxidation]'

    """

    return get_annotation_input(sequence=sequence, copy=True).condense_static_mods(inplace=False).serialize()


def pop_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
) -> tuple[str, dict[ModType, Any]]:
    """
    Removes all modifications from the given sequence, returning the unmodified sequence and a dictionary of the
    removed modifications.

    .. code-block:: python

        # Simply combines the functionality of strip_modifications and get_modifications
        >>> seq, mod_dict = pop_mods('PEP[phospho]TIDE')
        >>> seq
        'PEPTIDE'

    """
    annotation: ProFormaAnnotation = get_annotation_input(sequence=sequence, copy=True)
    mod_dict: dict[ModType, Any] = annotation.pop_mods(mod_types=mods)
    return (
        annotation.serialize(),
        mod_dict,
    )


def remove_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
) -> str:
    annotation = get_annotation_input(sequence=sequence, copy=True)

    return annotation.clear_mods(mods=mods, inplace=True).serialize()


def _strip_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
) -> str:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.clear_mods(mods=mods, inplace=True).serialize()


@overload
def strip_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
) -> str: ...


@overload
def strip_mods(
    sequence: Sequence[str | ProFormaAnnotation],
    mods: ModType | Iterable[ModType] | None = None,
) -> list[str]: ...


def strip_mods(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    mods: ModType | Iterable[ModType] | None = None,
) -> str | list[str]:
    """
    Strips all modifications from the given sequence or list of sequences, returning the unmodified sequence(s).

    .. code-block:: python

        # Removes internal modifications:
        >>> strip_mods('PEP[phospho]TIDE')
        'PEPTIDE'
    """

    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _strip_mods,
            sequence,
            mods=mods,
        )
    else:
        return _strip_mods(
            sequence=sequence,
            mods=mods,
        )


def filter_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
) -> str:
    """
    Keeps only the specified modifications in the sequence, removing all others.

    .. code-block:: python

        # Keeps only internal modifications:
        >>> filter_mods('PEP[phospho]TIDE', mods='internal')
        'PEP[phospho]TIDE'

    """
    return get_annotation_input(sequence=sequence, copy=True).filter_mods(mods=mods, inplace=True).serialize()


def _to_ms2_pip_single(
    sequence: ProFormaAnnotation | str,
) -> tuple[str, str]:
    annot = get_annotation_input(sequence=sequence, copy=True)
    return annot.to_ms2_pip(inplace=True)


@overload
def to_ms2_pip(
    sequence: ProFormaAnnotation | str,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> tuple[str, str]: ...


@overload
def to_ms2_pip(
    sequence: Sequence[ProFormaAnnotation | str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[tuple[str, str]]: ...


def to_ms2_pip(
    sequence: ProFormaAnnotation | str | Sequence[ProFormaAnnotation | str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> tuple[str, str] | list[tuple[str, str]]:
    """
    Convert a peptide sequence to MS2PIP format by condensing modifications.

    .. code-block:: python

        # Single sequence
        >>> to_ms2_pip('PEP[Phospho]TIDE')
        ('PEPTIDE', '3|Phospho')

        # Batch processing
        >>> sequences = ['PEP[Phospho]TIDE', 'PROT[Oxidation]EIN']
        >>> to_ms2_pip(sequences)
        [('PEPTIDE', '3|Phospho'), ('PROTEIN', '4|Oxidation')]

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _to_ms2_pip_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _to_ms2_pip_single(sequence=sequence)


def _from_ms2_pip_single(
    item: tuple[str, str],
    static_mods: Mapping[str, float] | None = None,
) -> str:
    sequence, modifications = item
    return ProFormaAnnotation.from_ms2_pip(
        sequence=sequence,
        mod_str=modifications,
        static_mods=static_mods,
    ).serialize()


@overload
def from_ms2_pip(
    sequence: tuple[str, str],
    static_mods: Mapping[str, float] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> str: ...


@overload
def from_ms2_pip(
    sequence: Sequence[tuple[str, str]],
    static_mods: Mapping[str, float] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[str]: ...


def from_ms2_pip(
    sequence: tuple[str, str] | Sequence[tuple[str, str]],
    static_mods: Mapping[str, float] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Convert MS2PIP format to ProForma string(s).

    .. code-block:: python

        # Single sequence
        >>> from_ms2_pip(('PEPTIDE', '3|Phospho'))
        'PEP[Phospho]TIDE'

        # Batch processing
        >>> items = [('PEPTIDE', '3|Phospho'), ('PROTEIN', '4|Oxidation')]
        >>> from_ms2_pip(items)
        ['PEP[Phospho]TIDE', 'PROT[Oxidation]EIN']

        # Empty modifications
        >>> from_ms2_pip(('PEPTIDE', ''))
        'PEPTIDE'
    """
    # Check if batch processing (list of tuples)
    if isinstance(sequence, Sequence) and not isinstance(sequence, tuple):
        # Validate that all items are tuples
        if not all(isinstance(item, tuple) and len(item) == 2 for item in sequence):
            raise ValueError("All items in sequence must be tuples of (sequence, modifications)")

        return parallel_apply_internal(
            _from_ms2_pip_single,
            sequence,
            static_mods=static_mods,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        # Single tuple processing
        if not isinstance(sequence, tuple) or len(sequence) != 2:
            raise ValueError("sequence must be a tuple of (sequence, modifications) for single processing")

        item = cast(tuple[str, str], sequence)
        return _from_ms2_pip_single(
            item=item,
            static_mods=static_mods,
        )


def _condense_to_peptidoform(
    sequence: str | ProFormaAnnotation,
) -> str:
    return get_annotation_input(sequence=sequence, copy=True).condense_to_peptidoform(inplace=False).serialize()


@overload
def condense_to_peptidoform(
    sequence: str | ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str: ...


@overload
def condense_to_peptidoform(
    sequence: Sequence[str | ProFormaAnnotation],
    include_plus: bool = False,
    precision: int | None = None,
) -> list[str]: ...


def condense_to_peptidoform(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
) -> str | list[str]:
    """
    Condenses all modifications into a peptidoform representation for a sequence or list of sequences.
    """

    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
        return parallel_apply_internal(
            _condense_to_peptidoform,
            sequence,
        )
    else:
        return _condense_to_peptidoform(
            sequence=sequence,
        )

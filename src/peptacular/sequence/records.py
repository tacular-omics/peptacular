"""Digest and fragment results as lists of plain dicts, one per row.

The records hold only ``str``, ``int``, ``float``, ``bool`` and ``None``, one type per key, so
``pandas.DataFrame(records)`` or ``polars.DataFrame(records)`` turns them into a table with
clean column types. peptacular does not depend on pandas or polars.
"""

import re
from collections.abc import Iterable
from typing import Any, cast

from ..annotation import ProFormaAnnotation
from ..annotation.frag import Fragment
from ..diagnostics import PeptacularError
from ..proforma_components import ChargedFormula
from .digestion import _digest_input, _span_output
from .util import HasSequence

__all__ = [
    "DIGEST_RECORD_KEYS",
    "FRAGMENT_RECORD_KEYS",
    "digest_records",
    "fragment_records",
]

DIGEST_RECORD_KEYS: tuple[str, ...] = (
    "peptide",
    "stripped_sequence",
    "start",
    "end",
    "missed_cleavages",
    "semi",
    "accession",
)
"""Keys of every :func:`digest_records` row, in order."""

FRAGMENT_RECORD_KEYS: tuple[str, ...] = (
    "ion_type",
    "position",
    "end_position",
    "charge_state",
    "mz",
    "mass",
    "neutral_mass",
    "monoisotopic",
    "deltas",
    "isotopes",
    "sequence",
    "parent_sequence",
    "mzpaf",
)
"""Keys of every :func:`fragment_records` row, in order."""


def _accession(sequence: object) -> str | None:
    """The entry's ``accession`` (fastatacular), else its ``db_unique_id`` (pefftacular), else None."""
    accession = getattr(sequence, "accession", None)
    if accession is None:
        accession = getattr(sequence, "db_unique_id", None)
    if accession is not None and not isinstance(accession, str):
        accession = str(accession)
    return accession


def _digest_one(
    sequence: str | ProFormaAnnotation | HasSequence,
    enzyme: str | re.Pattern[str],
    missed_cleavages: int,
    semi: bool,
    min_len: int | None,
    max_len: int | None,
) -> list[dict[str, Any]]:
    accession = _accession(sequence)
    annot, plain = _digest_input(sequence)
    spans = list(annot.digest_spans(enzyme, missed_cleavages=missed_cleavages, semi=semi, min_len=min_len, max_len=max_len))
    boundaries: set[int] | None = None
    if semi:
        boundaries = {0, len(annot), *annot.cleavage_sites(enzyme)}
    records: list[dict[str, Any]] = []
    for peptide, span in _span_output(annot, plain, spans):
        stripped = plain[span.start : span.end] if plain is not None else annot.stripped_sequence[span.start : span.end]
        records.append(
            {
                "peptide": peptide,
                "stripped_sequence": stripped,
                "start": span.start,
                "end": span.end,
                "missed_cleavages": span.missed_cleavages,
                "semi": boundaries is not None and not (span.start in boundaries and span.end in boundaries),
                "accession": accession,
            }
        )
    return records


def _is_single_sequence(value: object) -> bool:
    return isinstance(value, (str, ProFormaAnnotation)) or hasattr(value, "sequence")


def digest_records(
    sequence: str | ProFormaAnnotation | HasSequence | Iterable[str | ProFormaAnnotation | HasSequence],
    enzyme: str | re.Pattern[str],
    *,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[dict[str, Any]]:
    """Digest ``sequence`` and return one plain dict per peptide.

    ``enzyme``, ``missed_cleavages``, ``semi``, ``min_len`` and ``max_len`` work as in
    :func:`peptacular.digest`; it runs in this process (no ``n_workers``). ``sequence`` is one
    protein or any iterable of them, such as a list or the generator fastatacular's reader
    returns; many proteins give one flat list of rows, in input order. Every row has these
    keys (``DIGEST_RECORD_KEYS``):

    - ``peptide`` (str): the peptide as ProForma, with any modifications of the protein.
    - ``stripped_sequence`` (str): the residues only.
    - ``start``, ``end`` (int): 0-based, half-open position in the protein (``protein[start:end]``).
    - ``missed_cleavages`` (int): as in the digest's :class:`~peptacular.Span`.
    - ``semi`` (bool): True when one end of the peptide is not an enzyme cleavage site or a
      protein terminus. Always False unless ``semi=True``.
    - ``accession`` (str | None): the input's ``accession`` attribute (a fastatacular entry), or
      else its ``db_unique_id`` (a pefftacular entry), otherwise None.

    :raises UnknownEnzymeError: If ``enzyme`` is a string that names no known protease.
    :raises TypeError: If ``sequence`` (or an item of it) is not a string, annotation or object
        with a ``.sequence`` string, as :func:`peptacular.digest` does.
    :return: One dict per peptide, in the order :func:`peptacular.digest` returns them.

    .. code-block:: python

        >>> rows = digest_records("MKVLATSAGERTIDEK", "trypsin", missed_cleavages=1)
        >>> rows[1]
        {'peptide': 'MKVLATSAGER', 'stripped_sequence': 'MKVLATSAGER', 'start': 0, 'end': 11, 'missed_cleavages': 1, 'semi': False, 'accession': None}
        >>> [r["peptide"] for r in digest_records("TIDEKTIDE", "trypsin", semi=True) if r["semi"]][:3]
        ['T', 'TI', 'TID']
    """
    if _is_single_sequence(sequence) or not isinstance(sequence, Iterable):
        single = cast("str | ProFormaAnnotation | HasSequence", sequence)
        return _digest_one(single, enzyme, missed_cleavages, semi, min_len, max_len)
    return [row for item in sequence for row in _digest_one(item, enzyme, missed_cleavages, semi, min_len, max_len)]


def _format_counts(items: Iterable[tuple[object, int]]) -> str:
    parts: list[str] = []
    for key, count in items:
        if isinstance(key, float):
            label = f"{key:+}"
        elif isinstance(key, ChargedFormula):  # from ``Fragment.deltas``
            label = key.serialize().removeprefix("Formula:")
        else:
            label = str(key)
        parts.append(label if count == 1 else f"{label}^{count}")
    return ",".join(parts)


def _without_charge(sequence: str, cache: dict[str, str]) -> str:
    """``sequence`` as ProForma without its ``/charge`` suffix (the row's ``charge_state`` has it)."""
    if sequence not in cache:
        cache[sequence] = ProFormaAnnotation.parse(sequence).serialize(exclude_charge=True)
    return cache[sequence]


def _iter_fragments(fragments: Iterable[Fragment] | Iterable[Iterable[Fragment]]) -> Iterable[Fragment]:
    for item in fragments:
        if isinstance(item, Fragment):
            yield item
        elif isinstance(item, Iterable) and not isinstance(item, (str, bytes)):
            for inner in item:
                if not isinstance(inner, Fragment):
                    raise PeptacularError(f"fragment_records expects Fragment objects, got {type(inner).__name__} inside a list")
                yield inner
        else:
            raise PeptacularError(f"fragment_records expects Fragment objects or lists of them, got {type(item).__name__}")


def fragment_records(fragments: Iterable[Fragment] | Iterable[Iterable[Fragment]]) -> list[dict[str, Any]]:
    """Turn :class:`~peptacular.Fragment` objects into one plain dict per ion.

    Pass the list :func:`peptacular.fragment` returns for one peptide, or the list of lists it
    returns for several; the rows come out flat, and ``parent_sequence`` says which peptide
    each ion belongs to. Every row has these keys (``FRAGMENT_RECORD_KEYS``), named after the
    :class:`~peptacular.Fragment` constructor arguments:

    - ``ion_type`` (str): the ion letter, e.g. ``"b"``, ``"y"``, ``"by"``.
    - ``position`` (int | None): the ion number (``3`` for b3); for an internal ion, its start.
      None for precursor and neutral ions.
    - ``end_position`` (int | None): the end of an internal ion (its ``Fragment.position`` is
      the ``(position, end_position)`` pair); None for every other ion.
    - ``charge_state`` (int), ``mz`` (float), ``mass`` (float, charged), ``neutral_mass`` (float),
      ``monoisotopic`` (bool).
    - ``deltas`` (str): the fragment's deltas joined by ``","``. Each is a signed formula or
      mass, added ``count`` times, with ``^count`` when the count is not one. Named losses
      such as H2O are stored as negative formulas: water loss is ``"H-2O-1"``, a water gain
      ``"H-2O-1^-1"``, while a plain formula such as ``"C2H2O"`` is a gain; a mass keeps its
      sign (``"-17.0^2"``). ``""`` when there are none.
    - ``isotopes`` (str): the isotope swaps in the same form (``"13C"``, ``"13C^2"``, ``"15N"``);
      ``""`` when there are none.
    - ``sequence`` (str | None): the fragment's own residues as ProForma; ``parent_sequence``
      (str | None): the peptide it came from. Both leave out the charge (see ``charge_state``)
      and are None when the fragment was built without its parent sequence.
    - ``mzpaf`` (str | None): the mzPAF label from :meth:`~peptacular.Fragment.to_mzpaf`, or None
      when mzPAF cannot write the ion: an ion type with no mzPAF form, or a formula delta with
      both positive and negative element counts (``"CH-2"``).

    :raises TypeError: If ``fragments`` is not iterable (``None``, a number), as :func:`peptacular.digest` does.
    :raises PeptacularError: If an item is not a :class:`~peptacular.Fragment` or a list of them.

    .. code-block:: python

        >>> import peptacular as pt
        >>> rows = fragment_records(pt.fragment("PEPTIDE", ion_types=("b",), charges=(1, 2)))
        >>> len(rows), rows[0]["ion_type"], rows[0]["position"], rows[0]["charge_state"], round(rows[0]["mz"], 4)
        (14, 'b', 1, 1, 98.06)
        >>> rows[1]["sequence"], rows[1]["mzpaf"]
        ('PE', 'b2{PE}')
    """
    records: list[dict[str, Any]] = []
    stripped: dict[str, str] = {}
    for fragment in _iter_fragments(fragments):
        raw_isotopes = fragment._isotopes
        if isinstance(raw_isotopes, int):
            isotopes = _format_counts([("13C", raw_isotopes)]) if raw_isotopes else ""
        else:
            isotopes = _format_counts((raw_isotopes or {}).items())
        has_parent = bool(fragment.parent_sequence) and fragment.parent_sequence_length is not None
        position, end_position = fragment.position if isinstance(fragment.position, tuple) else (fragment.position, None)
        sequence = parent = None
        if has_parent:
            assert fragment.parent_sequence is not None
            own = fragment.sequence
            sequence = _without_charge(own, stripped) if own is not None else None
            parent = _without_charge(fragment.parent_sequence, stripped)
        records.append(
            {
                "ion_type": fragment.ion_type.value,
                "position": position,
                "end_position": end_position,
                "charge_state": fragment.charge_state,
                "mz": fragment.mz,
                "mass": fragment.mass,
                "neutral_mass": fragment.neutral_mass,
                "monoisotopic": fragment.monoisotopic,
                "deltas": _format_counts(fragment.deltas.items()),
                "isotopes": isotopes,
                "sequence": sequence,
                "parent_sequence": parent,
                "mzpaf": _mzpaf_or_none(fragment, has_parent),
            }
        )
    return records


def _mzpaf_or_none(fragment: Fragment, include_sequence: bool) -> str | None:
    try:
        return fragment.to_mzpaf(include_sequence=include_sequence)
    except PeptacularError:
        return None

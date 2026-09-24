"""Digest and fragment results as lists of plain dicts, one per row.

The records hold only ``str``, ``int``, ``float``, ``bool`` and ``None`` (plus a ``(start, end)``
tuple for an internal ion's position), so ``pandas.DataFrame(records)`` or
``polars.DataFrame(records)`` turns them into a table. peptacular does not depend on pandas
or polars.
"""

import re
from collections.abc import Iterable, Sequence
from typing import Any

from ..annotation import ProFormaAnnotation
from ..annotation.frag import Fragment
from ..diagnostics import PeptacularError
from .digestion import _digest_input, _span_output
from .util import HasSequence

__all__ = [
    "DIGEST_RECORD_FIELDS",
    "FRAGMENT_RECORD_FIELDS",
    "digest_records",
    "fragment_records",
]

DIGEST_RECORD_FIELDS: tuple[str, ...] = (
    "peptide",
    "stripped_sequence",
    "start",
    "end",
    "missed_cleavages",
    "semi",
    "accession",
)
"""Keys of every :func:`digest_records` row, in order."""

FRAGMENT_RECORD_FIELDS: tuple[str, ...] = (
    "ion_type",
    "position",
    "charge_state",
    "mz",
    "mass",
    "neutral_mass",
    "monoisotopic",
    "losses",
    "isotopes",
    "sequence",
    "parent_sequence",
    "mzpaf",
)
"""Keys of every :func:`fragment_records` row, in order."""


def _digest_one(
    sequence: str | ProFormaAnnotation | HasSequence,
    enzyme: str | re.Pattern[str],
    missed_cleavages: int,
    semi: bool,
    min_len: int | None,
    max_len: int | None,
) -> list[dict[str, Any]]:
    accession = getattr(sequence, "accession", None)
    if accession is not None and not isinstance(accession, str):
        accession = str(accession)
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


def digest_records(
    sequence: str | ProFormaAnnotation | HasSequence | Sequence[str | ProFormaAnnotation | HasSequence],
    enzyme: str | re.Pattern[str],
    *,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[dict[str, Any]]:
    """Digest ``sequence`` and return one plain dict per peptide.

    Takes the same arguments as :func:`peptacular.digest`. A list of proteins gives one flat
    list of rows, in input order. Every row has these keys (``DIGEST_RECORD_FIELDS``):

    - ``peptide`` (str): the peptide as ProForma, with any modifications of the protein.
    - ``stripped_sequence`` (str): the residues only.
    - ``start``, ``end`` (int): 0-based, half-open position in the protein (``protein[start:end]``).
    - ``missed_cleavages`` (int): as in the digest's :class:`~peptacular.Span`.
    - ``semi`` (bool): True when one end of the peptide is not an enzyme cleavage site or a
      protein terminus. Always False unless ``semi=True``.
    - ``accession`` (str | None): the input's ``accession`` attribute (a fastatacular or
      pefftacular entry has one), otherwise None.

    :raises UnknownEnzymeError: If ``enzyme`` is a string that names no known protease.
    :return: One dict per peptide, in the order :func:`peptacular.digest` returns them.

    .. code-block:: python

        >>> rows = digest_records("MKVLATSAGERTIDEK", "trypsin", missed_cleavages=1)
        >>> rows[1]
        {'peptide': 'MKVLATSAGER', 'stripped_sequence': 'MKVLATSAGER', 'start': 0, 'end': 11, 'missed_cleavages': 1, 'semi': False, 'accession': None}
        >>> [r["peptide"] for r in digest_records("TIDEKTIDE", "trypsin", semi=True) if r["semi"]][:3]
        ['T', 'TI', 'TID']
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, (str, ProFormaAnnotation)):
        return [row for item in sequence for row in _digest_one(item, enzyme, missed_cleavages, semi, min_len, max_len)]
    return _digest_one(sequence, enzyme, missed_cleavages, semi, min_len, max_len)


def _format_counts(items: Iterable[tuple[object, int]]) -> str:
    parts: list[str] = []
    for key, count in items:
        label = f"{key:+}" if isinstance(key, float) else str(key)
        parts.append(label if count == 1 else f"{label}^{count}")
    return ",".join(parts)


def _mzpaf_or_none(fragment: Fragment, include_sequence: bool) -> str | None:
    try:
        return fragment.to_mzpaf(include_sequence=include_sequence)
    except PeptacularError:
        return None


def fragment_records(fragments: Iterable[Fragment]) -> list[dict[str, Any]]:
    """Turn :class:`~peptacular.Fragment` objects into one plain dict per ion.

    Pass the list :func:`peptacular.fragment` returns for one peptide (for a batch, one list per
    peptide). Keys match the :class:`~peptacular.Fragment` attributes (``FRAGMENT_RECORD_FIELDS``):

    - ``ion_type`` (str): the ion letter, e.g. ``"b"``, ``"y"``, ``"i"``.
    - ``position`` (int | tuple[int, int] | None): the ion number, a ``(start, end)`` tuple for
      internal ions, None for precursor and neutral ions.
    - ``charge_state`` (int), ``mz`` (float), ``mass`` (float, charged), ``neutral_mass`` (float),
      ``monoisotopic`` (bool).
    - ``losses`` (str): the deltas and neutral losses as ProForma formulas or signed masses,
      joined by ``","``, with ``^n`` for a count above one (``"H-2O-1"``, ``"-17.0^2"``); ``""``
      when there are none.
    - ``isotopes`` (str): the isotope swaps in the same form (``"13C"``, ``"13C^2"``, ``"15N"``);
      ``""`` when there are none.
    - ``sequence`` (str | None): the fragment's own residues as ProForma; ``parent_sequence``
      (str | None): the peptide it came from. Both are None when the fragment was built without
      its parent sequence.
    - ``mzpaf`` (str | None): the mzPAF label from :meth:`~peptacular.Fragment.to_mzpaf`, or None
      when mzPAF cannot express the ion (a bare mass delta such as ``-17.0``).

    :raises TypeError: If an item is not a :class:`~peptacular.Fragment`.

    .. code-block:: python

        >>> import peptacular as pt
        >>> rows = fragment_records(pt.fragment("PEPTIDE", ion_types=("b",), charges=(1, 2)))
        >>> len(rows), rows[0]["ion_type"], rows[0]["position"], rows[0]["charge_state"], round(rows[0]["mz"], 4)
        (14, 'b', 1, 1, 98.06)
        >>> rows[1]["sequence"], rows[1]["mzpaf"]
        ('PE/1', 'b2{PE}')
    """
    records: list[dict[str, Any]] = []
    for fragment in fragments:
        if not isinstance(fragment, Fragment):
            raise TypeError(f"fragment_records expects Fragment objects, got {type(fragment).__name__}")
        raw_isotopes = fragment._isotopes
        if isinstance(raw_isotopes, int):
            isotopes = _format_counts([("13C", raw_isotopes)]) if raw_isotopes else ""
        else:
            isotopes = _format_counts((raw_isotopes or {}).items())
        has_parent = bool(fragment.parent_sequence) and fragment.parent_sequence_length is not None
        records.append(
            {
                "ion_type": fragment.ion_type.value,
                "position": fragment.position,
                "charge_state": fragment.charge_state,
                "mz": fragment.mz,
                "mass": fragment.mass,
                "neutral_mass": fragment.neutral_mass,
                "monoisotopic": fragment.monoisotopic,
                "losses": _format_counts((fragment._losses or {}).items()),
                "isotopes": isotopes,
                "sequence": fragment.sequence if has_parent else None,
                "parent_sequence": fragment.parent_sequence if has_parent else None,
                "mzpaf": _mzpaf_or_none(fragment, has_parent),
            }
        )
    return records

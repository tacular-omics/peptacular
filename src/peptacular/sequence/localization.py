"""Functional API for modification localization isomers."""

from collections.abc import Iterable, Sequence
from typing import Any

from tacular import IonType, IonTypeLiteral, ToleranceUnit

from ..annotation import Fragment, ProFormaAnnotation
from ..annotation.localization import DEFAULT_MAX_ISOMERS
from ..annotation.localization import candidate_sites as _candidate_sites
from ..annotation.localization import localization_isomers as _localization_isomers
from ..annotation.localization import pairwise_site_determining_ions as _pairwise_site_determining_ions
from ..annotation.localization import site_determining_ions as _site_determining_ions
from .util import HasSequence, get_annotation_input

__all__ = [
    "localization_isomers",
    "candidate_sites",
    "site_determining_ions",
    "pairwise_site_determining_ions",
]


def localization_isomers(peptide: str | ProFormaAnnotation | HasSequence, *, max_isomers: int | None = DEFAULT_MAX_ISOMERS) -> list[ProFormaAnnotation]:
    """Expand every ambiguous modification position of ``peptide`` into its concrete placements.

    Candidate positions come from the ProForma string alone:

    - a ``#label`` group (``PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE``): each residue carrying the label;
    - a range (``PEP(ST)[Phospho]IDE``): each residue inside the range;
    - an unknown-position mod with no range (``[Phospho]?PEPTIDE``, ``[Phospho]^2?PEPTIDE``):
      **any residue**, whatever its letter. ``^n`` puts ``n`` copies on ``n`` different residues.

    A placed mod never goes on a residue that already carries one, or on a residue another
    ambiguity used in the same isomer (one mod per residue, as in :func:`candidate_sites`).
    A group's label stays on the placed mod with the score of the residue it landed on
    (``PEPS[Phospho#g1(0.8)]TIDE``); that mod string is where the score is kept. Identical
    isomers are returned once. The order is fixed: groups, then ranges, then unknown-position
    mods, each placed N- to C-terminal, the last one varying fastest.

    :param peptide: A ProForma string, an annotation, or an object with a ``.sequence`` string.
    :param max_isomers: Raise :class:`PeptacularError` if there would be more isomers than this.
        Defaults to 10,000; pass ``None`` for no limit.
    :return: The isomers as annotations (call ``.serialize()`` for strings).

    >>> import peptacular as pt
    >>> [a.serialize() for a in pt.localization_isomers("[Phospho]?PEST")]
    ['P[Phospho]EST', 'PE[Phospho]ST', 'PES[Phospho]T', 'PEST[Phospho]']
    >>> [a.serialize() for a in pt.localization_isomers("PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE")]
    ['PEPS[Phospho#g1(0.8)]TIDE', 'PEPST[Phospho#g1(0.2)]IDE']
    """
    return _localization_isomers(get_annotation_input(peptide, copy=False), max_isomers=max_isomers)


def candidate_sites(peptide: str | ProFormaAnnotation | HasSequence, mod: Any, *, residues: str) -> list[tuple[int, ProFormaAnnotation]]:
    """Place ``mod`` on each unmodified residue of ``peptide`` whose letter is in ``residues``.

    ``residues`` is required: peptacular ships no site table, so you say where the mod can go.

    :param peptide: A ProForma string, an annotation, or an object with a ``.sequence`` string.
    :param mod: The modification (``"Phospho"``, ``79.966``, ...).
    :param residues: One-letter codes that can carry ``mod``, e.g. ``"STY"``.
    :return: ``(position, isomer)`` pairs in sequence order, position 0-based.

    >>> import peptacular as pt
    >>> [(i, a.serialize()) for i, a in pt.candidate_sites("PEPS[Phospho]TYK", "Phospho", residues="STY")]
    [(4, 'PEPS[Phospho]T[Phospho]YK'), (5, 'PEPS[Phospho]TY[Phospho]K')]
    """
    return _candidate_sites(get_annotation_input(peptide, copy=False), mod, residues=residues)


def site_determining_ions(
    isomers: Iterable[str | ProFormaAnnotation | HasSequence],
    *,
    ion_types: Sequence[IonType | IonTypeLiteral] = (IonType.B, IonType.Y),
    charges: Sequence[int] = (1,),
    tolerance: float | None = None,
    unit: ToleranceUnit = "da",
) -> list[list[Fragment]]:
    """Per isomer, the fragment ions whose m/z no other isomer can explain.

    An ion is site-determining for its isomer when its m/z is more than ``tolerance`` away from
    every ion of **every** other isomer (same ``ion_types`` and ``charges``). The ions are the
    :class:`Fragment` objects ``fragment()`` returns, in the same order. With three or more
    adjacent candidate sites the middle isomers get empty lists; use
    :func:`pairwise_site_determining_ions` to compare isomers two at a time.

    :param isomers: Candidate isomers, e.g. ``pt.localization_isomers(...)``.
    :param ion_types: Ion types to generate (default b and y).
    :param charges: Fragment charges (default 1).
    :param tolerance: Match tolerance; None compares m/z exactly (within 1e-6 Da). Window edges count as a match.
    :param unit: ``"da"`` or ``"ppm"``.
    :return: One list of fragments per isomer, in input order.

    >>> import peptacular as pt
    >>> ions = pt.site_determining_ions(pt.localization_isomers("PEP(ST)[Phospho]IDE"))
    >>> [[f"{f.ion_type}{f.position}" for f in per_isomer] for per_isomer in ions]
    [['b4', 'y4'], ['b4', 'y4']]
    """
    annotations = [get_annotation_input(isomer, copy=False) for isomer in isomers]
    return _site_determining_ions(annotations, ion_types=ion_types, charges=charges, tolerance=tolerance, unit=unit)


def pairwise_site_determining_ions(
    isomers: Iterable[str | ProFormaAnnotation | HasSequence],
    *,
    ion_types: Sequence[IonType | IonTypeLiteral] = (IonType.B, IonType.Y),
    charges: Sequence[int] = (1,),
    tolerance: float | None = None,
    unit: ToleranceUnit = "da",
) -> dict[tuple[int, int], list[Fragment]]:
    """For each ordered pair of isomers ``(i, j)``, the ions of ``i`` that ``j`` cannot explain.

    ``result[(i, j)]`` lists the fragments of isomer ``i`` whose m/z is more than ``tolerance``
    away from every ion of isomer ``j``: evidence for ``i`` over ``j``, the comparison Ascore
    and PhosphoRS score. Every ordered pair with ``i != j`` is a key.

    :param isomers: Candidate isomers, e.g. ``pt.localization_isomers(...)``.
    :param ion_types: Ion types to generate (default b and y).
    :param charges: Fragment charges (default 1).
    :param tolerance: Match tolerance; None compares m/z exactly (within 1e-6 Da). Window edges count as a match.
    :param unit: ``"da"`` or ``"ppm"``.
    :return: ``{(i, j): [Fragment, ...]}``, indices in input order.

    >>> import peptacular as pt
    >>> isomers = pt.localization_isomers("PEP(STY)[Phospho]IDEK")
    >>> pairs = pt.pairwise_site_determining_ions(isomers)
    >>> [f"{f.ion_type}{f.position}" for f in pairs[(1, 0)]], [f"{f.ion_type}{f.position}" for f in pairs[(1, 2)]]
    (['b4', 'y6'], ['b5', 'y5'])
    """
    annotations = [get_annotation_input(isomer, copy=False) for isomer in isomers]
    return _pairwise_site_determining_ions(annotations, ion_types=ion_types, charges=charges, tolerance=tolerance, unit=unit)

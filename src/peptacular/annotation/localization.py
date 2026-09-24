"""Expand ambiguous modification localisation into concrete isomers.

Candidate positions come only from the ProForma sequence itself: the residues a
``#label`` group names, the residues inside a range, or (for an unknown-position mod
with no range) every residue. No site data is looked up anywhere.
"""

from __future__ import annotations

import bisect
import itertools
import re
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

from tacular import IonType, IonTypeLiteral

from ..diagnostics import PeptacularError, UnsupportedOperationError

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation
    from .frag import Fragment

__all__ = [
    "localization_isomers",
    "candidate_sites",
    "site_determining_ions",
]

# ``#label`` or ``#label(score)`` inside a modification string.
_GROUP_TAG = re.compile(r"#([A-Za-z0-9_]+)(?:\(([^()]*)\))?")
# Default tolerance when none is given: only numerically identical m/z values match.
_EXACT_TOLERANCE_DA = 1e-6


def _is_cross_link_label(label: str) -> bool:
    upper = label.upper()
    return upper.startswith("XL") or upper == "BRANCH"


def _group_tag(mod_str: str) -> re.Match[str] | None:
    """Return the ``#label(score)`` match of a localisation group, or None (no label or a cross-link)."""
    match = _GROUP_TAG.search(mod_str)
    if match is None or _is_cross_link_label(match.group(1)):
        return None
    return match


@dataclass
class _Group:
    label: str
    mod: str | None = None
    count: int = 0
    positions: list[int] | None = None
    scores: dict[int, str] | None = None


@dataclass(frozen=True)
class _Placement:
    """Place ``count`` copies of a mod on distinct positions chosen from ``positions``."""

    positions: tuple[int, ...]
    count: int
    mod: str | None  # None: a group mod, built per position from ``group``
    group: _Group | None = None

    def mod_at(self, position: int) -> str:
        if self.group is None:
            assert self.mod is not None
            return self.mod
        assert self.group.mod is not None and self.group.scores is not None
        score = self.group.scores.get(position)
        suffix = f"({score})" if score is not None else ""
        return f"{self.group.mod}#{self.group.label}{suffix}"

    def choices(self) -> Iterator[tuple[int, ...]]:
        return itertools.combinations(self.positions, self.count)


def _reject_group_label(mods: dict[str, int] | None, where: str) -> None:
    if not mods:
        return
    for mod_str in mods:
        match = _group_tag(mod_str)
        if match is not None:
            raise UnsupportedOperationError(
                f"Localisation group #{match.group(1)} on {where} is not supported; groups are expanded only over residue positions."
            )


def _split_ambiguity(annotation: ProFormaAnnotation) -> tuple[ProFormaAnnotation, list[_Placement]]:
    """Return the annotation with every ambiguity removed, plus the placements to make."""
    base = annotation.copy()
    placements: list[_Placement] = []

    _reject_group_label(base._nterm_mods, "the N-terminus")
    _reject_group_label(base._cterm_mods, "the C-terminus")
    _reject_group_label(base._labile_mods, "a labile modification")
    _reject_group_label(base._unknown_mods, "an unknown-position modification")
    for interval in base.intervals:
        _reject_group_label(dict(interval.mods._mods or {}), "a range")

    # 1. #label groups on residues: collect members, then strip them from the residues.
    groups: dict[str, _Group] = {}
    internal = base._internal_mods or {}
    for position in sorted(internal):
        for mod_str, count in list(internal[position].items()):
            match = _group_tag(mod_str)
            if match is None:
                continue
            label, score = match.group(1), match.group(2)
            group = groups.setdefault(label, _Group(label=label, positions=[], scores={}))
            assert group.positions is not None and group.scores is not None
            mod_text = (mod_str[: match.start()] + mod_str[match.end() :]).strip()
            if mod_text:
                if group.mod is not None:
                    raise PeptacularError(f"Localisation group #{label} has more than one modification: {group.mod!r} and {mod_text!r}")
                group.mod = mod_text
                group.count = count
            if position not in group.positions:
                group.positions.append(position)
            if score is not None:
                group.scores[position] = score
            del internal[position][mod_str]
        if not internal[position]:
            del internal[position]
    for group in groups.values():
        if group.mod is None:
            raise PeptacularError(f"Localisation group #{group.label} has no modification, only #{group.label} references")
        assert group.positions is not None
        placements.append(_Placement(positions=tuple(group.positions), count=group.count, mod=None, group=group))

    # 2. Ranges: each range mod goes on residues inside the range.
    kept_intervals = []
    for interval in base.intervals:
        for mod_str, count in (interval.mods._mods or {}).items():
            placements.append(_Placement(positions=tuple(range(interval.start, interval.end)), count=count, mod=mod_str))
        if interval.ambiguous:
            kept_intervals.append(interval.update(mods=None))

    # 3. Unknown-position mods: any residue.
    all_positions = tuple(range(len(base.sequence)))
    for mod_str, count in (base._unknown_mods or {}).items():
        placements.append(_Placement(positions=all_positions, count=count, mod=mod_str))

    base._internal_mods = internal or None
    base.set_intervals(kept_intervals or None)
    base.clear_unknown_mods()

    for placement in placements:
        if placement.count > len(placement.positions):
            raise PeptacularError(f"Cannot place {placement.count} copies of a modification on {len(placement.positions)} candidate residues")
    return base, placements


def _canonical_key(annotation: ProFormaAnnotation) -> tuple[Any, ...]:
    internal = annotation._internal_mods or {}
    return tuple(sorted((position, tuple(sorted(mods.items()))) for position, mods in internal.items()))


def _iter_choices(placements: Sequence[_Placement]) -> Iterator[tuple[tuple[int, ...], ...]]:
    """Lazy ``itertools.product`` over each placement's choices, so ``max_isomers`` can stop a huge expansion early."""
    if not placements:
        yield ()
        return
    for first in placements[0].choices():
        for rest in _iter_choices(placements[1:]):
            yield (first, *rest)


def _iter_isomers(base: ProFormaAnnotation, placements: Sequence[_Placement]) -> Iterator[ProFormaAnnotation]:
    for choice in _iter_choices(placements):
        isomer = base.copy()
        for placement, positions in zip(placements, choice, strict=True):
            for position in positions:
                isomer.append_internal_mod_at_index(position, placement.mod_at(position), inplace=True)
        yield isomer


def localization_isomers(annotation: ProFormaAnnotation, *, max_isomers: int | None = None) -> list[ProFormaAnnotation]:
    """Expand every ambiguous modification position into its concrete placements.

    Three kinds of ambiguity are expanded, and positions come from the ProForma string only:

    - A ``#label`` group (``PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE``): the mod goes on each
      residue that carries the label.
    - A range (``PEP(ST)[Phospho]IDE``): the mod goes on each residue inside the range.
    - An unknown-position mod (``[Phospho]?PEPTIDE``, ``[Phospho]^2?PEPTIDE``): with no range
      to limit it, the mod can go on **any residue**, whatever its letter. A count of ``n``
      places ``n`` copies on ``n`` different residues. Termini are not candidates.

    Every combination of placements is returned. Isomers that come out identical (two copies
    of the same mod swapped between ranges, say) are returned once.

    **Group scores.** A group's label stays on the placed mod together with the score of the
    residue it landed on, so ``PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE`` gives
    ``PEPS[Phospho#g1(0.8)]TIDE`` and ``PEPST[Phospho#g1(0.2)]IDE``. Read the score back from the
    mod string; a residue with no score gets ``Phospho#g1``. Ranges and unknown-position mods
    carry no score.

    **Order.** Deterministic: ambiguities are taken in the order groups (by first residue),
    ranges, then unknown-position mods; each one's placements run N- to C-terminal, and the
    last ambiguity varies fastest. An annotation with no ambiguity returns a one-item list
    holding a copy of it.

    Cross-link (``#XL1``) and branch (``#BRANCH``) labels are left alone. An ambiguous-order
    range (``(?ST)``) keeps its range after its mods are placed.

    :param annotation: The annotation to expand. It is not modified.
    :param max_isomers: Raise if there would be more than this many isomers.
    :return: The isomers, each with only localised residue mods where the ambiguity was.
    :raises PeptacularError: If ``max_isomers`` is exceeded or not a positive int, or a group is
        malformed (no modification, or two different ones).
    :raises UnsupportedOperationError: If a group label is on a terminus, a range, a labile mod
        or an unknown-position mod.

    >>> import peptacular as pt
    >>> [a.serialize() for a in pt.parse("PEP(ST)[Phospho]IDE").localization_isomers()]
    ['PEPS[Phospho]TIDE', 'PEPST[Phospho]IDE']
    """
    if max_isomers is not None and (isinstance(max_isomers, bool) or not isinstance(max_isomers, int) or max_isomers < 1):
        raise PeptacularError(f"max_isomers must be a positive int or None, got {max_isomers!r}")

    base, placements = _split_ambiguity(annotation)
    isomers: list[ProFormaAnnotation] = []
    seen: set[tuple[Any, ...]] = set()
    for isomer in _iter_isomers(base, placements):
        key = _canonical_key(isomer)
        if key in seen:
            continue
        seen.add(key)
        isomers.append(isomer)
        if max_isomers is not None and len(isomers) > max_isomers:
            raise PeptacularError(f"{annotation.serialize()!r} has more than max_isomers={max_isomers} localisation isomers")
    return isomers


def candidate_sites(annotation: ProFormaAnnotation, mod: Any, *, residues: str) -> list[tuple[int, ProFormaAnnotation]]:
    """Place ``mod`` on each unmodified residue whose letter is in ``residues``.

    ``residues`` is required: peptacular has no built-in site table, so say which residues
    the mod can sit on (``residues="STY"`` for phosphorylation). A residue that already
    carries a modification is skipped.

    :param annotation: The annotation to place the mod on. It is not modified.
    :param mod: The modification (a ProForma mod string such as ``"Phospho"``, a mass, ...).
    :param residues: One-letter codes of the residues that can carry ``mod``.
    :return: ``(position, isomer)`` pairs, position 0-based, in sequence order.
    :raises PeptacularError: If ``residues`` is empty or not a string of letters.

    >>> import peptacular as pt
    >>> [(i, a.serialize()) for i, a in pt.parse("PEPSTIDE").candidate_sites("Phospho", residues="ST")]
    [(3, 'PEPS[Phospho]TIDE'), (4, 'PEPST[Phospho]IDE')]
    """
    if not isinstance(residues, str) or not residues or not residues.isalpha():
        raise PeptacularError(f"residues must be a non-empty string of one-letter residue codes, got {residues!r}")
    wanted = set(residues.upper())
    sites: list[tuple[int, ProFormaAnnotation]] = []
    for position, residue in enumerate(annotation.sequence):
        if residue not in wanted or annotation.has_internal_mods_at_index(position):
            continue
        sites.append((position, annotation.append_internal_mod_at_index(position, mod, inplace=False)))
    return sites


def site_determining_ions(
    isomers: Iterable[ProFormaAnnotation],
    *,
    ion_types: Sequence[IonType | IonTypeLiteral] = (IonType.B, IonType.Y),
    charges: Sequence[int] = (1,),
    tolerance: float | None = None,
    unit: Literal["da", "ppm"] = "da",
) -> list[list[Fragment]]:
    """Per isomer, the fragment ions no other isomer can explain.

    An ion of one isomer is site-determining when its m/z is more than ``tolerance`` away
    from **every** ion (of the requested types and charges) of every other isomer. A peak at
    that m/z is then evidence for this isomer alone.

    Fragments come from :meth:`ProFormaAnnotation.fragment`, so the result holds the same
    :class:`Fragment` objects ``fragment()`` returns, in its order.

    :param isomers: The candidate isomers, e.g. from :func:`localization_isomers`.
    :param ion_types: Ion types to generate.
    :param charges: Fragment charges.
    :param tolerance: Match tolerance. None compares m/z values exactly (within 1e-6 Da).
    :param unit: ``"da"`` or ``"ppm"`` (relative to the ion's own m/z).
    :return: One list per isomer, in input order. With a single isomer every ion is returned.
    :raises PeptacularError: If ``tolerance`` is negative or ``unit`` is not ``"da"``/``"ppm"``.
    """
    if unit not in ("da", "ppm"):
        raise PeptacularError(f"unit must be 'da' or 'ppm', got {unit!r}")
    if tolerance is not None and (isinstance(tolerance, bool) or not tolerance >= 0):
        raise PeptacularError(f"tolerance must be a non-negative number or None, got {tolerance!r}")

    fragments = [isomer.fragment(ion_types=ion_types, charges=charges) for isomer in isomers]
    pool = sorted((fragment.mz, owner) for owner, frags in enumerate(fragments) for fragment in frags)
    pool_mz = [mz for mz, _ in pool]

    def window(mz: float) -> float:
        if tolerance is None:
            return _EXACT_TOLERANCE_DA
        return tolerance if unit == "da" else mz * tolerance * 1e-6

    result: list[list[Fragment]] = []
    for owner, frags in enumerate(fragments):
        determining: list[Fragment] = []
        for fragment in frags:
            mz = fragment.mz
            tol = window(mz)
            index = bisect.bisect_left(pool_mz, mz - tol)
            shared = False
            while index < len(pool) and pool_mz[index] <= mz + tol:
                if pool[index][1] != owner:
                    shared = True
                    break
                index += 1
            if not shared:
                determining.append(fragment)
        result.append(determining)
    return result

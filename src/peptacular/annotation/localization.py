"""Expand ambiguous modification localization into concrete isomers.

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
from typing import TYPE_CHECKING, Any

from tacular import IonType, IonTypeLiteral, ToleranceUnit, tolerance_window

from ..diagnostics import PeptacularError, UnsupportedOperationError

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation
    from .frag import Fragment

__all__ = [
    "DEFAULT_MAX_ISOMERS",
    "localization_isomers",
    "candidate_sites",
    "site_determining_ions",
    "pairwise_site_determining_ions",
]

DEFAULT_MAX_ISOMERS = 10_000
"""Default ``max_isomers`` of :func:`localization_isomers`. Pass ``max_isomers=None`` for no limit."""

# ``#label`` or ``#label(score)`` at the end of one ``|``-separated part of a modification string.
_GROUP_TAG = re.compile(r"#([A-Za-z0-9_]+)(?:\(([^()]*)\))?\s*$")
# Default tolerance when none is given: only numerically identical m/z values match.
_EXACT_TOLERANCE_DA = 1e-6


def _is_cross_link_label(label: str) -> bool:
    upper = label.upper()
    return upper.startswith("XL") or upper == "BRANCH"


@dataclass(frozen=True)
class _GroupTag:
    """A ``#label(score)`` tag found in a mod string, with the text on either side of it."""

    label: str
    score: str | None
    prefix: str
    suffix: str

    @property
    def mod_text(self) -> str:
        return (self.prefix + self.suffix).strip()


def _group_tag(mod_str: str) -> _GroupTag | None:
    """Return the localization-group tag of ``mod_str``, or None (no label, or a cross-link).

    The tag is looked for at the end of each ``|``-separated part of the mod string, skipping
    ``INFO:`` parts, so a ``#`` inside free text is not read as a group label.
    """
    offset = 0
    for part in mod_str.split("|"):
        if not part.strip().upper().startswith("INFO:"):
            match = _GROUP_TAG.search(part)
            if match is not None and not _is_cross_link_label(match.group(1)):
                return _GroupTag(
                    label=match.group(1),
                    score=match.group(2),
                    prefix=mod_str[: offset + match.start()],
                    suffix=mod_str[offset + match.end() :],
                )
        offset += len(part) + 1
    return None


@dataclass
class _Group:
    label: str
    tag: _GroupTag | None = None  # the member that carries the modification
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
        tag = self.group.tag
        assert tag is not None and self.group.scores is not None
        score = self.group.scores.get(position)
        suffix = f"({score})" if score is not None else ""
        return f"{tag.prefix.rstrip()}#{self.group.label}{suffix}{tag.suffix}"

    def choices(self, used: frozenset[int]) -> Iterator[tuple[int, ...]]:
        free = [position for position in self.positions if position not in used]
        return itertools.combinations(free, self.count)


def _reject_group_label(mods: dict[str, int] | None, where: str) -> None:
    if not mods:
        return
    for mod_str in mods:
        tag = _group_tag(mod_str)
        if tag is not None:
            raise UnsupportedOperationError(f"Localization group #{tag.label} on {where} is not supported; groups are expanded only over residue positions.")


def _copies(count: int) -> str:
    return "1 copy" if count == 1 else f"{count} copies"


def _split_ambiguity(annotation: ProFormaAnnotation) -> tuple[ProFormaAnnotation, list[_Placement]]:
    """Return the annotation with every ambiguity removed, plus the placements to make.

    Residues that still carry a modification after the ambiguity is removed are not candidates
    for any placement (the same rule as :func:`candidate_sites`).
    """
    base = annotation.copy()
    raw: list[tuple[tuple[int, ...], int, str | None, _Group | None]] = []

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
            tag = _group_tag(mod_str)
            if tag is None:
                continue
            group = groups.setdefault(tag.label, _Group(label=tag.label, positions=[], scores={}))
            assert group.positions is not None and group.scores is not None
            if tag.mod_text:
                if group.tag is not None:
                    if group.tag.mod_text == tag.mod_text:
                        raise PeptacularError(
                            f"Localization group #{tag.label} holds one modification, but {tag.mod_text!r} is written with "
                            f"#{tag.label} on more than one residue; give each copy its own label (#g1, #g2)"
                        )
                    raise PeptacularError(f"Localization group #{tag.label} has more than one modification: {group.tag.mod_text!r} and {tag.mod_text!r}")
                group.tag = tag
                group.count = count
            if position not in group.positions:
                group.positions.append(position)
            if tag.score is not None:
                group.scores[position] = tag.score
            del internal[position][mod_str]
        if not internal[position]:
            del internal[position]
    for group in groups.values():
        if group.tag is None:
            raise PeptacularError(f"Localization group #{group.label} has no modification, only #{group.label} references")
        assert group.positions is not None
        raw.append((tuple(group.positions), group.count, None, group))

    # 2. Ranges: each range mod goes on residues inside the range. A range with no mods is dropped.
    kept_intervals = []
    for interval in base.intervals:
        for mod_str, count in (interval.mods._mods or {}).items():
            raw.append((tuple(range(interval.start, interval.end)), count, mod_str, None))
        if interval.ambiguous:
            kept_intervals.append(interval.update(mods=None))

    # 3. Unknown-position mods: any residue.
    all_positions = tuple(range(len(base.sequence)))
    for mod_str, count in (base._unknown_mods or {}).items():
        raw.append((all_positions, count, mod_str, None))

    base._internal_mods = internal or None
    base.set_intervals(kept_intervals or None)
    base.clear_unknown_mods()

    occupied = set(internal)
    placements: list[_Placement] = []
    for positions, count, mod, group in raw:
        if group is not None:
            taken = sorted(position for position in positions if position in occupied)
            if taken:
                raise PeptacularError(
                    f"Group #{group.label} lists residue {taken[0] + 1}, {base.sequence[taken[0]]} at position {taken[0]} "
                    "(0-based), which already carries a "
                    "modification; one mod per residue, so that placement is impossible. Remove the group tag from "
                    "that residue or drop its other mod."
                )
        free = tuple(position for position in positions if position not in occupied)
        if count > len(free):
            raise PeptacularError(f"Cannot place {_copies(count)} of a modification on {len(free)} unmodified candidate residue{'' if len(free) == 1 else 's'}")
        placements.append(_Placement(positions=free, count=count, mod=mod, group=group))
    return base, placements


def _canonical_key(annotation: ProFormaAnnotation) -> tuple[Any, ...]:
    internal = annotation._internal_mods or {}
    return tuple(sorted((position, tuple(sorted(mods.items()))) for position, mods in internal.items()))


def _iter_choices(placements: Sequence[_Placement], used: frozenset[int] = frozenset()) -> Iterator[tuple[tuple[int, ...], ...]]:
    """Lazy product over each placement's choices, never reusing a residue, so ``max_isomers`` can stop a huge expansion early."""
    if not placements:
        yield ()
        return
    for first in placements[0].choices(used):
        for rest in _iter_choices(placements[1:], used | frozenset(first)):
            yield (first, *rest)


def _iter_isomers(base: ProFormaAnnotation, placements: Sequence[_Placement]) -> Iterator[ProFormaAnnotation]:
    for choice in _iter_choices(placements):
        isomer = base.copy()
        for placement, positions in zip(placements, choice, strict=True):
            for position in positions:
                isomer.append_internal_mod_at_index(position, placement.mod_at(position), inplace=True)
        yield isomer


def localization_isomers(annotation: ProFormaAnnotation, *, max_isomers: int | None = DEFAULT_MAX_ISOMERS) -> list[ProFormaAnnotation]:
    """Expand every ambiguous modification position into its concrete placements.

    Three kinds of ambiguity are expanded, and positions come from the ProForma string only:

    - A ``#label`` group (``PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE``): the mod goes on each
      residue that carries the label.
    - A range (``PEP(ST)[Phospho]IDE``): the mod goes on each residue inside the range. A
      range with no modification (``P(ES)T``) has nothing to place and is dropped.
    - An unknown-position mod (``[Phospho]?PEPTIDE``, ``[Phospho]^2?PEPTIDE``): with no range
      to limit it, the mod can go on **any residue**, whatever its letter. A count of ``n``
      places ``n`` copies on ``n`` different residues. Termini are not candidates.

    **One mod per residue.** A placed mod never goes on a residue that already carries a
    modification, or on a residue another ambiguity placed a mod on in the same isomer, so
    ``[Phospho]?PES[Phospho]T`` never gives ``PES[Phospho][Phospho]T``. This is the same rule
    as :func:`candidate_sites`. A ``#label`` group that lists a residue which already carries
    another modification (``PS[Oxidation][Phospho#g1]T[#g1]``) raises
    :class:`~peptacular.PeptacularError` rather than dropping the placement the string names.

    Every combination of placements is returned. Isomers that come out identical (two copies
    of the same mod swapped between ranges, say) are returned once.

    **Group scores.** A group's label stays on the placed mod together with the score of the
    residue it landed on, so ``PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE`` gives
    ``PEPS[Phospho#g1(0.8)]TIDE`` and ``PEPST[Phospho#g1(0.2)]IDE``. Read the score back from the
    mod string; a residue with no score gets ``Phospho#g1``. Ranges and unknown-position mods
    carry no score. A ``#`` inside an ``INFO:`` tag is text, not a group label.

    **Order.** Deterministic: ambiguities are taken in the order groups (by first residue),
    ranges, then unknown-position mods; each one's placements run N- to C-terminal, and the
    last ambiguity varies fastest. An annotation with no ambiguity returns a one-item list
    holding a copy of it.

    Cross-link (``#XL1``) and branch (``#BRANCH``) labels are left alone. An ambiguous-order
    range (``(?ST)``) keeps its range after its mods are placed.

    :param annotation: The annotation to expand. It is not modified.
    :param max_isomers: Raise if there would be more than this many isomers. Defaults to
        ``DEFAULT_MAX_ISOMERS`` (10,000); pass ``None`` for no limit.
    :return: The isomers, each with only localized residue mods where the ambiguity was.
    :raises PeptacularError: If ``max_isomers`` is exceeded or not a positive int, a group is
        malformed (no modification, or more than one), or the mods cannot all be put on
        different unmodified residues.
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
            raise PeptacularError(
                f"{annotation.serialize()!r} has more than max_isomers={max_isomers} localization isomers; "
                "pass a larger max_isomers= or max_isomers=None for no limit"
            )
    if not isomers:
        raise PeptacularError(f"{annotation.serialize()!r}: the ambiguous modifications cannot all be put on different unmodified residues")
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
    :raises PeptacularError: If ``residues`` is empty or not a string of letters, or ``mod`` is
        not a valid modification.

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
        isomer = annotation.append_internal_mod_at_index(position, mod, inplace=False)
        if not sites and (error := isomer.get_internal_mods_at_index(position).validate()):
            raise PeptacularError(f"Invalid modification {mod!r}: {error}")
        sites.append((position, isomer))
    return sites


def _check_tolerance(tolerance: float | None, unit: str) -> None:
    if unit not in ("da", "ppm"):
        raise PeptacularError(f"unit must be 'da' or 'ppm', got {unit!r}")
    if tolerance is not None and (isinstance(tolerance, bool) or not tolerance >= 0):
        raise PeptacularError(f"tolerance must be a non-negative number or None, got {tolerance!r}")


def _window(mz: float, tolerance: float | None, unit: ToleranceUnit) -> tuple[float, float]:
    if tolerance is None:
        return tolerance_window(mz, _EXACT_TOLERANCE_DA, unit="da")
    return tolerance_window(mz, tolerance, unit=unit)


def _explained(mz: float, sorted_mz: Sequence[float], tolerance: float | None, unit: ToleranceUnit) -> bool:
    """Whether any value of ``sorted_mz`` lies in the tolerance window around ``mz`` (edges included)."""
    lo, hi = _window(mz, tolerance, unit)
    index = bisect.bisect_left(sorted_mz, lo)
    return index < len(sorted_mz) and sorted_mz[index] <= hi


def _isomer_fragments(isomers: Iterable[ProFormaAnnotation], ion_types: Sequence[IonType | IonTypeLiteral], charges: Sequence[int]) -> list[list[Fragment]]:
    return [isomer.fragment(ion_types=ion_types, charges=charges) for isomer in isomers]


def site_determining_ions(
    isomers: Iterable[ProFormaAnnotation],
    *,
    ion_types: Sequence[IonType | IonTypeLiteral] = (IonType.B, IonType.Y),
    charges: Sequence[int] = (1,),
    tolerance: float | None = None,
    unit: ToleranceUnit = "da",
) -> list[list[Fragment]]:
    """Per isomer, the fragment ions **no other isomer** can explain.

    An ion of one isomer is site-determining here when its m/z is more than ``tolerance`` away
    from **every** ion (of the requested types and charges) of **every** other isomer. A peak
    at that m/z is then evidence for this isomer alone.

    With three or more candidate sites next to each other, a middle isomer shares each of its
    b/y ions with one neighbour or the other, so its list is empty. To tell isomers apart two
    at a time (as Ascore and PhosphoRS do), use :func:`pairwise_site_determining_ions`.

    Fragments come from :meth:`ProFormaAnnotation.fragment`, so the result holds the same
    :class:`Fragment` objects ``fragment()`` returns, in its order.

    :param isomers: The candidate isomers, e.g. from :func:`localization_isomers`.
    :param ion_types: Ion types to generate.
    :param charges: Fragment charges.
    :param tolerance: Match tolerance. None compares m/z values exactly (within 1e-6 Da). The
        window is :func:`tacular.tolerance_window` and its edges count as a match.
    :param unit: ``"da"`` or ``"ppm"`` (relative to the ion's own m/z).
    :return: One list per isomer, in input order. With a single isomer every ion is returned.
    :raises PeptacularError: If ``tolerance`` is negative or ``unit`` is not ``"da"``/``"ppm"``.
    """
    _check_tolerance(tolerance, unit)
    fragments = _isomer_fragments(isomers, ion_types, charges)
    pool = sorted((fragment.mz, owner) for owner, frags in enumerate(fragments) for fragment in frags)
    pool_mz = [mz for mz, _ in pool]

    result: list[list[Fragment]] = []
    for owner, frags in enumerate(fragments):
        determining: list[Fragment] = []
        for fragment in frags:
            lo, hi = _window(fragment.mz, tolerance, unit)
            index = bisect.bisect_left(pool_mz, lo)
            shared = False
            while index < len(pool) and pool_mz[index] <= hi:
                if pool[index][1] != owner:
                    shared = True
                    break
                index += 1
            if not shared:
                determining.append(fragment)
        result.append(determining)
    return result


def pairwise_site_determining_ions(
    isomers: Iterable[ProFormaAnnotation],
    *,
    ion_types: Sequence[IonType | IonTypeLiteral] = (IonType.B, IonType.Y),
    charges: Sequence[int] = (1,),
    tolerance: float | None = None,
    unit: ToleranceUnit = "da",
) -> dict[tuple[int, int], list[Fragment]]:
    """For each ordered pair of isomers ``(i, j)``, the ions of ``i`` that ``j`` cannot explain.

    ``result[(i, j)]`` holds the fragments of isomer ``i`` whose m/z is more than
    ``tolerance`` away from every ion of isomer ``j``. A peak at one of them is evidence for
    ``i`` over ``j``; this is the pairwise comparison Ascore and PhosphoRS score. Unlike
    :func:`site_determining_ions`, a middle isomer of three adjacent sites still has ions
    against each neighbour.

    :param isomers: The candidate isomers, e.g. from :func:`localization_isomers`.
    :param ion_types: Ion types to generate.
    :param charges: Fragment charges.
    :param tolerance: Match tolerance. None compares m/z values exactly (within 1e-6 Da).
    :param unit: ``"da"`` or ``"ppm"`` (relative to the ion's own m/z).
    :return: A dict with a key for every ordered pair ``(i, j)``, ``i != j``, indices in input
        order, each value in ``fragment()`` order. One isomer gives an empty dict.
    :raises PeptacularError: If ``tolerance`` is negative or ``unit`` is not ``"da"``/``"ppm"``.
    """
    _check_tolerance(tolerance, unit)
    fragments = _isomer_fragments(isomers, ion_types, charges)
    sorted_mz = [sorted(fragment.mz for fragment in frags) for frags in fragments]
    result: dict[tuple[int, int], list[Fragment]] = {}
    for i, frags in enumerate(fragments):
        for j, other in enumerate(sorted_mz):
            if i != j:
                result[(i, j)] = [fragment for fragment in frags if not _explained(fragment.mz, other, tolerance, unit)]
    return result

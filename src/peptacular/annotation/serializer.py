from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


# Bracket constants
_MOD_BRACKETS = ("[", "]")
_LABILE_BRACKETS = ("{", "}")
_STATIC_BRACKETS = ("<", ">")
_ISOTOPE_BRACKETS = ("<", ">")
_UNKNOWN_BRACKETS = ("[", "]")

# Special characters
_TERM_MOD = "-"
_UNKNOWN = "?"
_INTERVAL_START = "("
_INTERVAL_END = ")"
_CHARGE_SEP = "/"


def _serialize_mod_dict(
    mods: dict[str, int] | None,
    brackets: tuple[str, str],
    comps: list[str],
    allow_multiplier: bool,
) -> None:
    """Serialize a modification dictionary into comps list."""
    if mods is None or len(mods) == 0:
        return

    for mod_str, count in mods.items():
        if count == 1:
            comps.append(f"{brackets[0]}{mod_str}{brackets[1]}")
        else:
            if allow_multiplier is False:
                s = f"{brackets[0]}{mod_str}{brackets[1]}"
                for _ in range(count):
                    comps.append(s)
            else:
                comps.append(f"{brackets[0]}{mod_str}{brackets[1]}^{count}")


def _serialize_labile(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize labile modifications into comps list."""
    comps.append(annotation.labile_mods.serialize())


def _serialize_static(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize static modifications into comps list."""
    comps.append(annotation.static_mods.serialize())


def _serialize_isotope(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize isotope modifications into comps list."""
    comps.append(annotation.isotope_mods.serialize())


def _serialize_unknown(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize unknown modifications into comps list."""
    comps.append(annotation.unknown_mods.serialize())


def _serialize_nterm(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize N-terminal modifications into comps list."""
    comps.append(annotation.nterm_mods.serialize())


def _serialize_cterm(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize C-terminal modifications into comps list."""
    comps.append(annotation.cterm_mods.serialize())


def _serialize_start(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize the start portion into comps list."""
    # Add names if present
    if annotation.has_compound_name:
        comps.append(f"(>>>{annotation.compound_name})")

    # Add modifications in proper order
    if annotation.has_isotope_mods:
        _serialize_isotope(annotation, comps)

    if annotation.has_static_mods:
        _serialize_static(annotation, comps)

    if annotation.has_ion_name:
        comps.append(f"(>>{annotation.ion_name})")

    if annotation.has_peptide_name:
        comps.append(f"(>{annotation.peptide_name})")

    if annotation.has_labile_mods:
        _serialize_labile(annotation, comps)

    if annotation.has_unknown_mods:
        _serialize_unknown(annotation, comps)

    if annotation.has_nterm_mods:
        _serialize_nterm(annotation, comps)


def _serialize_middle(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize the middle portion (sequence and internal mods) into comps list."""
    if not annotation.has_sequence:
        return

    # FAST PATH: No internal mods or intervals - just add the sequence
    if not annotation.has_internal_mods and not annotation.has_intervals:
        comps.append(annotation.sequence)
        return

    # COMPLEX PATH: Has modifications or intervals
    seq = annotation.sequence

    internal_mods = annotation.internal_mods
    intervals = annotation.intervals

    # Main sequence loop
    for i, aa in enumerate(seq):
        # Check for interval starts at this position
        if annotation.has_intervals:
            for interval in annotation.intervals:
                if interval.start == i:
                    comps.append(_INTERVAL_START)
                    if interval.ambiguous:
                        comps.append(_UNKNOWN)

        # Add amino acid
        comps.append(aa)

        # Add internal mods at this position
        if i in internal_mods:
            comps.append(internal_mods[i].serialize())

        # Check for interval ends at this position
        for interval in intervals:
            if interval.end == i + 1:  # Intervals end after the position
                comps.append(_INTERVAL_END)
                # Add interval mods
                if interval.has_mods:
                    comps.append(interval.mods.serialize())


def _serialize_charge(
    annotation: ProFormaAnnotation,
    comps: list[str],
) -> None:
    """Serialize charge and charge adducts into comps list."""
    charge_type = annotation.charge_type
    if charge_type == "none":
        return
    if charge_type == "int":
        c = annotation.charge_state
        if c == 0:
            return
        comps.append(f"{_CHARGE_SEP}{annotation.charge_state}")
    if charge_type == "adducts":
        comps.append(f"{_CHARGE_SEP}{annotation.charge_adducts.serialize()}")


def serialize_charge(
    annotation: ProFormaAnnotation,
) -> str:
    """Serialize only the charge and charge adducts of a ProForma annotation."""
    comps: list[str] = []
    _serialize_charge(annotation, comps)
    return "".join(comps)


def _serialize_end(
    annotation: ProFormaAnnotation,
    comps: list[str],
    exclude_charge: bool = False,
) -> None:
    """Serialize the end portion into comps list."""
    # C-term mods
    if annotation.has_cterm_mods:
        _serialize_cterm(annotation, comps)

    # Charge
    if not exclude_charge:
        _serialize_charge(annotation, comps)


def serialize_annotation(
    annotation: ProFormaAnnotation,
    exclude_charge: bool = False,
) -> str:
    """Serialize a complete ProForma annotation."""
    comps: list[str] = []

    # Start
    _serialize_start(annotation, comps)

    # Middle
    _serialize_middle(annotation, comps)

    # End
    _serialize_end(annotation, comps, exclude_charge=exclude_charge)

    return "".join(comps)

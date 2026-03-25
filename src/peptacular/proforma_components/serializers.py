"""
Serializer functions for converting component objects to strings.

This module contains all serialization logic (ProForma notation output).
"""

import sys
from functools import lru_cache
from typing import TYPE_CHECKING

from tacular import (
    Element,
    Monosaccharide,
)

from ..constants import (
    CV_TO_ACCESSION_PREFIX,
    CV_TO_MASS_PREFIX,
    CV_TO_NAME_PREFIX,
    Terminal,
)

if TYPE_CHECKING:
    from .comps import (
        MODIFICATION_TAG_TYPE,
        MODIFICATION_TYPE,
        ChargedFormula,
        CompoundPeptidoformIon,
        FixedModification,
        FormulaElement,
        GlobalChargeCarrier,
        GlycanComponent,
        GlycanTag,
        IsotopeReplacement,
        ModificationAmbiguousPrimary,
        ModificationAmbiguousSecondary,
        ModificationCrossLinker,
        ModificationTags,
        Peptidoform,
        PeptidoformIon,
        PositionRule,
        SequenceElement,
        SequenceRegion,
        TagAccession,
        TagCustom,
        TagInfo,
        TagMass,
        TagName,
    )


@lru_cache(maxsize=1024)
def get_element_key(fe: "FormulaElement") -> str:
    """Get a unique key for the formula element, including isotope if present"""
    if fe.isotope is not None:
        return sys.intern(f"[{fe.isotope}{fe.element}]")
    else:
        return sys.intern(f"{fe.element}")


def serialize_formula_element(fe: "FormulaElement") -> str:
    """
    Cached string representation of the formula element in ProForma notation.

    Follows ProForma 2.0 specification:
    - When isotope is specified, count MUST be inside brackets: [13C2]
    - When no isotope, count is after element: C2

    Args:
        fe: FormulaElement to serialize

    Returns:
        String representation like 'C2', '[13C2]', or 'H'
    """
    if fe.isotope is not None:
        # Format: [13C2] - count must be inside bracket
        if fe.occurance != 1:
            return sys.intern(f"[{fe.isotope}{fe.element}{fe.occurance}]")
        else:
            return sys.intern(f"[{fe.isotope}{fe.element}]")
    else:
        # Format: C2 - count after element
        if fe.occurance != 1:
            return sys.intern(f"{fe.element}{fe.occurance}")
        else:
            return sys.intern(f"{fe.element}")


@lru_cache(maxsize=1024)
def serialize_charged_formula(
    cf: "ChargedFormula",
    space: str = "",
    hill_order: bool = False,
    include_formula_prefix: bool = True,
) -> str:
    """
    Cached string representation of the charged formula in ProForma notation.

    Args:
        cf: ChargedFormula to serialize
        space: Separator between formula elements (default: no space)
        hill_order: If True, sort elements in Hill system order (C, H, then alphabetical)

    Returns:
        String representation like 'Formula:C2H6' or 'Formula:C2H6:z+2'
    """
    if hill_order:
        # Hill system: C first, H second, then alphabetical
        def hill_sort_key(fe: "FormulaElement") -> tuple[int, int | str, int | str]:
            # Get the element as a string
            element_str = fe.element.value

            # Handle isotopes - they come after their natural element
            if fe.isotope is not None:
                # C comes first, then isotopes of C
                if element_str == "C":
                    return (0, fe.isotope, element_str)
                # H comes second, then isotopes of H
                elif element_str == "H":
                    return (1, fe.isotope, element_str)
                # Others alphabetically, then by isotope number
                else:
                    return (2, element_str, fe.isotope)
            else:
                # Natural elements: C, H, then alphabetical
                if element_str == "C":
                    return (0, 0, element_str)
                elif element_str == "H":
                    return (1, 0, element_str)
                else:
                    return (2, element_str, 0)

        sorted_elements = sorted(cf.formula, key=hill_sort_key)
        formula_str = space.join(serialize_formula_element(fe) for fe in sorted_elements)
    else:
        formula_str = space.join(serialize_formula_element(fe) for fe in cf.formula)

    prefix = "Formula:" if include_formula_prefix else ""

    if cf.charge is not None:
        charge_str = f":z{cf.charge:+}"
        return sys.intern(f"{prefix}{formula_str}{charge_str}{cf.serialize_position_score()}")
    else:
        return sys.intern(f"{prefix}{formula_str}{cf.serialize_position_score()}")


@lru_cache(maxsize=512)
def serialize_position_rule(pr: "PositionRule") -> str:
    """
    Cached string representation of the position rule in ProForma notation.

    Args:
        pr: PositionRule to serialize

    Returns:
        String representation like 'N-term', 'C-term:K', or 'M'
    """
    if pr.terminal == Terminal.ANYWHERE:
        if pr.amino_acid is not None:
            return sys.intern(pr.amino_acid.value)
        raise ValueError("Amino acid must be specified for ANYWHERE position rule")

    elif pr.amino_acid is not None:
        return sys.intern(f"{pr.terminal.value}:{pr.amino_acid.value}")
    else:
        return sys.intern(f"{pr.terminal.value}")


@lru_cache(maxsize=512)
def serialize_tag_accession(ta: "TagAccession") -> str:
    """
    Cached string representation of the tag accession.

    Args:
        ta: TagAccession to serialize

    Returns:
        String representation like 'UNIMOD:35'
    """
    return sys.intern(f"{CV_TO_ACCESSION_PREFIX[ta.cv]}{ta.accession}{ta.serialize_position_score()}")


@lru_cache(maxsize=512)
def serialize_tag_mass(tm: "TagMass") -> str:
    """
    Cached string representation of the tag mass.

    Args:
        tm: TagMass to serialize

    Returns:
        String representation like '+15.995' or 'UNIMOD:+15.995'
    """
    mass = tm.mass
    mass_str = f"{int(mass):+}" if mass == int(mass) else f"{mass:+}"
    if tm.cv is not None:
        return sys.intern(f"{CV_TO_MASS_PREFIX[tm.cv]}{mass_str}{tm.serialize_position_score()}")
    else:
        return sys.intern(f"{mass_str}{tm.serialize_position_score()}")


@lru_cache(maxsize=512)
def serialize_tag_name(tn: "TagName") -> str:
    """
    Cached string representation of the tag name.

    Args:
        tn: TagName to serialize

    Returns:
        String representation like 'Oxidation' or 'U:Oxidation'
    """
    if tn.cv is not None:
        return sys.intern(f"{CV_TO_NAME_PREFIX[tn.cv]}{tn.name}{tn.serialize_position_score()}")
    else:
        return sys.intern(f"{tn.name}{tn.serialize_position_score()}")


@lru_cache(maxsize=256)
def serialize_tag_info(ti: "TagInfo") -> str:
    """
    Cached string representation of the tag info.

    Args:
        ti: TagInfo to serialize

    Returns:
        String representation like 'INFO:some text'
    """
    return sys.intern(f"INFO:{ti.info}")


@lru_cache(maxsize=256)
def serialize_tag_custom(tc: "TagCustom") -> str:
    """
    Cached string representation of the tag custom.

    Args:
        tc: TagCustom to serialize

    Returns:
        String representation like 'C:custom_name'
    """
    return sys.intern(f"C:{tc.name}{tc.serialize_position_score()}")


@lru_cache(maxsize=256)
def serialize_glycan_component(gc: "GlycanComponent") -> str:
    """
    Cached string representation of the glycan component.

    Args:
        gc: GlycanComponent to serialize

    Returns:
        String representation like 'Hex5' or 'HexNAc4'
    """
    if isinstance(gc.monosaccharide, Monosaccharide):
        mono_str = gc.monosaccharide.value
    else:
        # It's a ChargedFormula
        mono_str = f"({serialize_charged_formula(gc.monosaccharide)})"

    if gc.occurance == 1:
        return sys.intern(mono_str)
    else:
        return sys.intern(f"{mono_str}{gc.occurance}")


@lru_cache(maxsize=256)
def serialize_glycan_tag(gt: "GlycanTag") -> str:
    """
    Cached string representation of the glycan tag.

    Args:
        gt: GlycanTag to serialize

    Returns:
        String representation like 'Glycan:Hex5HexNAc4'
    """
    glycan_comps: list[str] = []
    for gc in gt.components:
        glycan_comps.append(serialize_glycan_component(gc))
    return sys.intern(f"Glycan:{''.join(glycan_comps)}{gt.serialize_position_score()}")


@lru_cache(maxsize=128)
def serialize_isotope_replacement(ir: "IsotopeReplacement") -> str:
    """
    Cached string representation of the isotope replacement.

    Args:
        ir: IsotopeReplacement to serialize

    Returns:
        String representation like '13C', '15N', or 'D'
    """

    # Special case: Deuterium can be written as D or 2H
    if ir.element == Element.H and ir.isotope == 2:
        return sys.intern("D")
    else:
        return sys.intern(f"{ir.isotope}{ir.element.value}")


@lru_cache(maxsize=256)
def serialize_global_charge_carrier(gcc: "GlobalChargeCarrier") -> str:
    """
    Cached string representation of the global charge carrier.

    Args:
        gcc: GlobalChargeCarrier to serialize

    Returns:
        String representation like 'Na:z+1' or 'H:z+1^2'
    """
    formula_str = serialize_charged_formula(gcc.charged_formula, include_formula_prefix=False)

    if gcc.occurance == 1.0:
        return sys.intern(formula_str)
    else:
        # Format occurance, remove trailing zeros for whole numbers
        if gcc.occurance == int(gcc.occurance):
            occ_str = f"^{int(gcc.occurance)}"
        else:
            occ_str = f"^{gcc.occurance:g}"
        return sys.intern(f"{formula_str}{occ_str}")


@lru_cache(maxsize=512)
def serialize_modification_ambiguous_primary(
    mod: "ModificationAmbiguousPrimary",
) -> str:
    """
    Serialize a ModificationAmbiguousPrimary to ProForma notation.
    Format: [Modification#label(score)|Position:...|Limit:...|CoMKP|CoMUP]
    """

    # Start with the modification tags
    tags_str = serialize_modification_tags(mod.tags)

    # Add the label
    result = tags_str + f"#{mod.label}"

    # Add score if present
    if mod.score is not None:
        result += f"({mod.score})"

    # Add position constraints if present
    if mod.position is not None:
        position_strs = [str(pr) for pr in mod.position]
        result += f"|Position:{','.join(position_strs)}"

    # Add limit if present
    if mod.limit is not None:
        result += f"|Limit:{mod.limit}"

    # Add CoMKP if present and True
    if mod.comkp:
        result += "|CoMKP"

    # Add CoMUP if present and True
    if mod.comup:
        result += "|CoMUP"

    return sys.intern(result)


@lru_cache(maxsize=512)
def serialize_modification_ambiguous_secondary(
    mod: "ModificationAmbiguousSecondary",
) -> str:
    """
    Serialize a ModificationAmbiguousSecondary to ProForma notation.
    Format: [#label(score)]
    """
    result = f"#{mod.label}"

    # Add score if present
    if mod.score is not None:
        result += f"({mod.score})"

    return sys.intern(result)


@lru_cache(maxsize=512)
def serialize_modification_cross_linker(mod: "ModificationCrossLinker") -> str:
    """
    Serialize a ModificationCrossLinker to ProForma notation.
    Three cases:
    1. Primary definition with tags and label: [XLMOD:02001#XL1]
    2. Secondary reference without tags: [#XL1]
    3. Branch: [MOD:00093#BRANCH] or [#BRANCH]
    """

    if mod.tags is not None:
        # Primary definition: has tags
        tags_str = serialize_modification_tags(mod.tags)

        # Add label if present
        if mod.label is not None:
            result = tags_str + (f"#XL{mod.label}" if mod.label != "BRANCH" else "#BRANCH")
        else:
            result = tags_str

        return sys.intern(result)
    else:
        # Secondary reference: no tags, just label
        if mod.label is not None:
            return sys.intern(f"#XL{mod.label}" if mod.label != "BRANCH" else "#BRANCH")
        else:
            # Edge case: no tags and no label (shouldn't happen)
            return sys.intern("")


@lru_cache(maxsize=256)
def serialize_fixed_modification(mod: "FixedModification") -> str:
    """
    Serialize a FixedModification to ProForma notation.
    Format: [mod]@position,position
    Example: [Oxidation]@M or [TMT6plex]@K,N-term
    """
    # Build modification string
    mod_str = serialize_modification_tags(mod.modifications)

    # Build position rules string
    if mod.position_rules:
        pos_strs = [str(pr) for pr in mod.position_rules]
        pos_str = ",".join(pos_strs)
        return sys.intern(f"[{mod_str}]@{pos_str}")
    else:
        # No position rules specified
        return sys.intern(f"[{mod_str}]")


def serialize_modification_tag(tag: "MODIFICATION_TAG_TYPE") -> str:
    """Serialize a single modification tag."""
    # Import here to avoid circular dependency at module level
    from .comps import (
        ChargedFormula,
        GlycanTag,
        TagAccession,
        TagCustom,
        TagInfo,
        TagMass,
        TagName,
    )

    match tag:
        case TagAccession() | TagMass() | TagName() | TagInfo() | TagCustom() | ChargedFormula():
            return str(tag)
        case GlycanTag():
            return serialize_glycan_tag(tag)
        case _:
            raise ValueError(f"Unsupported modification tag type: {type(tag)}")


def serialize_modification_tags(mod_tags: "ModificationTags") -> str:
    """Serialize ModificationTags to ProForma notation."""
    tags_list: list[str] = []
    for tag in mod_tags.tags:
        tags_list.append(serialize_modification_tag(tag))
    return sys.intern("|".join(tags_list))


def serialize_modification(mod: "MODIFICATION_TYPE") -> str:
    """Serialize a modification (could be ModificationTags or ambiguous/cross-linker)."""
    # Import here to avoid circular dependency at module level
    from .comps import (
        ModificationAmbiguousPrimary,
        ModificationAmbiguousSecondary,
        ModificationCrossLinker,
        ModificationTags,
    )

    match mod:
        case ModificationAmbiguousPrimary() | ModificationCrossLinker() | ModificationAmbiguousSecondary():
            return sys.intern(str(mod))
        case ModificationTags():
            return serialize_modification_tags(mod)
        case _:
            raise ValueError(f"Unsupported Modification type: {type(mod)}")


@lru_cache(maxsize=1024)
def serialize_sequence_element(se: "SequenceElement") -> str:
    """
    Serialize a SequenceElement to ProForma notation.

    Args:
        se: SequenceElement to serialize

    Returns:
        String representation like 'M[Oxidation]' or 'K'
    """

    result = se.amino_acid.value

    if se.modifications:
        result += "".join(f"[{serialize_modification(mod)}]" for mod in se.modifications)
    return sys.intern(result)


@lru_cache(maxsize=512)
def serialize_sequence_region(sr: "SequenceRegion") -> str:
    """
    Serialize a SequenceRegion to ProForma notation.

    Args:
        sr: SequenceRegion to serialize

    Returns:
        String representation like '(PEPTIDE)[Oxidation]'
    """

    # Format: (SEQUENCE)[mod][mod]...
    seq_str = "".join(str(se) for se in sr.sequence)
    result = f"({'' if not sr.ambiguous else '?'}{seq_str})"

    if sr.modifications:
        for mod in sr.modifications:
            result += f"[{serialize_modification(mod)}]"

    return sys.intern(result)


@lru_cache(maxsize=512)
def serialize_peptidoform(pf: "Peptidoform") -> str:  # pragma: no cover
    """
    Serialize a Peptidoform to ProForma notation.

    Args:
        pf: Peptidoform to serialize

    Returns:
        String representation like '[Oxidation]?{Glycan:Hex}-PEPTIDE-[Amidated]'
    """
    parts: list[str] = []

    # Add name if present: (>name)
    if pf.name:
        parts.append(f"(>{pf.name})")

    # Add unlocalised modifications with ? mark
    if pf.unlocalised_modifications:
        for mod in pf.unlocalised_modifications:
            parts.append(f"[{serialize_modification(mod)}]")
        parts.append("?")

    # Add labile modifications in curly braces
    if pf.labile_modifications:
        for labile in pf.labile_modifications:
            parts.append(f"{{{serialize_modification_tags(labile)}}}")

    # Add N-terminal modifications with dash
    if pf.n_term_modifications:
        for mod in pf.n_term_modifications:
            parts.append(f"[{serialize_modification(mod)}]")
        parts.append("-")

    # Add sequence - handle different Sequence types
    seq_parts: list[str] = []
    for s in pf.sequence:
        if isinstance(s, tuple):
            # tuple[SequenceElement, ...] - iterate through elements
            seq_parts.extend(str(elem) for elem in s)
        else:
            # SequenceElement or SequenceRegion - use __str__
            seq_parts.append(str(s))
    parts.append("".join(seq_parts))

    # Add C-terminal modifications with dash
    if pf.c_term_modifications:
        parts.append("-")
        for mod in pf.c_term_modifications:
            parts.append(f"[{serialize_modification(mod)}]")

    return sys.intern("".join(parts))


@lru_cache(maxsize=512)
def serialize_peptidoform_ion(pfi: "PeptidoformIon") -> str:  # pragma: no cover
    """
    Serialize a PeptidoformIon to ProForma notation.

    Args:
        pfi: PeptidoformIon to serialize

    Returns:
        String representation like '(>>name)PEPTIDE//SEQUENCE/2' or 'PEPTIDE/[Na:z+1,H:z+1]'
    """
    parts: list[str] = []

    # Add name if present: (>>name)
    if pfi.name:
        parts.append(f"(>>{pfi.name})")

    # Add peptidoforms separated by //
    peptidoform_strs = [serialize_peptidoform(pf) for pf in pfi.peptidoforms]
    parts.append("//".join(peptidoform_strs))

    # Add charge if present
    if pfi.charge is not None:
        if isinstance(pfi.charge, int):
            # Simple integer charge: /2
            parts.append(f"/{pfi.charge}")
        else:
            # Charge carriers: /[Na:z+1,H:z+1]
            carrier_strs = [serialize_global_charge_carrier(cc) for cc in pfi.charge]
            parts.append(f"/[{','.join(carrier_strs)}]")

    return sys.intern("".join(parts))


@lru_cache(maxsize=256)
def serialize_compound_peptidoform_ion(cpfi: "CompoundPeptidoformIon") -> str:  # pragma: no cover
    """
    Serialize a CompoundPeptidoformIon to ProForma notation.

    Args:
        cpfi: CompoundPeptidoformIon to serialize

    Returns:
        String representation like '(>>>name)<13C><[Oxidation]@M>PEPTIDE+SEQUENCE'
    """
    parts: list[str] = []

    # Add name if present: (>>>name)
    if cpfi.name:
        parts.append(f"(>>>{cpfi.name})")

    # Add isotope replacements: <13C><15N>
    if cpfi.isotope_replacement:
        for iso in cpfi.isotope_replacement:
            parts.append(f"<{serialize_isotope_replacement(iso)}>")

    # Add fixed modifications: <[mod]@position>
    if cpfi.fixed_modifications:
        for fixed_mod in cpfi.fixed_modifications:
            parts.append(f"<{serialize_fixed_modification(fixed_mod)}>")

    # Add peptidoform ions separated by +
    ion_strs = [serialize_peptidoform_ion(pfi) for pfi in cpfi.peptidoform_ions]
    parts.append("+".join(ion_strs))

    return sys.intern("".join(parts))

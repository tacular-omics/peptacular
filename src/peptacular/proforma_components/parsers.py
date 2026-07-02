"""
Parser functions for converting strings to component objects.

This module contains all parsing logic, independent of other modules.
"""

import re
from functools import lru_cache
from typing import TYPE_CHECKING

from tacular import AminoAcid, Element, Monosaccharide

from ..constants import CV

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
        IsotopeReplacement,
        ModificationAmbiguousPrimary,
        ModificationAmbiguousSecondary,
        ModificationCrossLinker,
        ModificationTags,
        Peptidoform,
        PeptidoformIon,
        PositionRule,
        PositionScore,
        SequenceElement,
        SequenceRegion,
        TagAccession,
        TagCustom,
        TagInfo,
        TagMass,
        TagName,
    )


# Try to match known monosaccharide names (longest first)
monosaccharide_names: list[str] = sorted(list(Monosaccharide), key=len, reverse=True)  # type: ignore

# ============================================================================
# Regex Patterns for Modification Types
# ============================================================================

# Accession patterns: MUST use full CV name (UNIMOD, MOD, RESID, GNO, XLMOD)
# Rule 1: Full CV name with colon, MUST NOT be followed with +/-
# e.g., "UNIMOD:35", "MOD:00719", "RESID:1234", "GNO:5678", "XLMOD:91011"
PATTERN_ACCESSION: re.Pattern[str] = re.compile(r"^(UNIMOD|MOD|RESID|GNO|XLMOD):(\w+)$", re.IGNORECASE)

# Custom pattern: Rule 3: MUST use shorthand "C:"
# e.g., "C:some custom data"
PATTERN_CUSTOM: re.Pattern[str] = re.compile(r"^C:(.+)$", re.IGNORECASE)
# Mass pattern: Rule 2: MUST include sign (+/-), can optionally use CV single letter prefix
# Format: [single-letter-CV:]?[Obs:]?<sign><number>[.number]
# e.g., "+15.995", "-18.010", "U:+15.995", "M:-18.010", "Obs:+17.05685", "U:Obs:+15.995"
PATTERN_MASS: re.Pattern[str] = re.compile(r"^(?:([UMRXGumrxg]):)?(?:(Obs):)?([+-]\d+(?:\.\d+)?)$", re.IGNORECASE)

# Named modification pattern: Rule 4: Can optionally use shorthand CVs (including C: for custom), MUST NOT be followed with +/-
# Format: [single-letter-CV:]?<name>
# e.g., "Oxidation", "Phospho", "U:Oxidation", "M:Phospho", "R:Acetyl", "C:CustomMod"
PATTERN_NAMED_MOD: re.Pattern[str] = re.compile(r"^(?:([UMRXGCumrxgc]):)?(?![+-]|Obs:)(.+)$", re.IGNORECASE)

# Regex pattern for parsing formula elements
# Matches two formats:
# 1. [isotope Element count]  (e.g., [13C2], [13C-2], [2H4])
#    Note: When isotope is specified, count MUST be inside the bracket
# 2. Element count            (e.g., C2, H, Ca3, N-1)
_FORMULA_ELEMENT_PATTERN: re.Pattern[str] = re.compile(
    r"^"
    r"(?:\[(\d+)([A-Z][a-z]?)([+-]?\d+)?\])|"  # Format 1: [isotope Element count]
    r"(?:([A-Z][a-z]?)([+-]?\d+)?)"  # Format 2: Element count
    r"$"
)


# Regex pattern for parsing charged formulas
# Matches: Formula:C2H6 or Formula:C2H6:z+2
_CHARGED_FORMULA_PATTERN: re.Pattern[str] = re.compile(r"^Formula:(.+?)(?::z([+-]?\d+))?$", re.IGNORECASE)


def parse_formula_element(s: str, allow_zero: bool = False) -> "FormulaElement":
    """
    Parse a formula element string like '[13C2]' or 'H2' or 'C'

    Follows the ProForma 2.0 specification for elemental formulas:
    - Formula Rule 1: Pairs of atoms and cardinality (C2).
    - Formula Rule 2: Cardinalities are positive or negative integers. Zero not supported.
                      Default is +1 if not specified.
    - Formula Rule 3: Isotopes prefixed with isotopic number in brackets: [13C2]
                      IMPORTANT: When isotope is specified, the count MUST be inside the bracket

    Valid formats:
        - Element only: "C", "H", "Ca"
        - Element with count: "C2", "H12", "N-1"
        - Isotope with element and count: "[13C2]", "[2H4]", "[15N-1]"
        - Isotope with element (count defaults to 1): "[13C]", "[2H]"

    INVALID formats:
        - "[13C]2" - count must be inside bracket when isotope is specified

    Examples:
        "C" -> FormulaElement(C, 1, None)
        "C2" -> FormulaElement(C, 2, None)
        "H-2" -> FormulaElement(H, -2, None)
        "[13C]" -> FormulaElement(C, 1, 13)
        "[13C2]" -> FormulaElement(C, 2, 13)
        "[13C-1]" -> FormulaElement(C, -1, 13)

    Args:
        s: String representation of a single formula element

    Returns:
        FormulaElement object

    Raises:
        ValueError: If the string cannot be parsed or is invalid
    """
    from .comps import FormulaElement

    if not s:
        raise ValueError("Empty formula element string")

    s = s.strip()
    match = _FORMULA_ELEMENT_PATTERN.match(s)

    if not match:
        raise ValueError(f"Invalid formula element: '{s}'")

    # Extract groups: [isotope Element count] or Element count
    isotope_str, element_bracketed, count_bracketed, element_plain, count_plain = match.groups()

    # Determine which format was used
    if isotope_str and element_bracketed:
        # Format: [13C2] or [13C]
        element_symbol = element_bracketed
        isotope = int(isotope_str)
        count = int(count_bracketed) if count_bracketed else 1
    elif element_plain:
        # Format: C2 or C
        element_symbol = element_plain
        isotope = None
        count = int(count_plain) if count_plain else 1
    else:
        raise ValueError(f"No element found in '{s}'")

    # Validate count (zero not allowed)
    if count == 0 and not allow_zero:
        raise ValueError(f"Zero cardinality not allowed in '{s}'")

    # Validate element symbol
    try:
        if element_symbol == "D" or element_symbol == "T":
            # Special cases for Deuterium and Tritium
            element_symbol = "H"
        element = Element(element_symbol)
    except ValueError as e:
        raise ValueError(f"Unknown element symbol: '{element_symbol}'") from e

    return FormulaElement(element=element, occurance=count, isotope=isotope)


@lru_cache(maxsize=1024)
def parse_charged_formula(s: str, allow_zero: bool = False, require_formula_prefix: bool = True, sep: str = "") -> "ChargedFormula":
    """
    Parse a charged formula string like 'Formula:C2H6' or 'Formula:C2H6:z+2'

    Follows the ProForma 2.0 specification for formula modifications.

    Valid formats:
        - Formula only: "Formula:C2H6"
        - Formula with charge: "Formula:C2H6:z+2", "Formula:H2O:z-1"

    Examples:
        "Formula:C2H6" -> ChargedFormula((C, 2), (H, 6), None)
        "Formula:C2H6:z+2" -> ChargedFormula((C, 2), (H, 6), 2)
        "Formula:H-2O-1:z-1" -> ChargedFormula((H, -2), (O, -1), -1)

    Args:
        s: String representation of a charged formula

    Returns:
        ChargedFormula object

    Raises:
        ValueError: If the string cannot be parsed or is invalid
    """
    from .comps import ChargedFormula

    if sep != "":
        s = s.replace(sep, "")

    if not s:
        raise ValueError("Empty charged formula string")

    s = s.strip()

    match = LOC_PATTERN.match(s)
    if match:
        s, pos, score = match.groups()
        if score is not None:
            score = float(score)
    else:
        s, pos, score = s, None, None

    if require_formula_prefix:
        match = _CHARGED_FORMULA_PATTERN.match(s)

        if not match:
            raise ValueError(f"Invalid charged formula: '{s}'")
    else:
        # Try to parse without the Formula: prefix
        match = re.match(r"^(.+?)(?::z([+-]?\d+))?$", s, re.IGNORECASE)
        if not match:
            raise ValueError(f"Invalid charged formula: '{s}'")

    formula_str, charge_str = match.groups()

    # Parse the formula string into FormulaElements
    formula_elements = _parse_formula_string(formula_str, allow_zero=allow_zero)

    # Parse charge if present
    charge = int(charge_str) if charge_str else None

    return ChargedFormula(formula=formula_elements, charge=charge, position_id=pos, score=score)


def _parse_formula_string(formula_str: str, allow_zero: bool = False) -> tuple["FormulaElement", ...]:
    """
    Parse a chemical formula string into FormulaElements.

    This is a helper function used by parse_charged_formula.
    For the full specification, see parse_formula_element.

    Follows ProForma 2.0 notation:
    - Isotopes must have count inside bracket: [13C2] not [13C]2
    - Spaces allowed between element pairs: C2 H6 or C2H6

    Examples:
        "C2H6" -> (FormulaElement(C, 2), FormulaElement(H, 6))
        "[13C2]H6" -> (FormulaElement(C, 2, 13), FormulaElement(H, 6))
        "[13C2][12C-2]H2N" -> isotope replacement notation
    """
    from .comps import FormulaElement

    elements: list[FormulaElement] = []
    i = 0

    while i < len(formula_str):
        # Skip whitespace (spaces are allowed between element pairs)
        if formula_str[i].isspace():
            i += 1
            continue

        # Try to parse an element starting at position i
        # We need to find the end of this element specification
        start = i

        # Check if it starts with isotope bracket [13C2]
        if formula_str[i] == "[":
            # Find the closing bracket - everything must be inside
            close_bracket = formula_str.find("]", i)
            if close_bracket == -1:
                raise ValueError(f"Unclosed isotope bracket in '{formula_str}'")
            i = close_bracket + 1
            # Extract the entire bracketed element including the brackets
            element_str = formula_str[start:i]
        # Check if it's a regular element (starts with uppercase letter)
        elif formula_str[i].isupper():
            i += 1
            # Add lowercase letter if present (for elements like Ca, Br)
            if i < len(formula_str) and formula_str[i].islower():
                i += 1

            # Now parse optional count (including sign)
            if i < len(formula_str) and (formula_str[i].isdigit() or formula_str[i] in "+-"):
                if formula_str[i] in "+-":
                    i += 1
                while i < len(formula_str) and formula_str[i].isdigit():
                    i += 1

            # Extract the element string
            element_str = formula_str[start:i]
        else:
            raise ValueError(f"Unexpected character '{formula_str[i]}' at position {i} in '{formula_str}'")

        # Parse the element string
        try:
            element = parse_formula_element(element_str, allow_zero=allow_zero)
            elements.append(element)
        except ValueError as e:
            raise ValueError(f"Failed to parse element '{element_str}' in formula '{formula_str}'") from e

    if not elements:
        raise ValueError(f"No elements found in formula '{formula_str}'")

    return tuple(elements)


def parse_position_score(s: str) -> "PositionScore":
    from .comps import PositionScore

    s = s.strip()

    # example : #1
    # example: #4(0.95)
    pos = None
    score = None
    match = LOC_PATTERN.match("X" + s)
    if match:
        s, pos, score_str = match.groups()
        if score_str is not None:
            score = float(score_str)
    if pos is None:
        raise ValueError(f"Invalid position id in position score string: '{s}'")

    if len(s) > 1:
        raise ValueError(f"Invalid position score string: '{s[1:]}'")

    return PositionScore(position_id=pos, score=score)


# ============================================================================
# Modification Tag Parsing
# ============================================================================


@lru_cache(maxsize=1024)
def parse_modification_tag(mod_str: str) -> "MODIFICATION_TAG_TYPE":
    """
    Parse a ProForma modification string into its corresponding tag object.

    Rules:
    1. Accession Tags: MUST use full CV name (UNIMOD, MOD, RESID, GNO, XLMOD) with colon
       MUST NOT be followed with +/-
       e.g., "UNIMOD:35", "MOD:00719", "RESID:1234", "GNO:5678", "XLMOD:91011"

    2. Mass Delta Tags: MUST include sign (+/-), can optionally use CV single letter prefix
       e.g., "+15.995", "-18.010", "U:+15.995", "M:-18.010"

    3. Custom Tags: MUST use shorthand "C:"
       e.g., "C:some custom data"

    4. Named Modifications: Can optionally use shorthand CVs (single letter), MUST NOT be followed with +/-
       e.g., "Oxidation", "Phospho", "U:Oxidation", "M:Phospho", "R:Acetyl"

    Also handles:
    - Formulas: "Formula:C2H6", "Formula:C2H6:z+2"
    - Glycans: "Glycan:Hex5HexNAc4" (prefix mandatory)
    - INFO tags: "INFO:custom text"

    Args:
        mod_str: The modification string to parse

    Returns:
        A ModificationTag object (TagAccession, TagMass, TagName, GlycanTag, etc.)

    Raises:
        ValueError: If the modification string cannot be parsed
    """
    from .comps import ComkpTag, ComupTag, GlycanTag, LimitTag, PositionTag, TagName

    mod_str = mod_str.strip()

    if not mod_str:
        raise ValueError("Empty modification string")

    mod_str_lower = mod_str.lower()

    if mod_str_lower.startswith("#"):
        return parse_position_score(mod_str)

    # Try to match accession patterns (Rule 1: full CV names only)
    if ":" in mod_str_lower:
        accession_term = mod_str_lower.split(":", 1)[0].upper()
        if accession_term in {"UNIMOD", "MOD", "RESID", "GNO", "XLMOD"}:
            return parse_tag_accession(mod_str)

    # Try to match mass pattern (Rule 2: must have sign, optional single letter CV prefix)
    if mod_str_lower.startswith(("+", "-")) or ":" in mod_str_lower and mod_str_lower.split(":", 1)[1].startswith(("+", "-")):
        return parse_tag_mass(mod_str)

    # Try to match formula pattern
    if mod_str_lower.startswith("formula:"):
        return parse_charged_formula(mod_str)

    # Try to match glycan pattern
    if mod_str_lower.startswith("glycan:"):
        match = LOC_PATTERN.match(mod_str)
        if match:
            mod_str, pos, score = match.groups()
            if score is not None:
                score = float(score)
        else:
            mod_str, pos, score = mod_str, None, None

        glycan_components = parse_glycan(mod_str)
        return GlycanTag(components=glycan_components, position_id=pos, score=score)

    # Try to match INFO pattern
    if mod_str_lower.startswith("info:"):
        return parse_tag_info(mod_str)

    if mod_str_lower.startswith("position:"):
        return PositionTag.from_string(mod_str)

    if mod_str_lower.startswith("limit:"):
        return LimitTag.from_string(mod_str)

    if mod_str_lower == "comkp":
        return ComkpTag.from_string(mod_str)

    if mod_str_lower == "comup":
        return ComupTag.from_string(mod_str)

    # Default: treat as named modification (Rule 4: optional single letter CV prefix, no +/-)
    try:
        return parse_tag_name(mod_str)
    except ValueError:
        return TagName(name=mod_str, cv=None)


# Regex pattern for parsing glycan compositions
# Matches: Glycan:Hex5HexNAc4
_GLYCAN_PATTERN = re.compile(r"^Glycan:(.+)$", re.IGNORECASE)


# TODO: Support mixed glycan and formulas?  Proforma seems to support this yet Peptacular only allows glycans.
@lru_cache(maxsize=1024)
def parse_glycan(s: str) -> tuple["GlycanComponent", ...]:
    """
    Parse a glycan string like 'Glycan:Hex5HexNAc4'

    Follows the ProForma 2.0 specification for glycan modifications.
    The Glycan: prefix is mandatory.

    Valid formats:
        - "Glycan:Hex5HexNAc4"
        - "Glycan:Hex5HexNAc4NeuAc2"

    Examples:
        "Glycan:Hex5HexNAc4" -> (GlycanComponent(Hex, 5), GlycanComponent(HexNAc, 4))
        "Glycan:Hex5HexNAc4NeuAc2" -> (GlycanComponent(Hex, 5), GlycanComponent(HexNAc, 4), GlycanComponent(NeuAc, 2))

    Args:
        s: String representation of a glycan composition

    Returns:
        Tuple of GlycanComponent objects

    Raises:
        ValueError: If the string cannot be parsed or is invalid
    """
    if not s:
        raise ValueError("Empty glycan string")

    s = s.strip()
    match = _GLYCAN_PATTERN.match(s)

    if not match:
        raise ValueError(f"Invalid glycan string (must start with 'Glycan:'): '{s}'")

    glycan_str = match.group(1)

    # Parse the glycan composition into GlycanComponents
    return _parse_glycan_composition(glycan_str)


def _parse_glycan_composition(glycan_str: str) -> tuple["GlycanComponent", ...]:
    """
    Parse a glycan composition string into GlycanComponents.

    This is a helper function used by parse_glycan.

    Examples:
        "Hex5HexNAc4" -> (GlycanComponent(Hex, 5), GlycanComponent(HexNAc, 4))
        "Hex5HexNAc4NeuAc2" -> (GlycanComponent(Hex, 5), GlycanComponent(HexNAc, 4), GlycanComponent(NeuAc, 2))
    """
    from .comps import GlycanComponent

    components: list[GlycanComponent] = []
    i = 0

    while i < len(glycan_str):
        # Try to match a monosaccharide name
        matched = False
        for mono_name in monosaccharide_names:
            if glycan_str[i:].startswith(mono_name):
                # Found a match
                i += len(mono_name)

                # Parse count
                count_start = i
                while i < len(glycan_str) and glycan_str[i].isdigit():
                    i += 1

                if i > count_start:
                    count = int(glycan_str[count_start:i])
                else:
                    count = 1

                try:
                    monosaccharide = Monosaccharide(mono_name)
                except ValueError as e:
                    raise ValueError(f"Unknown monosaccharide: {mono_name}") from e

                components.append(GlycanComponent(monosaccharide=monosaccharide, occurance=count))
                matched = True
                break

        if not matched:
            raise ValueError(f"Could not parse glycan composition at position {i}: '{glycan_str[i:]}'")

    return tuple(components)


@lru_cache(maxsize=1024)
def parse_modification_tags(mod_str: str) -> "ModificationTags":
    """
    Parse a modification string that may contain multiple tags separated by pipe (|).

    ProForma allows multiple representations of the same modification to be specified
    together, separated by the pipe character. For example:
    - [Oxidation|INFO:probable]
    - [UNIMOD:35|+15.995|Formula:O]

    Args:
        mod_str: The modification string to parse (may contain | separators)

    Returns:
        ModificationTags object containing the parsed tags
    """
    from .comps import ModificationTags

    # Split on pipe character to handle multiple tag definitions
    if "|" in mod_str:
        parts = mod_str.split("|")
        tags: list[MODIFICATION_TAG_TYPE] = []
        for part in parts:
            part = part.strip()
            tags.append(parse_modification_tag(part))
        return ModificationTags(tags=tuple(tags))
    else:
        # Single tag
        tag = parse_modification_tag(mod_str)
        return ModificationTags(tags=(tag,))


# Regex pattern for parsing position rules
# Matches: N-term, C-term, N-term:K, C-term:K, or just K
_POSITION_RULE_PATTERN = re.compile(r"^(?:(N-term|C-term|Protein N-term|Protein C-term):)?([A-Z])?$", re.IGNORECASE)


@lru_cache(maxsize=1024)
def parse_position_rule(s: str) -> "PositionRule":
    """
    Parse a position rule string like 'N-term', 'C-term:K', or 'M'

    Follows the ProForma 2.0 specification for position rules.

    Valid formats:
        - Terminal only: "N-term", "C-term", "Protein N-term", "Protein C-term"
        - Terminal with amino acid: "N-term:K", "C-term:K"
        - Amino acid only: "K", "M" (implies ANYWHERE)

    Examples:
        "N-term" -> PositionRule(N-term, None)
        "C-term:K" -> PositionRule(C-term, K)
        "M" -> PositionRule(ANYWHERE, M)

    Args:
        s: String representation of a position rule

    Returns:
        PositionRule object

    Raises:
        ValueError: If the string cannot be parsed or is invalid
    """
    from ..constants import Terminal
    from .comps import PositionRule

    if not s:
        raise ValueError("Empty position rule string")

    s = s.strip()

    # Check if it contains a colon (terminal:amino_acid)
    if ":" in s:
        terminal_str, aa_str = s.split(":", 1)
        terminal_str = terminal_str.strip()
        aa_str = aa_str.strip()

        try:
            terminal = Terminal.from_str(terminal_str)
        except ValueError as e:
            raise ValueError(f"Unknown terminal: '{terminal_str}'") from e

        if aa_str:
            try:
                amino_acid = AminoAcid(aa_str)
            except ValueError as e:
                raise ValueError(f"Unknown amino acid: '{aa_str}'") from e
            return PositionRule(terminal=terminal, amino_acid=amino_acid)
        else:
            return PositionRule(terminal=terminal, amino_acid=None)
    else:
        # No colon - check if it's a terminal or amino acid
        try:
            terminal = Terminal.from_str(s)
            return PositionRule(terminal=terminal, amino_acid=None)
        except ValueError:
            # Must be an amino acid (ANYWHERE position)
            try:
                amino_acid = AminoAcid.from_str(s)
                return PositionRule(terminal=Terminal.ANYWHERE, amino_acid=amino_acid)
            except ValueError as e:
                raise ValueError(f"Unknown terminal or amino acid: '{s}'") from e


LOC_PATTERN = re.compile(r"^([^#(]+)(?:#([^(]+))?(?:\(([^)]+)\))?$")


@lru_cache(maxsize=512)
def parse_tag_accession(s: str) -> "TagAccession":
    """
    Parse an accession string using FULL CV names only.

    Rule 1: MUST use full CV name (UNIMOD, MOD, RESID, GNO, XLMOD) with colon
    MUST NOT be followed with +/-

    Supported formats:
    - UNIMOD:35
    - MOD:00719 (PSI-MOD)
    - RESID:AA0038
    - GNO:G12345
    - XLMOD:02001

    Args:
        s: String representation of an accession tag

    Returns:
        TagAccession object

    Raises:
        ValueError: If the string is not a valid accession format
    """

    from .comps import TagAccession

    s = s.strip()

    cv, accession = s.split(":", 1)
    match = LOC_PATTERN.match(accession)
    if match:
        name, pos, score = match.groups()
        if score is not None:
            score = float(score)
    else:
        name, pos, score = accession, None, None

    match cv_upper := cv.upper():
        case "UNIMOD":
            return TagAccession(accession=name, cv=CV.UNIMOD, position_id=pos, score=score)
        case "MOD":
            return TagAccession(accession=name, cv=CV.PSI_MOD, position_id=pos, score=score)
        case "RESID":
            return TagAccession(accession=name, cv=CV.RESID, position_id=pos, score=score)
        case "GNO":
            return TagAccession(accession=name, cv=CV.GNOME, position_id=pos, score=score)
        case "XLMOD":
            return TagAccession(accession=name, cv=CV.XL_MOD, position_id=pos, score=score)
        case _:
            raise ValueError(f"Unknown CV: {cv_upper} (must be UNIMOD, MOD, RESID, GNO, or XLMOD)")

    return TagAccession(accession=accession, cv=cv)


@lru_cache(maxsize=512)
def parse_tag_mass(s: str) -> "TagMass":
    """
    Parse a mass delta string according to ProForma 2.1 specification.

    Rule 2: MUST include sign (+/-), can optionally use CV single letter prefix

    Valid formats:
        - Basic: "+15.995", "-18.010", "+79.9663"
        - Observed: "Obs:+15.995", "Obs:-18.01"
        - CV-tagged (single letter): "U:+15.995", "M:-18.01", "R:+10.0", "X:+138.0", "G:+12.0", "C:+100.0"
        - Combined: "U:Obs:+15.995"

    Single letter CV prefixes:
        - U: UNIMOD
        - M: PSI-MOD (MOD)
        - R: RESID
        - X: XLMOD
        - G: GNO
        - C: Custom (CV.CUSTOM)

    Examples:
        "+15.995" -> TagMass(mass=15.995, cv=None)
        "-18.010" -> TagMass(mass=-18.010, cv=None)
        "Obs:+79.9663" -> TagMass(mass=79.9663, cv=CV.OBSERVED)
        "U:+15.995" -> TagMass(mass=15.995, cv=CV.UNIMOD)
        "M:-18.01" -> TagMass(mass=-18.01, cv=CV.PSI_MOD)
        "C:+100.0" -> TagMass(mass=100.0, cv=CV.CUSTOM)
        "U:Obs:+15.995" -> TagMass(mass=15.995, cv=CV.UNIMOD)

    Args:
        s: String representation of a mass delta

    Returns:
        TagMass object

    Raises:
        ValueError: If the string cannot be parsed or lacks required sign
    """
    from .comps import TagMass

    match = LOC_PATTERN.match(s)
    if match:
        s, pos, score = match.groups()
        if score is not None:
            score = float(score)
    else:
        s, pos, score = s, None, None

    if s.startswith(("+", "-")):
        # No CV prefix, just mass
        return TagMass(mass_str=s, cv=None, position_id=pos, score=score)

    cv, mass_str = s.split(":", 1)
    match cv.upper():
        case "U":
            return TagMass(mass_str=mass_str, cv=CV.UNIMOD, position_id=pos, score=score)
        case "M":
            return TagMass(mass_str=mass_str, cv=CV.PSI_MOD, position_id=pos, score=score)
        case "R":
            return TagMass(mass_str=mass_str, cv=CV.RESID, position_id=pos, score=score)
        case "X":
            return TagMass(mass_str=mass_str, cv=CV.XL_MOD, position_id=pos, score=score)
        case "G":
            return TagMass(mass_str=mass_str, cv=CV.GNOME, position_id=pos, score=score)
        case "C":
            return TagMass(mass_str=mass_str, cv=CV.CUSTOM, position_id=pos, score=score)
        case "OBS":
            return TagMass(mass_str=mass_str, cv=CV.OBSERVED, position_id=pos, score=score)
        case _:
            raise ValueError(f"Invalid CV prefix: {cv} (must be U, M, R, X, G, or C)")


@lru_cache(maxsize=512)
def parse_tag_name(s: str) -> "TagName | TagCustom":
    """
    Parse a named modification string.

    Rule 4: Can optionally use shorthand CVs (single letter, including C: for custom), MUST NOT be followed with +/-

    Valid formats:
        - "Oxidation", "Phospho", "Acetyl"
        - "U:Oxidation", "M:Phospho", "R:Acetyl", "C:CustomMod"

    Single letter CV prefixes:
        - U: UNIMOD (Optional)
        - M: PSI-MOD (MOD) (Optional)
        - R: RESID
        - X: XLMOD
        - G: GNO
        - C: Custom (CV.CUSTOM)

    Args:
        s: String representation of a tag name

    Returns:
        TagName object

    Raises:
        ValueError: If the string contains +/- or is invalid
    """
    from .comps import TagCustom, TagName

    s = s.strip()

    match = LOC_PATTERN.match(s)
    if match:
        s, pos, score = match.groups()
        if score is not None:
            score = float(score)
    else:
        s, pos, score = s, None, None

    if ":" not in s:
        # No CV prefix, just name
        return TagName(name=s, cv=None, position_id=pos, score=score)

    cv, name = s.split(":", 1)
    match cv.upper():
        case "U":
            return TagName(name=name, cv=CV.UNIMOD, position_id=pos, score=score)
        case "M":
            return TagName(name=name, cv=CV.PSI_MOD, position_id=pos, score=score)
        case "R":
            return TagName(name=name, cv=CV.RESID, position_id=pos, score=score)
        case "X":
            return TagName(name=name, cv=CV.XL_MOD, position_id=pos, score=score)
        case "G":
            return TagName(name=name, cv=CV.GNOME, position_id=pos, score=score)
        case "C":
            return TagCustom(name=name, position_id=pos, score=score)
        case _:
            raise ValueError(f"Invalid CV prefix: {cv} (must be U, M, R, X, G, or C)")


@lru_cache(maxsize=512)
def parse_tag_info(s: str) -> "TagInfo":
    """
    Parse an INFO string like 'INFO:some text'

    Args:
        s: String representation of an INFO tag

    Returns:
        TagInfo object

    Raises:
        ValueError: If the string is not a valid INFO tag
    """
    from .comps import TagInfo

    if not s.lower().startswith("info:"):
        raise ValueError(f"Invalid INFO string: {s}")
    info = s[5:]  # Remove 'INFO:' prefix
    return TagInfo(info=info)


@lru_cache(maxsize=512)
def parse_tag_custom(s: str) -> "TagCustom":
    """
    Parse a custom tag string

    Args:
        s: String representation of a custom tag

    Returns:
        TagCustom object

    Raises:
        ValueError: If the string doesn't have the required C: prefix
    """
    from .comps import TagCustom

    # Custom tags must be in format "C:name"
    if not s.startswith("C:"):
        raise ValueError(f"Custom tag must start with 'C:': {s}")

    name = s[2:]  # Remove 'C:' prefix
    return TagCustom(name=name)


@lru_cache(maxsize=512)
def parse_glycan_component(s: str) -> "GlycanComponent":
    """
    Parse a glycan component string like 'Hex5' or 'HexNAc4'

    Args:
        s: String representation of a glycan component

    Returns:
        GlycanComponent object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import GlycanComponent

    s = s.strip()

    for mono_name in monosaccharide_names:
        if s.startswith(mono_name):
            # Found a match
            count_str = s[len(mono_name) :]
            if count_str:
                try:
                    count = int(count_str)
                except ValueError:
                    continue  # Not a valid count, try next monosaccharide
            else:
                count = 1

            try:
                monosaccharide = Monosaccharide(mono_name)
            except ValueError:
                continue

            return GlycanComponent(monosaccharide=monosaccharide, occurance=count)

    # If no monosaccharide match, it might be a formula in parentheses
    if s.startswith("(") and ")" in s:
        close_paren = s.find(")")
        formula_str = s[1:close_paren]
        count_str = s[close_paren + 1 :]

        # Parse the formula
        formula = parse_charged_formula(f"Formula:{formula_str}")

        # Parse count
        count = int(count_str) if count_str else 1

        return GlycanComponent(monosaccharide=formula, occurance=count)

    raise ValueError(f"Could not parse glycan component: '{s}'")


@lru_cache(maxsize=512)
def parse_isotope_replacement(s: str) -> "IsotopeReplacement":
    """
    Parse an isotope replacement string like '13C' or '15N' or 'D'

    Args:
        s: String representation of an isotope replacement

    Returns:
        IsotopeReplacement object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import IsotopeReplacement

    s = s.strip()

    # Special case for Deuterium
    if s == "D":
        return IsotopeReplacement(element=Element.H, isotope=2)

    # Parse isotope number and element symbol
    # Format: <number><element>
    i = 0
    while i < len(s) and s[i].isdigit():
        i += 1

    if i == 0:
        raise ValueError(f"Expected isotope number at start of string: {s}")

    isotope = int(s[:i])
    element_str = s[i:]

    if not element_str:
        raise ValueError(f"Missing element symbol in: {s}")

    try:
        element = Element(element_str)
    except ValueError as e:
        raise ValueError(f"Unknown element symbol: {element_str}") from e

    return IsotopeReplacement(element=element, isotope=isotope)


@lru_cache(maxsize=512)
def parse_global_charge_carrier(s: str) -> "GlobalChargeCarrier":
    """
    Parse a charge carrier string like 'Na:z+1' or 'H:z+1^2'.

    A charge carrier is a *bare* charged formula (``<formula>:z<charge>``) per ProForma
    2.1 section 11.5, optionally with an occurrence specifier ``^n``. A CV/type prefix such
    as ``Formula:`` or ``Glycan:`` is only valid for a localised residue modification, not
    for a charge carrier, and is rejected here.

    Args:
        s: String representation of a global charge carrier

    Returns:
        GlobalChargeCarrier object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import GlobalChargeCarrier

    s = s.strip()

    # Check for occurrence specification (^number at the end)
    occurance = 1
    if "^" in s:
        formula_part, occ_str = s.rsplit("^", 1)
        occurance = int(occ_str)
    else:
        formula_part = s

    # A charge carrier must be a bare charged formula (e.g. 'Na:z+1', 'C2H6:z+2',
    # '[15N1]H4:z+1'). The formula portion (before the ':z<charge>' suffix) never contains a
    # ':' — a colon there means an illegal 'Formula:'/'Glycan:'-style prefix was used.
    formula_only = formula_part.split(":z", 1)[0]
    if ":" in formula_only:
        prefix = formula_only.split(":", 1)[0]
        raise ValueError(
            f"Invalid charge carrier '{s}': a charge carrier must be a bare charged formula "
            f"like 'Na:z+1' or 'C2H6:z+2' (ProForma 2.1 section 11.5), not a '{prefix}:'-prefixed value."
        )

    # Parse the formula part
    charged_formula = parse_charged_formula(formula_part, require_formula_prefix=False)

    return GlobalChargeCarrier(charged_formula=charged_formula, occurance=occurance)


@lru_cache(maxsize=512)
def parse_modification_ambiguous_primary(s: str) -> "ModificationAmbiguousPrimary":
    """
    Parse an ambiguous primary modification string.
    Format: [Modification#label(score)|Position:...|Limit:...|CoMKP|CoMUP]

    Args:
        s: String representation of an ambiguous primary modification

    Returns:
        ModificationAmbiguousPrimary object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import ModificationAmbiguousPrimary

    s = s.strip()

    # Parse the different parts separated by |
    parts = s.split("|")

    # Find the part with the label (contains #)
    label_part = None
    score = None
    mod_tags_parts: list[str] = []
    position = None
    limit = None
    comkp = False
    comup = False

    for part in parts:
        if "#" in part:
            # This is the label part
            mod_part, label_part = part.split("#", 1)
            mod_tags_parts.append(mod_part)

        elif part.startswith("Position:"):
            # Position rules
            pos_str = part[9:]  # Remove "Position:" prefix
            pos_strs = pos_str.split(",")
            position = tuple(parse_position_rule(p.strip()) for p in pos_strs)
        elif part.startswith("Limit:"):
            # Limit
            limit = int(part[6:])
        elif part == "CoMKP":
            comkp = True
        elif part == "CoMUP":
            comup = True
        else:
            # This is a modification tag part
            mod_tags_parts.append(part)

    if label_part is None:
        raise ValueError(f"No label found in ambiguous primary modification: {s}")

    # Parse the label and optional score
    if "(" in label_part:
        label_score, score_str = label_part.split("(", 1)
        score = float(score_str.rstrip(")"))
        label = label_score.lstrip("#")
    else:
        label = label_part.lstrip("#")

    # Parse modification tags
    if mod_tags_parts:
        mod_tags_str = "|".join(mod_tags_parts)
        modification_tags = parse_modification_tags(mod_tags_str)
    else:
        from .comps import ModificationTags

        modification_tags = ModificationTags(tags=())

    return ModificationAmbiguousPrimary(
        label=label,
        tags=modification_tags,
        score=score,
        position=position,
        limit=limit,
        comkp=comkp,
        comup=comup,
    )


@lru_cache(maxsize=512)
def parse_modification_ambiguous_secondary(s: str) -> "ModificationAmbiguousSecondary":
    """
    Parse an ambiguous secondary modification string.
    Format: [#label(score)]

    Args:
        s: String representation of an ambiguous secondary modification

    Returns:
        ModificationAmbiguousSecondary object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import ModificationAmbiguousSecondary

    s = s.strip()

    if not s.startswith("#"):
        raise ValueError(f"Ambiguous secondary modification must start with #: {s}")

    # Parse label and optional score
    if "(" in s:
        label_part, score_str = s.split("(", 1)
        score = float(score_str.rstrip(")"))
        label = label_part.lstrip("#")
    else:
        label = s.lstrip("#")
        score = None

    return ModificationAmbiguousSecondary(label=label, score=score)


@lru_cache(maxsize=512)
def parse_modification_cross_linker(s: str) -> "ModificationCrossLinker":
    """
    Parse a cross-linker modification string.
    Three cases:
    1. Primary definition with tags and label: [XLMOD:02001#XL1]
    2. Secondary reference without tags: [#XL1]
    3. Branch: [MOD:00093#BRANCH] or [#BRANCH]

    Args:
        s: String representation of a cross-linker modification

    Returns:
        ModificationCrossLinker object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import ModificationCrossLinker

    s = s.strip()

    if "#" not in s:
        raise ValueError(f"Cross-linker modification must contain #: {s}")

    # Split on # to get tags and label
    parts = s.split("#", 1)

    if len(parts) == 1:
        # Just a label (shouldn't happen based on check above)
        raise ValueError(f"Invalid cross-linker format: {s}")

    tags_str = parts[0]
    label_str = parts[1]

    # Extract label (remove XL prefix if present)
    if label_str.startswith("XL"):
        label = label_str[2:]
    else:
        label = label_str

    # Parse tags if present
    if tags_str:
        modification_tags = parse_modification_tags(tags_str)
    else:
        modification_tags = None

    return ModificationCrossLinker(label=label, tags=modification_tags)


@lru_cache(maxsize=256)
def parse_fixed_modification(s: str) -> "FixedModification":
    """
    Parse a fixed modification string.
    Format: [mod]@position,position
    Example: [Oxidation]@M or [TMT6plex]@K,N-term

    Args:
        s: String representation of a fixed modification

    Returns:
        FixedModification object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import FixedModification

    s = s.strip()

    # Check if there are position rules
    if "@" in s:
        # Split on @ to get modification and position rules
        if not s.startswith("[") or "]@" not in s:
            raise ValueError(f"Invalid fixed modification format: {s}")

        # Find the ]@ separator
        bracket_end = s.index("]@")
        mod_str = s[1:bracket_end]  # Remove [ and ]
        pos_str = s[bracket_end + 2 :]  # Everything after ]@

        # Parse modification tags
        modification_tags = parse_modification_tags(mod_str)

        # Parse position rules
        pos_strs = pos_str.split(",")
        position_rules = tuple(parse_position_rule(p.strip()) for p in pos_strs)
    else:
        # No position rules
        if not s.startswith("[") or not s.endswith("]"):
            raise ValueError(f"Invalid fixed modification format: {s}")

        mod_str = s[1:-1]  # Remove [ and ]
        modification_tags = parse_modification_tags(mod_str)
        position_rules = ()

    return FixedModification(modifications=modification_tags, position_rules=position_rules)


def parse_modification(s: str) -> "MODIFICATION_TYPE":
    """
    Parse a modification string into its corresponding modification object.

    Determines the type of modification (standard tags, ambiguous primary/secondary, cross-linker)
    and delegates to the appropriate parsing function.

    Args:
        s: The modification string to parse

    Returns:
        A modification object (ModificationTags, ModificationAmbiguousPrimary, etc.)

    Raises:
        ValueError: If the modification string cannot be parsed
    """
    s = s.strip()

    if not s:
        raise ValueError("Empty modification string")

    # Check for ambiguous or cross-linker modifications
    if "#" in s:
        # The label follows the '#'; strip any (score) suffix and pipe-delimited tag
        # alternatives so we can inspect the bare label. ProForma reserves the ``XL``
        # prefix (e.g. #XL1) and the special ``BRANCH`` label for cross-links; every
        # other label denotes an ambiguous modification group.
        label = s.split("#", 1)[1].split("|")[0].split("(")[0].strip()
        is_crosslink = label.upper().startswith("XL") or label.upper() == "BRANCH"

        if is_crosslink:
            # Cross-linker definition ([XLMOD:02001#XL1]) or reference ([#XL1], [#BRANCH])
            return parse_modification_cross_linker(s)
        elif s.startswith("#"):
            # Ambiguous secondary reference ([#g1])
            return parse_modification_ambiguous_secondary(s)
        else:
            # Ambiguous primary definition ([Phospho#g1])
            return parse_modification_ambiguous_primary(s)
    else:
        # Standard modification tags
        return parse_modification_tags(s)


def _extract_bracketed_modifications(s: str) -> tuple["MODIFICATION_TYPE", ...]:
    """
    Helper function to extract modifications from a string with bracketed modifications.

    Parses strings like '[Oxidation][Phospho]' into a tuple of modification objects.

    Args:
        s: String containing bracketed modifications

    Returns:
        Tuple of modification objects

    Raises:
        ValueError: If brackets are unmatched or invalid characters are found
    """
    from .comps import MODIFICATION_TYPE

    if not s:
        return ()

    modifications: list[MODIFICATION_TYPE] = []
    i = 0

    while i < len(s):
        if s[i] == "[":
            # Find matching closing bracket
            depth = 1
            j = i + 1
            while j < len(s) and depth > 0:
                if s[j] == "[":
                    depth += 1
                elif s[j] == "]":
                    depth -= 1
                j += 1

            if depth != 0:
                raise ValueError(f"Unmatched bracket in modifications: {s}")

            # Extract modification string (without brackets)
            mod_str = s[i + 1 : j - 1]

            # Use the new parse_modification function to determine the type
            modifications.append(parse_modification(mod_str))
            i = j
        else:
            raise ValueError(f"Unexpected character '{s[i]}' at position {i} in modifications")

    return tuple(modifications)


@lru_cache(maxsize=1024)
def parse_sequence_element(s: str) -> "SequenceElement":
    """
    Parse a sequence element string like 'M[Oxidation]' or 'K'

    Args:
        s: String representation of a sequence element

    Returns:
        SequenceElement object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import SequenceElement

    if not s:
        raise ValueError("Empty sequence element string")

    # First character should be the amino acid
    aa = AminoAcid(s[0])

    if len(s) == 1:
        # No modifications
        return SequenceElement(amino_acid=aa, modifications=())

    # Parse modifications (everything after the amino acid)
    mods_str = s[1:]
    modifications = _extract_bracketed_modifications(mods_str)

    return SequenceElement(amino_acid=aa, modifications=modifications)


@lru_cache(maxsize=512)
def parse_sequence_region(s: str) -> "SequenceRegion":
    """
    Parse a sequence region string like '(PEPTIDE)[Oxidation]' or '(PEPTIDE)'

    Format: (SEQUENCE)[mod][mod]...
    - Sequence is enclosed in parentheses
    - Modifications are optional and appear after the closing parenthesis

    Examples:
        "(PEPTIDE)" -> SequenceRegion with sequence, no modifications
        "(PEPTIDE)[Oxidation]" -> SequenceRegion with sequence and one modification
        "(PEPTIDE)[Oxidation][Phospho]" -> SequenceRegion with sequence and two modifications

    Args:
        s: String representation of a sequence region

    Returns:
        SequenceRegion object

    Raises:
        ValueError: If the string cannot be parsed
    """
    from .comps import SequenceRegion

    if not s:
        raise ValueError("Empty sequence region string")

    s = s.strip()

    # Check if it starts with opening parenthesis
    if not s.startswith("("):
        raise ValueError(f"Sequence region must start with '(': {s}")

    ambiguous = s[1] == "?"

    # Find the closing parenthesis
    close_paren = s.find(")")
    if close_paren == -1:
        raise ValueError(f"Sequence region missing closing ')': {s}")

    # Extract the sequence part (between parentheses)
    seq_str: str = s[1 + ambiguous : close_paren]

    # Parse the sequence into SequenceElements
    # We need to split the sequence string into individual amino acid + modification pairs
    sequence_elements: list[SequenceElement] = []
    i = 0
    while i < len(seq_str):
        # Each element starts with an amino acid letter
        if not seq_str[i].isupper():
            raise ValueError(f"Expected amino acid at position {i} in sequence '{seq_str}'")

        # Find the extent of this sequence element (amino acid + any modifications)
        start = i
        i += 1  # Move past the amino acid

        # Skip any modifications (enclosed in square brackets)
        while i < len(seq_str) and seq_str[i] == "[":
            # Find matching closing bracket
            depth = 1
            j = i + 1
            while j < len(seq_str) and depth > 0:
                if seq_str[j] == "[":
                    depth += 1
                elif seq_str[j] == "]":
                    depth -= 1
                j += 1

            if depth != 0:
                raise ValueError(f"Unmatched bracket in sequence: {seq_str}")

            i = j

        # Extract this sequence element string and parse it using helper
        elem_str = seq_str[start:i]
        sequence_elements.append(parse_sequence_element(elem_str))

    # Parse modifications after the closing parenthesis using helper
    mods_str = s[close_paren + 1 :]
    modifications = _extract_bracketed_modifications(mods_str)

    return SequenceRegion(sequence=tuple(sequence_elements), modifications=tuple(modifications), ambiguous=ambiguous)


@lru_cache(maxsize=512)
def parse_peptidoform(s: str) -> "Peptidoform":
    """
    Parse a ProForma peptidoform string.

    Args:
        s: String representation of a peptidoform

    Returns:
        Peptidoform object

    Raises:
        NotImplementedError: Parser not yet implemented
    """
    raise NotImplementedError("Peptidoform parsing not yet implemented")


@lru_cache(maxsize=512)
def parse_peptidoform_ion(s: str) -> "PeptidoformIon":
    """
    Parse a ProForma peptidoform ion string.

    Args:
        s: String representation of a peptidoform ion

    Returns:
        PeptidoformIon object

    Raises:
        NotImplementedError: Parser not yet implemented
    """
    raise NotImplementedError("PeptidoformIon parsing not yet implemented")


@lru_cache(maxsize=256)
def parse_compound_peptidoform_ion(s: str) -> "CompoundPeptidoformIon":
    """
    Parse a ProForma compound peptidoform ion string.

    Args:
        s: String representation of a compound peptidoform ion

    Returns:
        CompoundPeptidoformIon object

    Raises:
        NotImplementedError: Parser not yet implemented
    """
    raise NotImplementedError("CompoundPeptidoformIon parsing not yet implemented")

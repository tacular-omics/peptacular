import sys
from collections.abc import Generator
from typing import NoReturn

from .mod import VALID_AMINO_ACIDS, Interval

_VALID_AA_SET = frozenset(VALID_AMINO_ACIDS)
_TERMINATOR_SET = frozenset(["/", "+"])
_DIGIT_SET = frozenset("0123456789")
_SIGN_SET = frozenset("+-")


class ProFormaParser:
    """
    A recursive descent parser for ProForma 2.1.

    Implements a state-machine approach to handle the hierarchy:
    1. Compound Ion (Chimeras separated by '+')
    2. Peptidoform Ion (Crosslinks separated by '//' sharing a charge)
    3. Peptidoform (The specific sequence and local mods)
    """

    __slots__ = (
        "original_sequence",
        "cursor",
        "length",
        "amino_acids",
        "compound_name",
        "ion_name",
        "peptide_name",
        "global_mods",
        "labile_mods",
        "nterm_mods",
        "cterm_mods",
        "internal_mods",
        "unknown_mods",
        "intervals",
        "charge",
        "charge_adducts",
        "is_chimeric",
        "is_crosslinked",
    )

    def __init__(self, proforma_sequence: str):
        self.original_sequence = proforma_sequence
        self.cursor = 0
        self.length = len(proforma_sequence)

        # --- Current Parse State Attributes ---
        self.amino_acids: list[str] = []

        # Identifiers
        self.compound_name: str | None = None  # (>>>Name)
        self.ion_name: str | None = None  # (>>Name)
        self.peptide_name: str | None = None  # (>Name)

        # Modifications
        self.global_mods: dict[str, int] | None = None  # <Mod>
        self.labile_mods: dict[str, int] | None = None  # {Mod}
        self.nterm_mods: dict[str, int] | None = None  # [Mod]-
        self.cterm_mods: dict[str, int] | None = None  # -[Mod]
        self.internal_mods: dict[int, dict[str, int]] | None = None  # S[Mod]
        self.unknown_mods: dict[str, int] | None = None  # [Mod]?
        self.intervals: list[Interval] = []

        # Properties
        self.charge: int | None = None
        self.charge_adducts: dict[str, int] | None = None
        self.is_chimeric: bool = False
        self.is_crosslinked: bool = False

    def _peek_char(self) -> str:
        """Get current character without bounds checking."""
        return self.original_sequence[self.cursor]

    def _advance(self, n: int = 1) -> None:
        """Advance cursor without bounds checking."""
        self.cursor += n

    def parse(self) -> Generator[tuple["ProFormaParser", bool | None], None, None]:
        """
        Main Entry Point.
        Parses a Compound Peptidoform Ion string.
        Yields (self, connection_type) where connection_type is:
           False: Chimeric connection (+)
           True:  Crosslink connection (//)
           None:  End of chain
        """
        self.cursor = 0

        if self.length == 0:
            yield self, None
            return

        # 1. Parse Compound Header (>>>Name) and <Global Mods>
        # These apply to the entire string context
        self._parse_compound_header()

        while self.cursor < self.length:
            # We are now parsing a specific Peptidoform Ion (chain or set of crosslinks)
            # This handles the Chimeric split (+)

            # We must buffer crosslinked peptides (A//B/2) because the charge
            # at the end applies to A and B.
            peptidoform_group = self._parse_peptidoform_ion_group()

            # Determine if there are more chimeric parts coming
            has_more_chimeras = False
            if self.cursor < self.length and self.original_sequence[self.cursor] == "+":
                self.cursor += 1
                has_more_chimeras = True

            # Yield the buffered group
            for i, parser_state in enumerate(peptidoform_group):
                # Determine connection type for the yield
                if i < len(peptidoform_group) - 1:
                    connection = True  # Connected by //
                else:
                    connection = False if has_more_chimeras else None

                # Apply the global context we parsed earlier to this instance
                parser_state.compound_name = self.compound_name
                if self.global_mods:
                    if parser_state.global_mods is None:
                        parser_state.global_mods = {}
                    for mod, count in self.global_mods.items():
                        parser_state.global_mods[mod] = parser_state.global_mods.get(mod, 0) + count

                yield parser_state, connection

    # =========================================================================
    # Level 1: Compound Ion & Headers
    # =========================================================================

    def _parse_compound_header(self):
        """Parses (>>>Name) and <GlobalMods>"""
        # Check for Compound Name (>>>Name)
        if self._peek_startswith("(>>>"):
            compound_name = self._parse_header_name(">>>")

            if compound_name.startswith(">"):
                self._raise_parse_error(
                    "Invalid compound name: cannot start with '>'",
                    self.cursor - len(compound_name) - 2,
                )

            self.compound_name = compound_name

        # Check for Global Modifications <Mod>
        # Two types:
        # 1. Isotope modifications: <13C>, <15N>
        # 2. Static modifications: <[Thr->Xle]@V>
        while self.cursor < self.length and self.original_sequence[self.cursor] == "<":
            mod_start = self.cursor
            self.cursor += 1  # Skip

            # Parse content until matching >
            # Need to track nested brackets [ ] to handle <[...]@X>
            depth = 0
            content_start = self.cursor

            while self.cursor < self.length:
                c = self.original_sequence[self.cursor]
                if c == "[":
                    depth += 1
                elif c == "]":
                    depth -= 1
                elif c == ">" and depth == 0:
                    # Found the closing > for the global mod
                    break
                self.cursor += 1

            if self.cursor >= self.length:
                self._raise_parse_error("Unclosed global modification bracket", mod_start)

            content = self.original_sequence[content_start : self.cursor]
            self.cursor += 1  # Skip closing >

            if self.global_mods is None:
                self.global_mods = {}

            mod_key = sys.intern(content)
            self.global_mods[mod_key] = self.global_mods.get(mod_key, 0) + 1

    def _parse_peptidoform_ion_group(self) -> list["ProFormaParser"]:
        """
        Parses a group of peptides connected by // (Crosslinks).
        Handles the shared charge state at the end.
        Returns a list of populated parser objects.
        """
        group: list[ProFormaParser] = []

        # 1. Parse Peptides separated by //
        while True:
            # Create a new parser instance for this specific peptide
            # (We assume the caller will inject global state later)
            p = ProFormaParser("")
            # We hijack the internal state of 'p' using our cursor
            self._parse_single_peptidoform_into(p)
            group.append(p)

            if self.cursor < self.length and self._peek_startswith("//"):
                self.cursor += 2  # Skip //
                p.is_crosslinked = True
                continue
            break

        # 2. Parse Shared Charge State (Applies to the whole group)
        # The charge sits at the end of the last peptide in the group
        charge, adducts = self._parse_charge_state()

        # 3. Distribute charge info to all peptides in the crosslink group
        for p in group:
            p.charge = charge
            p.charge_adducts = adducts

        return group

    # =========================================================================
    # Level 2: Single Peptidoform Parsing
    # =========================================================================

    def _parse_single_peptidoform_into(self, target: "ProFormaParser"):
        """
        Parses a single linear peptide sequence into a target Parser object.
        State: Expects cursor to be at the start of a peptide.
        Stops at: //, +, or End of String.
        """

        # 1. Parse Ion Name (>>Name)
        if self._peek_startswith("(>>"):
            ion_name = self._parse_header_name(">>")

            if ion_name.startswith(">"):
                self._raise_parse_error(
                    "Invalid ion name: cannot start with '>'",
                    self.cursor - len(ion_name) - 2,
                )

            target.ion_name = ion_name

        # 2. Parse Peptide Name (>Name)
        if self._peek_startswith("(>"):
            peptide_name = self._parse_header_name(">")

            if peptide_name.startswith(">"):
                self._raise_parse_error(
                    "Invalid peptide name: cannot start with '>'",
                    self.cursor - len(peptide_name) - 2,
                )

            target.peptide_name = peptide_name

        # Check for empty header ()
        if self._peek_startswith("()"):
            self._raise_parse_error("Expected 1-3 '>' characters")

        # 3. Parse Prefix: Unlocalized Mods [?], Labile {}, N-Term [ ]-
        self._parse_sequence_prefixes(target)

        # 4. Parse Amino Acid Sequence and Inline Mods
        self._parse_sequence_body(target)

        # 5. Parse C-Terminal Mods -[ ]
        # Note: We must distinguish C-term mods from C-term Sequence.
        # C-term mods must be preceded by '-' and followed by End, /, //, or +
        if self.cursor < self.length and self.original_sequence[self.cursor] == "-":
            # Lookahead to ensure this is a mod, not part of a sequence like A-B (if that were allowed)
            # ProForma 2.1 requires -[Mod] for C-term
            if self.cursor + 1 < self.length and self.original_sequence[self.cursor + 1] == "[":
                self.cursor += 1  # Consume '-'
                mods = self._parse_bracket_content("[", "]", allow_multiplier=False)
                for mod in mods:
                    if target.cterm_mods is None:
                        target.cterm_mods = {}
                    target.cterm_mods[mod] = target.cterm_mods.get(mod, 0) + 1

    def _parse_sequence_prefixes(self, target: "ProFormaParser"):
        """Parses [Mod]?, {Mod}, and [Mod]-"""
        while self.cursor < self.length:
            char = self.original_sequence[self.cursor]

            # Check for stop characters (Start of sequence or end of parsing)
            if char in VALID_AMINO_ACIDS or char == "(":
                # '(' could be start of interval (Range) or (ambiguous AA)
                if self._is_start_of_sequence_char():
                    break

            if char == "{":
                # Labile Mod
                mods = self._parse_bracket_content("{", "}", allow_multiplier=False)
                for mod in mods:
                    if target.labile_mods is None:
                        target.labile_mods = {}
                    target.labile_mods[mod] = target.labile_mods.get(mod, 0) + 1

            elif char == "[":
                # Could be Unlocalized [Mod]? or N-Term [Mod]-
                # Need to lookahead to determine if multipliers are allowed
                start_pos = self.cursor

                # Peek ahead to find the suffix (- or ?)
                temp_cursor = self.cursor
                depth = 0
                while temp_cursor < self.length:
                    c = self.original_sequence[temp_cursor]
                    if c == "[":
                        depth += 1
                    elif c == "]":
                        depth -= 1
                        if depth == 0:
                            temp_cursor += 1
                            # Check for multiplier
                            if temp_cursor < self.length and self.original_sequence[temp_cursor] == "^":
                                while temp_cursor < self.length and (
                                    self.original_sequence[temp_cursor].isdigit() or self.original_sequence[temp_cursor] == "^"
                                ):
                                    temp_cursor += 1
                            # Now check what follows
                            if temp_cursor < self.length and self.original_sequence[temp_cursor] == "[":
                                # Another bracket follows, continue
                                continue
                            break
                    temp_cursor += 1

                # Determine suffix type
                is_unknown = temp_cursor < self.length and self.original_sequence[temp_cursor] == "?"
                is_nterm = temp_cursor < self.length and self.original_sequence[temp_cursor] == "-"

                # Parse with appropriate multiplier setting
                mods = self._parse_bracket_content("[", "]", allow_multiplier=is_unknown)

                if is_nterm:
                    # N-Terminal
                    self.cursor += 1  # Skip -
                    for mod in mods:
                        if target.nterm_mods is None:
                            target.nterm_mods = {}
                        target.nterm_mods[mod] = target.nterm_mods.get(mod, 0) + 1
                elif is_unknown:
                    # Unlocalized / Unknown position
                    self.cursor += 1  # Skip ?
                    for mod in mods:
                        if target.unknown_mods is None:
                            target.unknown_mods = {}
                        target.unknown_mods[mod] = target.unknown_mods.get(mod, 0) + 1
                else:
                    # If neither - nor ?, this might be an orphaned mod or error.
                    self._raise_parse_error(
                        "Modification must be followed by '-' (N-term) or '?' (Unknown pos)",
                        start_pos,
                    )

            elif char == "?" and self.cursor + 1 < self.length and self.original_sequence[self.cursor + 1] == "[":
                # Handling rare case if ? comes before (unlikely in standard ProForma 2.1 but robust to check)
                pass
            else:
                # If we hit something else, assume sequence start if valid, else error
                if char in VALID_AMINO_ACIDS:
                    break
                # It might be a format error or we reached valid sequence start
                break

    def _get_context_snippet(self, position: int, window: int = 20) -> str:
        """Get a snippet of the sequence around the error position"""
        start = max(0, position - window)
        end = min(self.length, position + window)
        snippet = self.original_sequence[start:end]

        # Calculate where the pointer should go
        pointer_offset = position - start
        pointer = " " * pointer_offset + "^"

        return f"{snippet}\n{pointer}"

    def _raise_parse_error(self, message: str, position: int | None = None) -> NoReturn:
        """Raise a ValueError with position context"""
        if position is None:
            position = self.cursor

        context = self._get_context_snippet(position)
        raise ValueError(f"{message}\nPosition {position}:\n{context}")

    def _parse_sequence_body(self, target: "ProFormaParser"):
        """Iterate over amino acids, intervals, and internal mods - optimized version"""
        seq = self.original_sequence
        length = self.length

        while self.cursor < length:
            char = seq[self.cursor]

            # Fast terminator check
            if char in _TERMINATOR_SET:
                break
            if char == "-" and self._peek_is_cterm():
                break

            # Standard Amino Acid - hot path
            if char in _VALID_AA_SET:
                target.amino_acids.append(char)
                self.cursor += 1
                # Only check for mods if '[' follows
                if self.cursor < length and seq[self.cursor] == "[":
                    self._parse_inline_mods(target, len(target.amino_acids) - 1)
            elif char == "(":
                self._parse_interval(target)
            else:
                hint = " (amino acids must be uppercase)" if char.islower() else ""
                self._raise_parse_error(
                    f"Unexpected character '{char}'{hint}: expected an uppercase amino acid (A-Z) or a ProForma token "
                    "('[' for a modification, '(' for an interval, '/' for charge, or '-' for a terminal modification)"
                )

    def _parse_interval(self, target: "ProFormaParser"):
        """Parses (StartSeq-EndSeq), (Seq), or (?Seq)"""
        start_pos = len(target.amino_acids)
        self.cursor += 1  # Skip (

        # Check for ambiguity flag '?' at the start of the interval
        # e.g., (?DQ) or (?EPTI)
        is_ambiguous = False
        if self.cursor < self.length and self.original_sequence[self.cursor] == "?":
            is_ambiguous = True
            self.cursor += 1  # Consume the ?

        # Parse the content of the interval
        while self.cursor < self.length and self.original_sequence[self.cursor] != ")":
            char = self.original_sequence[self.cursor]

            if char in VALID_AMINO_ACIDS:
                target.amino_acids.append(char)
                self.cursor += 1
                # CRITICAL FIX: Check for mods on this specific AA inside the interval
                # e.g. The [12] in P(E[12]PTI)DE
                self._parse_inline_mods(target, len(target.amino_acids) - 1)
            else:
                # If we hit something invalid inside parens (that isn't a mod caught above)
                hint = " (amino acids must be uppercase)" if char.islower() else ""
                self._raise_parse_error(
                    f"Unexpected character '{char}' inside interval '(...)'{hint}: expected an uppercase amino acid (A-Z) or ')' to close the interval"
                )

        if self.cursor >= self.length:
            interval_start = self.cursor
            # Find the start of the interval for better error reporting
            temp = start_pos
            while temp > 0 and self.original_sequence[temp - 1] != "(":
                temp -= 1
            if temp > 0:
                interval_start = temp - 1
            self._raise_parse_error("Unclosed interval parenthesis", interval_start)

        self.cursor += 1  # Skip )
        end_pos = len(target.amino_acids)

        # Check for modifications applied to the ENTIRE interval
        # e.g. The [Oxidation] in P(EPTI)[Oxidation]DE
        interval_mods: dict[str, int] | None = None
        while self.cursor < self.length and self.original_sequence[self.cursor] == "[":
            mods = self._parse_bracket_content("[", "]", allow_multiplier=False)
            for mod in mods:
                if interval_mods is None:
                    interval_mods = {}
                interval_mods[mod] = interval_mods.get(mod, 0) + 1

        # Note: Some ProForma versions allow '?' at the end for unlocalized ranges
        # If encountered here, we mark the interval as ambiguous
        if self.cursor < self.length and self.original_sequence[self.cursor] == "?":
            is_ambiguous = True
            self.cursor += 1

        target.intervals.append(Interval(start_pos, end_pos, is_ambiguous, interval_mods))

    def _parse_inline_mods(self, target: "ProFormaParser", aa_index: int):
        """Checks for [Mod] immediately following an AA"""
        while self.cursor < self.length and self.original_sequence[self.cursor] == "[":
            mods = self._parse_bracket_content("[", "]", allow_multiplier=False)
            if target.internal_mods is None:
                target.internal_mods = {}
            for mod in mods:
                if aa_index not in target.internal_mods:
                    target.internal_mods[aa_index] = {}
                target.internal_mods[aa_index][mod] = target.internal_mods[aa_index].get(mod, 0) + 1

    # =========================================================================
    # Level 3: Charge and Utilities
    # =========================================================================

    # Optimization 6: Faster charge state parsing
    def _parse_charge_state(self) -> tuple[int | None, dict[str, int] | None]:
        """Optimized charge parsing"""
        if self.cursor >= self.length or self.original_sequence[self.cursor] != "/":
            return None, None

        slash_pos = self.cursor
        self.cursor += 1  # Skip /
        if self.cursor >= self.length:
            self._raise_parse_error("Expected a charge value after '/' (e.g. '/2' or '/2[+H]')", slash_pos)

        seq = self.original_sequence
        char = seq[self.cursor]

        # Adduct case
        if char == "[":
            adducts = {}
            self.cursor += 1
            start = self.cursor
            depth = 1

            while self.cursor < self.length:
                c = seq[self.cursor]
                if c == "[":
                    depth += 1
                elif c == "]":
                    depth -= 1
                    if depth == 0:
                        break
                self.cursor += 1

            if depth > 0:
                self._raise_parse_error("Unclosed adduct bracket", start)

            content = seq[start : self.cursor]
            self.cursor += 1

            # Parse adducts
            for part in content.split(","):
                part = part.strip()
                if "^" in part:
                    base, mult = part.rsplit("^", 1)
                    count = int(mult) if mult.isdigit() else 1
                    part = base
                else:
                    count = 1

                key = sys.intern(part)
                adducts[key] = adducts.get(key, 0) + count

            return None, adducts

        # Integer charge case
        start = self.cursor
        if char in _SIGN_SET:
            self.cursor += 1

        while self.cursor < self.length and seq[self.cursor] in _DIGIT_SET:
            self.cursor += 1

        try:
            charge = int(seq[start : self.cursor])
            return charge, None
        except ValueError:
            self._raise_parse_error(
                f"Invalid charge after '/': expected an integer (optionally signed) or an adduct in '[...]', got {seq[slash_pos + 1 :] or '(nothing)'!r}",
                slash_pos,
            )

    # Optimization 7: Reduce property overhead
    @property
    def unmod_sequence(self) -> str:
        """Cache this if called multiple times"""
        return "".join(self.amino_acids)

    # =========================================================================
    # Helpers
    # =========================================================================

    def _parse_header_name(self, prefix_str: str) -> str:
        """Parses (>Name) style headers. Assumes cursor is at (."""
        start_pos = self.cursor
        # Skip the (>>> part
        self.cursor += 1 + len(prefix_str)
        content_start = self.cursor
        content = self._read_until_balanced("(", ")")

        if self.cursor >= self.length:
            # No balanced ')' closes the name -- e.g. a name containing a stray unbalanced
            # '(' like '(>a(b)PEPTIDE'. Fall back to the first ')' so such names stay
            # parseable (the pre-3.1.2 behavior); only a name with no ')' at all is a real
            # error. Balanced names like '(>my (special) peptide)' still take the branch above.
            self.cursor = content_start
            content = self._read_until_first(")")
            if self.cursor >= self.length:
                self._raise_parse_error("Unmatched '(' for name", start_pos)

        self.cursor += 1  # Skip )

        if content.startswith(">"):
            self._raise_parse_error("Expected 1-3 '>' characters", start_pos)

        return sys.intern(content)

    def _parse_bracket_content(self, open_char: str, close_char: str, allow_multiplier: bool) -> list[str]:
        """Optimized bracket parsing - reduces string operations"""
        items: list[str] = []
        seq = self.original_sequence
        length = self.length

        while self.cursor < length and seq[self.cursor] == open_char:
            self.cursor += 1  # Skip open
            start = self.cursor
            depth = 1

            # Fast bracket matching
            while self.cursor < length:
                c = seq[self.cursor]
                if c == open_char:
                    depth += 1
                elif c == close_char:
                    depth -= 1
                    if depth == 0:
                        break
                self.cursor += 1

            if depth != 0:
                # Reached end of sequence without a matching close bracket.
                self._raise_parse_error(
                    f"Unclosed '{open_char}': reached end of sequence before finding a matching '{close_char}'",
                    start - 1,
                )

            content = sys.intern(seq[start : self.cursor])
            if not content:
                self._raise_parse_error(
                    f"Empty modification '{open_char}{close_char}' is not allowed; put a modification name, accession, formula, or mass inside the brackets",
                    start - 1,
                )
            self.cursor += 1  # Skip close

            # Check for multiplier
            multiplier = 1
            if self.cursor < length and seq[self.cursor] == "^":
                m_start = self.cursor
                self.cursor += 1

                # Fast digit parsing
                while self.cursor < length and seq[self.cursor] in _DIGIT_SET:
                    self.cursor += 1

                if self.cursor > m_start + 1:
                    multiplier = int(seq[m_start + 1 : self.cursor])
                    if not allow_multiplier:
                        self._raise_parse_error(
                            "Multipliers not allowed for this bracketed content",
                            m_start,
                        )

            # Extend list once instead of repeated appends
            if multiplier == 1:
                items.append(content)
            else:
                items.extend([content] * multiplier)

        return items

    def _read_until_first(self, terminator: str) -> str:
        """Scan to the first ``terminator`` without tracking nesting; used as a lenient
        fallback for names whose parentheses don't balance."""
        start = self.cursor
        seq = self.original_sequence
        while self.cursor < self.length and seq[self.cursor] != terminator:
            self.cursor += 1
        return seq[start : self.cursor]

    def _read_until_balanced(self, open_char: str, close_char: str) -> str:
        """Scan for the matching ``close_char``, tracking nesting depth so balanced
        occurrences of ``open_char``/``close_char`` inside the content don't
        terminate the scan early. Assumes one ``open_char`` has already been consumed."""
        start = self.cursor
        depth = 1
        seq = self.original_sequence
        while self.cursor < self.length:
            c = seq[self.cursor]
            if c == open_char:
                depth += 1
            elif c == close_char:
                depth -= 1
                if depth == 0:
                    break
            self.cursor += 1
        return seq[start : self.cursor]

    def _peek_startswith(self, s: str) -> bool:
        return self.original_sequence.startswith(s, self.cursor)

    def _peek_is_cterm(self) -> bool:
        """
        Check if current position is start of C-term mod (-[Mod]).
        Used to distinguish from sequences if hyphens were allowed (they aren't really, but good for safety).
        """
        if self.cursor + 1 < self.length and self.original_sequence[self.cursor + 1] == "[":
            return True
        return False

    def _is_start_of_sequence_char(self) -> bool:
        """
        Determines if the current char is the start of the AA sequence.
        This effectively ends the Prefix parsing loop.
        """
        c = self.original_sequence[self.cursor]
        if c in VALID_AMINO_ACIDS:
            return True
        if c == "(":
            return True  # Range or Ambiguous
        return False

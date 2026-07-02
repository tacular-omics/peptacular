from __future__ import annotations

import re
from collections.abc import Generator, Sequence

from tacular import PROTEASE_LOOKUP

from ..regex_utils import get_regex_match_indices
from ..spans import (
    Span,
    build_left_semi_spans,
    build_non_enzymatic_spans,
    build_right_semi_spans,
    build_spans,
)
from .types import (
    DigestProtocol,
    EnzymeConfig,
)


def left_semi_spans(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[Span]:
    """Get left semi-enzymatic sequences from a ProForma annotation."""
    span = Span(0, len(annotation), 0)
    return build_left_semi_spans(span=span, min_len=min_len, max_len=max_len)


def right_semi_spans(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[Span]:
    """Get right semi-enzymatic sequences from a ProForma annotation."""
    span = Span(0, len(annotation), 0)
    return build_right_semi_spans(span=span, min_len=min_len, max_len=max_len)


def semi_spans(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[Span]:
    """Get both left and right semi-enzymatic sequences from a ProForma annotation."""
    yield from left_semi_spans(annotation, min_len=min_len, max_len=max_len)
    yield from right_semi_spans(annotation, min_len=min_len, max_len=max_len)


def nonspecific_spans(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[Span]:
    """Get non-enzymatic sequences from a ProForma annotation."""
    span = Span(0, len(annotation), 0)
    return build_non_enzymatic_spans(span=span, min_len=min_len, max_len=max_len)


def get_cleavage_sites(annotation: DigestProtocol, enzyme: str | re.Pattern[str]) -> Generator[int]:
    """Get cleavage sites for a given enzyme (name, regex string, or compiled pattern)."""

    # Normalize to compiled pattern
    if isinstance(enzyme, str):
        # Try to look up by name first
        protease_info = PROTEASE_LOOKUP.get(enzyme)
        pattern = protease_info.pattern if protease_info else re.compile(enzyme)
    else:
        pattern = enzyme

    # Handle non-specific cleavage. An empty pattern (from an empty/None cleave_on)
    # is treated the same as the explicit non-specific sentinel "()".
    if pattern.pattern in ("()", ""):
        return (i for i in range(len(annotation.stripped_sequence) + 1))

    return get_regex_match_indices(input_str=annotation.stripped_sequence, regex_str=pattern)


def _convert_to_aa_set(aa_keys: str | None) -> set[str]:
    X = set("GASPVTCILJNDQKEMHFRYWUO")
    B = set("DN")
    J = set("IL")
    Z = set("EQ")

    def make_aa_set(aa_keys: str) -> set[str]:
        aas: set[str] = set()
        for char in aa_keys:
            if char == "X":
                aas.update(X)
            elif char == "B":
                aas.update(B)
            elif char == "J":
                aas.update(J)
            elif char == "Z":
                aas.update(Z)
            else:
                aas.add(char)

        return aas

    if aa_keys is None:
        return set()

    if "-" in aa_keys:
        if aa_keys.count("-") != 1:
            raise ValueError("Amino acid keys with '-' must contain exactly one '-' character.")
        first_part, second_part = aa_keys.split("-")
        return make_aa_set(first_part) - make_aa_set(second_part)

    return make_aa_set(aa_keys)


def generate_regex(
    cleave_on: str | None = None,
    restrict_before: str | None = None,
    restrict_after: str | None = None,
    cterminal: bool = True,
) -> re.Pattern[str]:
    """Generate regex for proteolytic cleavage sites."""

    if cleave_on is None:
        return re.compile("")

    if not cleave_on:
        return re.compile("")

    # Escape special regex characters
    escaped_cleave = "".join(re.escape(char) for char in _convert_to_aa_set(cleave_on))

    parts: list[str] = []

    if cterminal:
        # Cleave after residue: lookbehind for cleavage residue
        parts.append(f"(?<=[{escaped_cleave}])")

        # Restrict after: what follows cleavage point
        if restrict_after:
            escaped_restrict = "".join(re.escape(char) for char in _convert_to_aa_set(restrict_after))
            parts.append(f"(?=[^{escaped_restrict}])")

        # Restrict before: what precedes cleavage residue
        if restrict_before:
            escaped_restrict = "".join(re.escape(char) for char in _convert_to_aa_set(restrict_before))
            # Need to check what comes before the cleavage residue
            parts.insert(0, f"(?<=[^{escaped_restrict}][{escaped_cleave}])")
            # Remove the simple lookbehind since we now have the combined one
            parts = [parts[0]] + parts[2:]
    else:
        # Cleave before residue: lookahead for cleavage residue
        parts.append(f"(?=[{escaped_cleave}])")

        # Restrict before: what precedes cleavage point
        if restrict_before:
            escaped_restrict = "".join(re.escape(char) for char in _convert_to_aa_set(restrict_before))
            parts.insert(0, f"(?<=[^{escaped_restrict}])")

        # Restrict after: what follows cleavage residue
        if restrict_after:
            escaped_restrict = "".join(re.escape(char) for char in _convert_to_aa_set(restrict_after))
            # Replace simple lookahead with one that checks after the residue
            parts[-1] = f"(?=[{escaped_cleave}][^{escaped_restrict}])"

    regex_str = "".join(parts)

    return re.compile(regex_str)


def digest_annotation_by_aa(
    annotation: DigestProtocol,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    complete_digestion: bool = True,
    sort_output: bool = True,
) -> Generator[Span]:
    """Digest annotation by amino acid cleavage rules."""
    return digest_annotation_by_regex(
        annotation=annotation,
        enzyme_regex=generate_regex(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        ),
        missed_cleavages=missed_cleavages,
        semi=semi,
        min_len=min_len,
        max_len=max_len,
        complete_digestion=complete_digestion,
        sort_output=sort_output,
    )


def digest_annotation_by_regex(
    annotation: DigestProtocol,
    enzyme_regex: str | re.Pattern[str],
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
) -> Generator[Span]:
    """Digest annotation using a regex pattern."""
    all_spans: set[Span] = set()
    if not complete_digestion:
        all_spans.add(Span(0, len(annotation), 0))

    cleavage_sites_list: list[int] = list(get_cleavage_sites(annotation, enzyme=enzyme_regex))

    spans = build_spans(
        max_index=len(annotation),
        enzyme_sites=cleavage_sites_list,
        missed_cleavages=missed_cleavages,
        min_len=min_len,
        max_len=max_len,
        semi=semi,
    )
    all_spans.update(spans)

    if sort_output:
        sorted_spans = sorted(all_spans, key=lambda x: (x[0], x[1], x[2]))
        return (span for span in sorted_spans)

    return (span for span in all_spans)


def sequential_digest_annotation(
    annotation: DigestProtocol,
    enzyme_configs: Sequence[EnzymeConfig],
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[Span]:
    """Perform sequential digestion with multiple enzymes."""
    digested_spans: list[Span] = []

    for i, enzyme_config in enumerate(enzyme_configs):
        if i == 0:
            digested_spans = list(
                digest_annotation_by_regex(
                    annotation=annotation,
                    enzyme_regex=enzyme_config.enzyme_regex,
                    missed_cleavages=enzyme_config.missed_cleavages,
                    semi=enzyme_config.semi_enzymatic,
                    min_len=min_len,
                    max_len=None,
                    sort_output=True,
                    complete_digestion=enzyme_config.complete_digestion,
                )
            )
        else:
            if len(digested_spans) == 0:
                break

            sequential_digested_spans: list[Span] = []

            for span in digested_spans:
                sub_annotation = annotation.slice(span[0], span[1], inplace=False)

                _digested_spans = list(
                    digest_annotation_by_regex(
                        annotation=sub_annotation,
                        enzyme_regex=enzyme_config.enzyme_regex,
                        missed_cleavages=enzyme_config.missed_cleavages,
                        semi=enzyme_config.semi_enzymatic,
                        min_len=min_len,
                        max_len=None,
                        sort_output=True,
                        complete_digestion=enzyme_config.complete_digestion,
                    )
                )

                # Adjust spans to be relative to original annotation
                for digested_span in _digested_spans:
                    fixed_span = Span(
                        span[0] + digested_span[0],
                        span[0] + digested_span[1],
                        span[2],
                    )
                    sequential_digested_spans.append(fixed_span)

            digested_spans = sequential_digested_spans

    if max_len is not None:
        digested_spans = [span for span in digested_spans if span[1] - span[0] <= max_len]

    return (span for span in digested_spans)

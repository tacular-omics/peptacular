"""Every digest input form agrees, and the plain-sequence fast path agrees with slicing (5.0).

Slicing an annotation keeps the parent's charge and its labile and unlocalised
modifications; digesting drops them, because a digestion product is a different molecule.
So the slice reference holds for every input that carries none of those three, and
``ANALYTE_DECORATED`` covers the inputs where the two contracts part company.
"""

from dataclasses import dataclass

import pytest
from hypothesis import given
from hypothesis import strategies as st

import peptacular as pt
from peptacular.sequence.digestion import _digest_input

SEQUENCES = [
    "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK",
    "PEPTIDEKAAARPPKRRKK",
    "K",
    "R",
    "PEPTIDE",
    "UOXBZJKRPEPK",
]
DECORATED = [
    "PEP[Phospho]TIDEKAAARPPK",
    "[Acetyl]-PEPTIDEKAAAR",
    "PEPTIDEKAAAR-[Amidated]",
    "(>name)PEPTIDEKAAAR",
    "<13C>PEPTIDEKAAAR",
    "<[Carbamidomethyl]@C>PEPCTIDEKAAAR",
    "PEP(TIDEK)[Phospho]AAAR",
]


def _reference(sequence: str, kind: str) -> list[tuple[str, pt.Span]]:
    annot = pt.parse(sequence)
    match kind:
        case "digest":
            spans = annot.digest_spans(enzyme="trypsin", missed_cleavages=2, min_len=2, max_len=30)
        case "semi":
            spans = annot.digest_spans(enzyme="trypsin", missed_cleavages=1, semi=True)
        case "simple":
            spans = annot.simple_digest_spans(cleave_on="KR", restrict_before="P", missed_cleavages=1)
        case "left":
            spans = annot.left_semi_spans(min_len=2, max_len=8)
        case "right":
            spans = annot.right_semi_spans(min_len=2, max_len=8)
        case "semi_spans":
            spans = annot.semi_spans(min_len=2, max_len=8)
        case _:
            spans = annot.nonspecific_spans(min_len=2, max_len=6)
    return [(annot[span].serialize(), span) for span in spans]


def _functional(sequence, kind: str) -> list[tuple[str, pt.Span]]:
    match kind:
        case "digest":
            return pt.digest(sequence, "trypsin", missed_cleavages=2, min_len=2, max_len=30)
        case "semi":
            return pt.digest(sequence, "trypsin", missed_cleavages=1, semi=True)
        case "simple":
            return pt.simple_digest(sequence, cleave_on="KR", restrict_before="P", missed_cleavages=1)
        case "left":
            return pt.left_semi_digest(sequence, min_len=2, max_len=8)
        case "right":
            return pt.right_semi_digest(sequence, min_len=2, max_len=8)
        case "semi_spans":
            return pt.semi_digest(sequence, min_len=2, max_len=8)
        case _:
            return pt.nonspecific_digest(sequence, min_len=2, max_len=6)


# Annotations that belong to the parent analyte rather than to any one digestion product.
ANALYTE_DECORATED = [
    "PEPTIDEKAAAR/2",
    "PEPTIDEKAAAR/[Na:z+1]",
    "{Glycan:Hex}PEPTIDEKAAAR",
    "[Phospho]?PEPTIDEKAAAR",
]

KINDS = ["digest", "semi", "simple", "left", "right", "semi_spans", "nonspecific"]


@dataclass
class Entry:
    sequence: str


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("sequence", SEQUENCES + DECORATED)
def test_same_output_as_slicing(sequence, kind):
    try:
        expected = _reference(sequence, kind)
    except pt.PeptacularError as exc:
        # Unsliceable input (an interval cut by a span) fails the same way on every input type.
        for value in (sequence, pt.parse(sequence), Entry(sequence)):
            with pytest.raises(type(exc)):
                _functional(value, kind)
        return
    assert _functional(sequence, kind) == expected
    assert _functional(pt.parse(sequence), kind) == expected
    assert _functional(Entry(sequence), kind) == expected


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("sequence", ANALYTE_DECORATED)
def test_analyte_annotations_are_dropped_not_sliced(sequence, kind):
    """Digest differs from slicing by exactly the parent's charge, labile and unknown mods."""
    expected = []
    for peptide, span in _reference(sequence, kind):
        product = pt.parse(peptide)
        product.clear_charge(inplace=True)
        product.clear_labile_mods(inplace=True)
        product.clear_unknown_mods(inplace=True)
        expected.append((product.serialize(), span))
    assert _functional(sequence, kind) == expected
    assert _functional(pt.parse(sequence), kind) == expected
    assert _functional(Entry(sequence), kind) == expected


def test_plain_detection():
    annot, plain = _digest_input("PEPTIDEK")
    assert plain == "PEPTIDEK"
    assert annot.stripped_sequence == "PEPTIDEK"
    assert _digest_input(pt.parse("PEPTIDEK"))[1] == "PEPTIDEK"
    assert _digest_input(Entry("PEPTIDEK"))[1] == "PEPTIDEK"
    for decorated in DECORATED + ANALYTE_DECORATED:
        assert _digest_input(decorated)[1] is None


def test_plain_string_is_not_parsed_but_invalid_strings_still_raise():
    with pytest.raises(pt.ProFormaFormatError):
        pt.digest("pepTIDEK", "trypsin")
    with pytest.raises(pt.ProFormaFormatError):
        pt.digest("PEP[", "trypsin")


def test_list_input_matches_single():
    seqs = SEQUENCES + DECORATED + ANALYTE_DECORATED
    assert pt.digest(seqs, "trypsin", missed_cleavages=1) == [pt.digest(s, "trypsin", missed_cleavages=1) for s in seqs]


@given(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=60))
def test_property_plain_matches_slicing(sequence):
    for kind in ("digest", "semi", "simple"):
        assert _functional(sequence, kind) == _reference(sequence, kind)

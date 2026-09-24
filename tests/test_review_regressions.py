"""Regressions for bugs found in the post-4.1.0 code review."""

import re

import pytest

import peptacular as pt
from peptacular.digestion.core import generate_regex, get_cleavage_sites


@pytest.mark.parametrize(
    "sequence",
    [
        "<[Carbamidomethyl]@N-term>PEPTIDE",
        "<[Amidated]@C-term>PEPTIDE",
        "<[Carbamidomethyl]@C>PEPCTIDEC",
        "<[Acetyl]@N-term><[Amidated]@C-term><[Oxidation]@M>PEMTIDEM",
        "<[Acetyl]@N-term:P>PEPTIDE",
        "<[Acetyl]@N-term:K>PEPTIDE",
        "<[Amidated]@C-term>[Acetyl]-PEPTIDE-[Methyl]",
    ],
)
@pytest.mark.parametrize("ion", ["a", "b", "c", "x", "y", "z", "p"])
@pytest.mark.parametrize("charge", [1, 2])
def test_fast_fragment_matches_frag_with_static_mods(sequence, ion, charge):
    annot = pt.parse(sequence)
    fast = annot.fast_fragment(ion_types=[ion], charges=[charge])[(pt.IonType(ion), charge)]
    if ion == "p":
        expected = [annot.frag(ion_type=ion, charge=charge).mz] * len(annot)
    else:
        expected = [annot.frag(ion_type=ion, charge=charge, position=i).mz for i in range(1, len(annot) + 1)]
    assert fast == pytest.approx(expected, abs=1e-6)


def test_immonium_neutral_loss_sites_use_the_residue():
    # Only T can lose water; P immonium ions must not raise or get a water loss.
    frags = pt.parse("PPPPTPP").fragment(ion_types=["i"], charges=[1], neutral_deltas=["H2O"])
    lossy = [f for f in frags if f.losses]
    assert {f.position for f in lossy} == {5}


def test_internal_neutral_loss_sites_use_the_subsequence():
    frags = pt.parse("APPPPPPT").fragment(ion_types=["by"], charges=[1], neutral_deltas=["H2O"])
    assert frags
    assert not [f for f in frags if f.losses]


@pytest.mark.parametrize(
    ("sequence", "kwargs", "expected"),
    [
        ("KAAAKAAA", dict(cleave_on="K", restrict_before="P"), [1, 5]),
        ("PKAAKAAA", dict(cleave_on="K", restrict_before="P"), [5]),
        ("AAADAAAD", dict(cleave_on="D", restrict_after="P", cterminal=False), [3, 7]),
        ("AAADPAAD", dict(cleave_on="D", restrict_after="P", cterminal=False), [7]),
        ("KAAAKPAA", dict(cleave_on="K", restrict_after="P"), [1]),
        ("DAAPDAAD", dict(cleave_on="D", restrict_before="P", cterminal=False), [7]),
    ],
)
def test_generate_regex_restrictions_at_sequence_ends(sequence, kwargs, expected):
    annot = pt.parse(sequence)
    assert list(get_cleavage_sites(annot, generate_regex(**kwargs))) == expected


def _config(regex, mc):
    return pt.EnzymeConfig(enzyme=re.compile(regex), missed_cleavages=mc, semi_enzymatic=False, complete_digestion=True)


def test_sequential_digest_counts_second_enzyme_missed_cleavages():
    annot = pt.parse("AAKAADAA")
    spans = set(annot.sequential_digest_spans([_config("(?<=K)", 0), _config("(?=D)", 1)]))
    assert spans == {(0, 3, 0), (3, 5, 0), (3, 8, 1), (5, 8, 0)}


def test_sequential_digest_counts_only_first_enzyme_sites_inside_span():
    annot = pt.parse("AAKAADAA")
    spans = set(annot.sequential_digest_spans([_config("(?<=K)", 1), _config("(?=D)", 0)]))
    # (5, 8) comes from both the (0, 8) and (3, 8) parents; neither contains the K site
    assert spans == {(0, 3, 0), (0, 5, 1), (5, 8, 0), (3, 5, 0)}


def test_is_subsequence_unordered_missing_residue():
    assert pt.is_subsequence("W", "PEPTIDE", order=False) is False
    assert pt.is_subsequence("EP", "PEPTIDE", order=False) is True


def test_annotation_hash_is_order_independent():
    a = pt.parse("[Acetyl][Formula:C2]-PEPTIDE")
    b = pt.parse("[Formula:C2][Acetyl]-PEPTIDE")
    assert a == b
    assert hash(a) == hash(b)
    assert len({a, b}) == 1


def test_annotation_eq_with_other_type():
    annot = pt.parse("PEPTIDE")
    assert (annot == "PEPTIDE") is False
    assert annot != "PEPTIDE"
    assert annot not in ["PEPTIDE", 1, None]

"""mzPAF neutral-loss labels use the canonical names, and impossible losses are skipped."""

import json
from pathlib import Path

import pytest

import peptacular as pt

REFERENCE = Path(__file__).parent / "reference" / "mzpaf_neutral_losses_4_2_0.json"


@pytest.mark.parametrize(
    ("loss", "label"),
    [("NH3", "-NH3"), ("H2O", "-H2O"), ("H3PO4", "-H3PO4"), ("CO", "-CO"), ("HCONH2", "-HCONH2"), ("HCOOH", "-HCOOH"), ("HPO3", "-HPO3")],
)
def test_named_losses_use_mzpaf_names(loss, label):
    # tacular 2.0 stores compositions H-first (H3N); mzPAF forbids "H3N" for ammonia.
    frags = pt.fragment("S[Phospho]EQNKDR", ["b", "y"], [1], neutral_deltas=[loss])
    labels = {f.to_mzpaf(include_sequence=False) for f in frags}
    assert any(x.endswith(label) for x in labels)
    assert not any(x.endswith(("-H3N", "-H3CON", "-H2CO2", "-H3O4P", "-HO3P")) for x in labels)


def test_named_delta_given_as_a_delta_keeps_the_name():
    frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"NH3": 1})
    assert frag.to_mzpaf(include_sequence=False) == "b3-NH3"


@pytest.mark.parametrize(("delta", "label"), [("C2H2O", "+C2H2O"), ("OC2H2", "+C2H2O"), ("NO", "+NO"), ("SNa", "+NaS"), ("[13C2]H2", "+[13C2]H2")])
def test_unnamed_formulas_use_hill_order(delta, label):
    # a plain formula is a gain, written in Hill order (C, H, then alphabetical)
    frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={delta: 1})
    assert frag.to_mzpaf(include_sequence=False) == f"b3{label}"


def test_labels_match_4_2_0():
    data = json.loads(REFERENCE.read_text())
    renamed = {"-H3CON": "-HCONH2", "-H2CO2": "-HCOOH"}  # 4.2.0 wrote these in composition order
    assert len(data["labels"]) >= 30
    for key, expected in data["labels"].items():
        peptide, loss = key.split("|")
        for old, new in renamed.items():
            expected = [label.replace(old, new) for label in expected]
        labels = sorted(f.to_mzpaf() for f in pt.fragment(peptide, ["b", "y"], [1, 2], neutral_deltas=[loss]))
        assert labels == expected, key


@pytest.mark.parametrize("loss", ["H3PO4", "HPO3", "SO3"])
def test_impossible_loss_is_skipped_not_raised(loss):
    # 4.x aborted the whole call: the loss is offered on any S/T/Y, but an unmodified residue
    # has no phosphorus or sulfur to lose.
    frags = pt.fragment("PEPSTYDEK", ["b", "y"], [1], neutral_deltas=[loss])
    assert len(frags) == len(pt.fragment("PEPSTYDEK", ["b", "y"], [1]))
    assert all(not f.deltas for f in frags)


def test_possible_loss_still_produced_next_to_impossible_ones():
    frags = pt.fragment("PEPS[Phospho]TIDEK", ["b"], [1], neutral_deltas=["H3PO4"])
    labels = {f.to_mzpaf(include_sequence=False) for f in frags}
    assert "b4-H3PO4" in labels and "b3-H3PO4" not in labels


def test_user_delta_that_is_impossible_still_raises():
    with pytest.raises(pt.InvalidAdjustmentError):
        pt.fragment("PEPTIDE", ["b"], [1], deltas=[{"H3PO4": 1}], calculate_with_composition=True)

"""A slice or digest must not separate a localization group's label from its primary tag (#85)."""

import pytest

import peptacular as pt
from peptacular.diagnostics import PeptacularError


def test_digest_splitting_localization_group_raises():
    with pytest.raises(PeptacularError, match="#g1"):
        pt.digest("EM[#g1(0.10)]EVTKS[Phospho#g1(0.90)]ESPEK", enzyme="trypsin", missed_cleavages=0)


def test_digest_splitting_localization_group_primary_in_first_half_raises():
    with pytest.raises(PeptacularError, match="#g1"):
        pt.digest("EM[Phospho#g1(0.90)]EVTKS[#g1(0.10)]ESPEK", enzyme="trypsin", missed_cleavages=0)


def test_slice_splitting_localization_group_raises():
    with pytest.raises(PeptacularError, match="#g1"):
        pt.parse("EM[#g1(0.10)]EVTKS[Phospho#g1(0.90)]ESPEK")[0:6]


def test_group_wholly_inside_one_product_digests_and_keeps_score():
    products = pt.digest("PEPS[Phospho#g1(0.90)]T[#g1(0.10)]IDEK", enzyme="trypsin", missed_cleavages=0)
    serialized = [seq for seq, _ in products]
    assert any("Phospho#g1(0.90)" in seq and "#g1(0.10)" in seq for seq in serialized)

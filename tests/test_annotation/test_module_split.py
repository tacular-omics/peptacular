"""The annotation split into private modules must not leak into pickles, reprs or errors."""

import pickle

import pytest

import peptacular as pt
from peptacular.annotation.annotation import ChargeType, ProFormaAnnotation


@pytest.mark.parametrize("member", list(ChargeType))
def test_charge_type_pickles_under_its_pre_split_module(member: ChargeType) -> None:
    assert ChargeType.__module__ == "peptacular.annotation.annotation"
    data = pickle.dumps(member)
    # older peptacular only has ChargeType in peptacular.annotation.annotation
    assert b"peptacular.annotation.annotation" in data
    assert b"_mod_access" not in data
    assert pickle.loads(data) is member


def test_charge_type_pickles_inside_a_container() -> None:
    annot = pt.parse("PEPTIDE/2")
    value = {"annot": annot, "charge_type": annot.charge_type}
    restored = pickle.loads(pickle.dumps(value))
    assert restored["charge_type"] is ChargeType.INT
    assert restored["annot"].serialize() == "PEPTIDE/2"


def test_mod_accessors_are_named_after_proforma_annotation() -> None:
    annot = pt.parse("PEM[Oxidation]TIDE")
    assert "_ModAccessMixin" not in repr(annot.strip_mods)
    assert ProFormaAnnotation.strip_mods.__qualname__ == "ProFormaAnnotation.strip_mods"
    assert ProFormaAnnotation.has_internal_mods.fget.__qualname__ == "ProFormaAnnotation.has_internal_mods"  # type: ignore[attr-defined]
    with pytest.raises(TypeError, match=r"^ProFormaAnnotation\.strip_mods\(\)"):
        annot.strip_mods(1, 2, 3)  # type: ignore

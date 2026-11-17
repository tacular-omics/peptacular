"""Fragmenting a peptide whose localization group spans two residues.

#95 added a guard that refuses a slice cutting a localization group, and #112 exempted
``split`` because a split keeps every residue. The fragment engine slices the annotation
once per ion for the same reason and needs the same exemption: without it the internal and
satellite ion types raised ``PeptacularError`` instead of returning ions, while a, b, c,
x, y, z and p went on working.

Every pair below was measured on v5.1.0 (``8253412``) and again after the fix, and the two
runs agree to the last digit.
"""

import pytest

import peptacular as pt

SEQ = "PEPT[Phospho#g1]IDEK[#g1]"

# ion type -> (number of ions, sum of the singly charged m/z rounded to 6 places).
V510 = {
    "a": (8, 4134.803022),
    "b": (8, 4358.762339),
    "c": (8, 4494.974732),
    "x": (8, 4913.968959),
    "y": (8, 4706.134843),
    "z": (8, 4569.922449),
    "i": (8, 773.509589),
    "by": (21, 7219.637874),
    "ay": (21, 6631.744664),
    "w": (5, 2300.065833),
    "d": (5, 2712.151559),
    "v": (8, 4156.770658),
    "p": (1, 1008.428535),
}

# The immonium m/z in full, because this is the ion type the review's reproducer used.
I_MZ_510 = (70.065126, 70.065126, 86.096426, 88.039305, 101.107325, 102.054955, 102.054955, 154.026371)


def _mz(ion_type: str) -> list[float]:
    return sorted(round(f.mz, 6) for f in pt.fragment(SEQ, ion_types=(ion_type,), charges=(1,)))


@pytest.mark.parametrize("ion_type", sorted(V510))
def test_every_ion_type_keeps_its_v510_ions(ion_type: str):
    expected_count, expected_sum = V510[ion_type]
    mz = _mz(ion_type)
    assert len(mz) == expected_count
    assert sum(mz) == pytest.approx(expected_sum, abs=1e-6)


def test_immonium_mz_values_are_unchanged():
    assert tuple(_mz("i")) == pytest.approx(I_MZ_510)


def test_fragment_arrays_and_fast_fragment_agree_with_fragment():
    assert len(pt.fragment_arrays(SEQ, ion_types=("i",), charges=(1,))["mz"]) == 8
    assert sum(len(v) for v in pt.fast_fragment(SEQ, charges=(1,)).values()) == len(pt.fragment(SEQ, charges=(1,)))


def test_the_default_ion_types_still_fragment():
    assert len(pt.fragment(SEQ, charges=(1,))) == 16


def test_a_user_slice_that_cuts_the_group_still_raises():
    # The exemption is for the engine's own per-ion slices, not for slicing the annotation.
    annot = pt.parse(SEQ)
    with pytest.raises(pt.PeptacularError, match="Localization group"):
        annot.slice(0, 5)
    with pytest.raises(pt.PeptacularError, match="Localization group"):
        _ = annot[0:5]
    with pytest.raises(pt.PeptacularError, match="Localization group"):
        pt.coverage(SEQ, ["PEPT[Phospho#g1]"])


def test_split_still_keeps_the_group():
    assert pt.join(pt.split(SEQ)) == SEQ

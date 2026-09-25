"""Residue composition via residue counts equals the per-residue sum it replaced (5.0)."""

from collections import Counter

import pytest
from hypothesis import given
from hypothesis import strategies as st
from tacular import AA_LOOKUP

import peptacular as pt

DEFINED = "".join(aa for aa, info in AA_LOOKUP.items() if len(aa) == 1 and info.composition is not None)


def _per_residue(sequence: str) -> Counter:
    total: Counter = Counter()
    for aa in sequence:
        for element, count in AA_LOOKUP[aa].composition.items():
            total[element] += count
    return total


@given(st.text(alphabet=DEFINED, min_size=0, max_size=80))
def test_matches_per_residue_sum(sequence):
    assert pt.ProFormaAnnotation(sequence=sequence).get_sequence_composition() == _per_residue(sequence)


def test_every_defined_residue():
    for aa in DEFINED:
        assert pt.ProFormaAnnotation(sequence=aa * 3).get_sequence_composition() == _per_residue(aa * 3)


@pytest.mark.parametrize("sequence", ["PEPTIDE", "[Acetyl]-PEPM[Oxidation]TIDEKC[Carbamidomethyl]LLSGR/2", "AS[Phospho]DFGHIK/3"])
def test_comp_and_mass_still_agree(sequence):
    annot = pt.parse(sequence).set_charge(None)
    comp = annot.comp()
    mass = sum(element.get_mass(monoisotopic=True) * count for element, count in comp.items())
    assert mass == pytest.approx(annot.mass(), abs=1e-6)


def test_undefined_residue_raises():
    with pytest.raises(pt.CompositionError):
        pt.ProFormaAnnotation(sequence="PEPBIDE").get_sequence_composition()

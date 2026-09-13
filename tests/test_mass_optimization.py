"""Scalar mass optimizations must preserve fragment and composition semantics."""

import pytest
from tacular import AA_LOOKUP, IonType

import peptacular as pt
from peptacular.annotation.cached_comps import _get_isotopes


@pytest.mark.parametrize(
    "sequence",
    [
        "PEPTIDEK",
        "ACDEFGHIKLMNPQRSTVWY" * 25,
        "[Acetyl]-PEPM[Oxidation]TIDEK-[Amidated]",
        "{Glycan:Hex}PEPTIDEK",
        "[Phospho]?PEPTIDEK",
        "PEP(TID)[Phospho]EK",
        "<[Carbamidomethyl]@C>ACDCCK",
        "PEP[Formula:Na:z+1]TIDEK",
        "PEP[Formula:H:z-1]TIDEK",
        "PEPTIDEK/3",
        "PEPTIDEK/[Na:z+1^2]",
        "PEPTIDEK[+12.345]",
    ],
)
@pytest.mark.parametrize("charge", [None, 0, 1, 2])
@pytest.mark.parametrize("monoisotopic", [True, False])
@pytest.mark.parametrize("ion_type", [IonType.PRECURSOR, IonType.NEUTRAL])
def test_scalar_matches_fragment(sequence, charge, monoisotopic, ion_type):
    annot = pt.parse(sequence)
    before = annot.serialize()
    kwargs = dict(charge=charge, monoisotopic=monoisotopic, ion_type=ion_type)
    fragment = annot.frag(**kwargs)
    assert annot.mass(**kwargs) == fragment.mass
    assert annot.mz(**kwargs) == fragment.mz
    if "+12.345" not in sequence:
        # Vocabulary masses have less precision than elemental reference masses.
        tolerance = 1e-6 if monoisotopic else 1e-3
        assert annot.mass(**kwargs) == pytest.approx(annot.mass(**kwargs, calculate_with_composition=True), abs=tolerance, rel=0)
    assert annot.serialize() == before


@pytest.mark.parametrize(
    "sequence,kwargs",
    [
        ("<13C>PEPTIDEK", {}),
        ("PEPTIDEK", {"isotopes": 2}),
        ("PEPTIDEK", {"deltas": {"H2O": -1}}),
        ("PEPTIDEK", {"deltas": 10.0}),
        ("PEPTIDEK", {"charge": -2}),
        ("PEPTIDEK", {"charge": "Na:z+1"}),
        ("{Glycan:Hex}PEPTIDEK", {"ion_type": "b", "charge": 2}),
        ("{Glycan:Hex}PEPTIDEK", {"ion_type": "y", "charge": 2}),
    ],
)
def test_fallback_matches_fragment(sequence, kwargs):
    annot = pt.parse(sequence)
    fragment = annot.frag(**kwargs)
    assert pt.mass(annot, **kwargs) == fragment.mass
    assert pt.mz(annot, **kwargs) == fragment.mz


@pytest.mark.parametrize("aa", list(AA_LOOKUP.one_letter_to_info))
@pytest.mark.parametrize("monoisotopic", [True, False])
def test_residue_lookup_preserves_reference_masses(aa, monoisotopic):
    annot = pt.ProFormaAnnotation(aa)
    expected = AA_LOOKUP.one_letter_to_info[aa].get_mass(monoisotopic)
    if expected is None:
        with pytest.raises(ValueError, match="Mass not available"):
            annot.mass(monoisotopic=monoisotopic)
    else:
        assert annot.mass(ion_type="n", monoisotopic=monoisotopic) == expected
        assert annot._build_mass_vector(monoisotopic) == [expected]


@pytest.mark.parametrize("func", [pt.mass, pt.mz, pt.comp])
def test_empty_sequence_has_clear_error(func):
    with pytest.raises(ValueError, match="empty sequence"):
        func("")


@pytest.mark.parametrize("charge", [True, False, 1.5])
def test_invalid_charge_still_rejected(charge):
    with pytest.raises(ValueError):
        pt.mass("PEPTIDE", charge=charge)


@pytest.mark.parametrize("count", [True, False, 1.0, 0.0])
@pytest.mark.parametrize("warm_cache", [False, True])
def test_invalid_isotope_count_rejected_independent_of_cache(count, warm_cache):
    _get_isotopes.cache_clear()
    if warm_cache:
        pt.mass("PEPTIDE", isotopes={"13C": int(count)})
    with pytest.raises(pt.InvalidAdjustmentError, match="Isotope count must be an integer"):
        pt.mass("PEPTIDE", isotopes={"13C": count})


def test_reused_annotation_observes_mutation_and_charge_override():
    annot = pt.parse("PEPTIDEK/[Na:z+1^2]")
    assert annot.neutral_mass() == pt.mass("PEPTIDEK")
    pt.mass(annot, charge=2)
    annot.sequence = "ACDEK"
    annot.set_internal_mods({1: "Oxidation"})
    assert pt.mass(annot, charge=2) == pt.mass("AC[Oxidation]DEK", charge=2)
    assert annot.charge_state == 2

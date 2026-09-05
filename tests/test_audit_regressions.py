"""Scientific consistency and mutation regressions from the repository audit."""

import pickle
from collections import Counter

import pytest

import peptacular as pt
from peptacular.annotation.cached_comps import ChargeCarrierInfo, DeltaInfo


@pytest.mark.parametrize("label", ["13C", "15N", "18O", "2H"])
@pytest.mark.parametrize("ion", ["p", "a", "b", "c", "x", "y", "z"])
@pytest.mark.parametrize("mono", [True, False])
def test_labeled_mass_agrees_with_composition(label, ion, mono):
    peptide = pt.parse(f"<{label}>PEM[Oxidation]TIDE")
    kwargs = dict(ion_type=ion, charge=2, monoisotopic=mono)
    assert peptide.mass(**kwargs) == pytest.approx(peptide.mass(**kwargs, calculate_with_composition=True), abs=1e-7)


@pytest.mark.parametrize("label", ["13C", "18O"])
def test_labeled_mass_only_mods_and_deltas_remain_additive(label):
    baseline = pt.parse(f"<{label}>PEPTIDE").mass(ion_type="a", deltas={"H-2O-1": 1})
    peptide = pt.parse(f"<{label}>PEP[+42]TIDE")
    assert peptide.mass(ion_type="a", deltas={"H-2O-1": 1, 17.0: 1}) == pytest.approx(baseline + 59)
    with pytest.raises(pt.CompositionError):
        peptide.comp()


@pytest.mark.parametrize(
    "sequence",
    [
        "PEPTIDE",
        "<13C>PEPTIDE",
        "<18O>PEPTIDE",
        "{Glycan:Hex}PEPTIDE",
        "PEP[Formula:Zn1:z+2]TIDE",
        "<[Formula:Zn1:z+2]@P>PEPTIDE",
    ],
)
@pytest.mark.parametrize("charge", [1, 2, -1])
@pytest.mark.parametrize("mono", [True, False])
def test_fast_fragment_agrees_at_every_position(sequence, charge, mono):
    peptide = pt.parse(sequence)
    result = peptide.fast_fragment(ion_types=["a", "b", "y", "p"], charges=[charge], monoisotopic=mono)
    for (ion, _), values in result.items():
        for position, value in enumerate(values, 1):
            expected = peptide.frag(ion_type=ion, charge=charge, monoisotopic=mono, position=None if ion == pt.IonType.PRECURSOR else position).mz
            assert value == pytest.approx(expected, abs=1e-7)


@pytest.mark.parametrize("charge", [0, True, 1.5])
def test_fast_fragment_rejects_invalid_charge(charge):
    with pytest.raises(ValueError, match="nonzero integers"):
        pt.parse("PEPTIDE").fast_fragment(charges=[charge])


def test_fast_fragment_rejects_unsupported_series():
    with pytest.raises(pt.UnsupportedOperationError, match="Use fragment"):
        pt.fast_fragment("PEPTIDE", ion_types=["by"])


@pytest.mark.parametrize("keep", [True, False])
@pytest.mark.parametrize("inplace", [True, False])
def test_filter_modes(keep, inplace):
    sequence = "[Acetyl]-PEM[Oxidation]TIDE"
    annotation = pt.parse(sequence)
    result = annotation.filter_mods("nterm", keep=keep, inplace=inplace)
    assert result.serialize() == ("[Acetyl]-PEMTIDE" if keep else "PEM[Oxidation]TIDE")
    assert (result is annotation) == inplace
    if not inplace:
        assert annotation.serialize() == sequence


def test_delta_cache_is_defensive_and_pickleable():
    source = {pt.ChargedFormula.from_string("H-2O-1", require_formula_prefix=False): 1}
    direct = DeltaInfo(source)
    source.clear()
    assert direct.deltas
    delta = DeltaInfo.from_input({"H-2O-1": 1})
    expected_mass = pt.mass("PEPTIDE", deltas={"H-2O-1": 1})
    expected_comp = pt.comp("PEPTIDE", deltas={"H-2O-1": 1})
    dict_view = delta.deltas
    dict_view.clear()
    assert pt.mass("PEPTIDE", deltas={"H-2O-1": 1}) == expected_mass
    assert pt.comp("PEPTIDE", deltas={"H-2O-1": 1}) == expected_comp
    restored = pickle.loads(pickle.dumps(delta))
    assert restored.composition == delta.composition
    assert restored.get_mass_delta() == delta.get_mass_delta()


def test_charge_cache_preserves_negative_atoms_and_returns_fresh_mapping():
    assert ChargeCarrierInfo.from_input(-2).composition == Counter({pt.ELEMENT_LOOKUP["H"]: -2})
    sodium = ChargeCarrierInfo.from_input("Na:z+1")
    sodium.to_proforma_charge.clear()
    assert sodium.to_proforma_charge == {"Na": 1}


@pytest.mark.parametrize(
    "kwargs",
    [
        {"isotopes": 1000},
        {"isotopes": {"13C": -1}},
        {"isotopes": {"13C": 1.5}},
        {"isotopes": {"C": 1}},
        {"isotopes": True},
        {"deltas": {"H-999": 1}},
        {"deltas": float("nan")},
        {"deltas": float("inf")},
        {"deltas": {"H2O": 1.5}},
    ],
)
@pytest.mark.parametrize("composition", [True, False])
def test_invalid_adjustments_rejected_in_both_paths(kwargs, composition):
    with pytest.raises(ValueError):
        pt.parse("PEPTIDE").mass(**kwargs, calculate_with_composition=composition)


def test_adjustments_include_precursor_terminal_atoms():
    peptide = pt.parse("G")
    for kwargs in [{"isotopes": {"18O": 2}}, {"deltas": {"H-4": 1}}]:
        assert peptide.mass(**kwargs) == pytest.approx(peptide.mass(**kwargs, calculate_with_composition=True))


@pytest.mark.parametrize("composition", [True, False])
def test_excessive_deprotonation_rejected(composition):
    with pytest.raises(pt.InvalidAdjustmentError, match="Negative element counts"):
        pt.parse("PEPTIDE").mass(charge=-100, calculate_with_composition=composition)


def test_integer_delta_direct_constructor_normalizes_mass():
    delta = DeltaInfo({42: 1})
    assert delta.has_floats
    with pytest.raises(pt.CompositionError):
        _ = delta.composition

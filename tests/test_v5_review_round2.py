"""PR #29 re-review items: parse input types, Interval.append_mod, EnzymeConfig.semi, the proton
carrier constant and single-value fragment options."""

import pytest
from tacular.constants import ELECTRON_MASS, HYDROGEN_MASS, PROTON_MASS

import peptacular as pt


class Entry:
    sequence = "PEM[Oxidation]TIDE"


def test_parse_accepts_objects_with_a_sequence():
    assert pt.parse(Entry()).serialize() == "PEM[Oxidation]TIDE"
    assert [a.serialize() for a in pt.parse([Entry(), "PEPTIDE"])] == ["PEM[Oxidation]TIDE", "PEPTIDE"]


@pytest.mark.parametrize("bad", [b"PEPTIDE", None, 5, [None], object()])
def test_parse_rejects_other_types_naming_the_accepted_ones(bad):
    with pytest.raises(TypeError, match="ProForma str, an object with a str 'sequence'"):
        pt.parse(bad)  # type: ignore[arg-type]


def test_interval_append_mod_returns_the_interval():
    interval = pt.Interval(1, 3)
    copy = interval.append_mod("Oxidation", inplace=False)
    assert isinstance(copy, pt.Interval) and copy is not interval
    assert copy.has_mods and not interval.has_mods
    assert interval.append_mod("Phospho") is interval and interval.has_mods


def test_enzyme_config_semi():
    config = pt.EnzymeConfig("trypsin", missed_cleavages=1, semi=True)
    assert config.semi and not hasattr(config, "semi_enzymatic")
    with pytest.raises(TypeError):
        pt.EnzymeConfig(enzyme="trypsin", semi_enzymatic=True)  # type: ignore[call-arg]
    annot = pt.parse("PEPTIDEKAAR")
    semi = {annot[sp].serialize() for sp in annot.sequential_digest_spans([pt.EnzymeConfig("trypsin", semi=True)])}
    full = {annot[sp].serialize() for sp in annot.sequential_digest_spans([pt.EnzymeConfig("trypsin")])}
    assert full < semi and "PEPTIDE" in semi


def test_proton_carrier_constant():
    assert pt.PROTON_CARRIER_MASS == HYDROGEN_MASS - ELECTRON_MASS
    assert pt.PROTON_CARRIER_MASS == pytest.approx(PROTON_MASS, abs=2e-8)
    annot = pt.parse("PEPTIDE")
    assert annot.mass(charge=2) - annot.mass() == pytest.approx(2 * pt.PROTON_CARRIER_MASS, abs=1e-9)
    b1, b2 = (annot.frag(ion_type="b", charge=z, position=3).mass for z in (1, 2))
    assert b2 - b1 == pytest.approx(pt.PROTON_CARRIER_MASS, abs=1e-9)
    fast = annot.fast_fragment(["b"], [2])[(pt.IonType.B, 2)][2]
    assert fast == pytest.approx(b2 / 2, abs=1e-9)


def test_single_fragment_options_need_no_list():
    assert pt.fragment("PEPTIDE", "b", 2) == pt.fragment("PEPTIDE", ["b"], [2])
    assert pt.parse("PEPTIDE").fast_fragment("y", 1) == pt.parse("PEPTIDE").fast_fragment(["y"], [1])
    # a multi-letter ion type is one ion type, not one per letter
    assert {f.ion_type for f in pt.fragment("PEPTIDE", "by", 1)} == {pt.IonType("by")}
    assert pt.fragment("PEPTIDE", "b", 1, neutral_deltas="H2O") == pt.fragment("PEPTIDE", ["b"], [1], neutral_deltas=["H2O"])
    assert pt.fragment("PEPTIDE", "b", 1, isotopes=1) == pt.fragment("PEPTIDE", ["b"], [1], isotopes=[1])


def _roundtrip_fragments():
    yield pt.parse("PEPTIDE").frag(ion_type="b", charge=2, position=3, deltas={"NH3": 1}, isotopes={"15N": 1})
    yield pt.fragment("PEPTIDE", ["b"], ["Na:z+1"])[3]
    yield pt.fragment("PEPTIDE", ["y"], [1])[2]


@pytest.mark.parametrize("field", ["isotopes", "deltas", "charge_adducts"])
def test_fragment_replace_accepts_its_own_property_values(field):
    # the properties return parsed objects (ChargedFormula keys, Mods); replace and the
    # constructor must take them back without breaking str/mzPAF or equality
    for frag in _roundtrip_fragments():
        copy = frag.replace(**{field: getattr(frag, field)})
        assert copy == frag
        assert copy.to_mzpaf() == frag.to_mzpaf() and str(copy) == str(frag)


def test_fragment_constructor_accepts_formula_delta_keys():
    frag = pt.parse("PEPTIDE").frag(ion_type="b", charge=1, position=3, deltas={"NH3": 1})
    rebuilt = pt.Fragment(frag.ion_type, frag.position, frag.mass, frag.monoisotopic, frag.charge_state, deltas=frag.deltas)
    assert rebuilt.to_mzpaf(include_sequence=False) == "b3-NH3"
    assert rebuilt.deltas == frag.deltas


def test_fragment_constructor_expands_a_mods_of_adducts():
    frag = pt.fragment("PEPTIDE", ["b"], ["Na:z+1"])[3]
    rebuilt = pt.Fragment(frag.ion_type, frag.position, frag.mass, frag.monoisotopic, frag.charge_state, charge_adducts=frag.charge_adducts)
    assert rebuilt._charge_adducts == ("Na:z+1",)


def test_fragment_replace_keeps_a_13c_count():
    frag = pt.fragment("PEPTIDE", ["b"], [1], isotopes=[1])[3]
    assert frag.is_c13
    copy = frag.replace(isotopes=frag.isotopes)
    assert copy.is_c13 and copy == frag
    assert pt.Fragment(frag.ion_type, frag.position, frag.mass, True, 1, isotopes={"15N": 1}).replace(isotopes={}).isotopes == {}

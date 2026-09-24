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

"""Regressions for bugs found by the reference-validation and property-test sweep."""

import pytest

import peptacular as pt
from peptacular.property.data import pk_cterminal

# --------------------------------------------------------------------------- pKa table


def test_cterminal_pka_glu_gln_match_source_table():
    # peptideweb pKa table (cited in property/data.py): Glu 2.19, Gln 2.17.
    assert pk_cterminal["E"] == 2.19
    assert pk_cterminal["Q"] == 2.17


# --------------------------------------------------------------------------- names ending in (...)

# Monoisotopic deltas from Unimod: Label:13C(6) = 6.020129, HexNAc(2) = 406.158745,
# Hex(1)HexNAc(1) = 365.132196.
PAREN_NAMES = [
    ("K", "U:Label:13C(6)", 6.020129),
    ("K", "Label:13C(6)", 6.020129),
    ("N", "HexNAc(2)", 406.158745),
    ("N", "U:HexNAc(2)", 406.158745),
    ("S", "Hex(1)HexNAc(1)", 365.132196),
]


@pytest.mark.parametrize(("aa", "tag", "delta"), PAREN_NAMES)
def test_trailing_parenthesis_is_part_of_the_name(aa, tag, delta):
    annot = pt.parse(f"PEPTIDE{aa}[{tag}]")
    assert annot.mass() - pt.parse(f"PEPTIDE{aa}").mass() == pytest.approx(delta, abs=1e-5)
    assert annot.serialize() == f"PEPTIDE{aa}[{tag}]"
    assert pt.parse(annot.serialize()) == annot


def test_localisation_score_still_follows_group_label():
    annot = pt.parse("EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK")
    assert annot.mass() == pytest.approx(pt.parse("EM[Oxidation]EVTSES[Phospho]PEK").mass(), abs=1e-9)
    assert pt.parse(annot.serialize()) == annot


def test_name_with_parenthesis_and_group_label():
    # ProForma 2.0 spec section 4.2.3: a PSI-MOD name containing "(...)" followed by a
    # cross-link group label. Mass is two half-cystines (-2 H in total).
    annot = pt.parse("EVTSEKC[L-cystine (cross-link)#XL1]LEMSC[#XL1]EFD")
    ref = pt.parse("EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD")
    assert annot.mass() == pytest.approx(ref.mass(), abs=1e-6)


# --------------------------------------------------------------------------- isotopes on labelled ions


def test_isotope_offset_on_fully_labelled_fragment():
    # y1 of a 13C6-lysine has no 12C left to swap for 13C: the M+1 peak must still be
    # reported, one 13C-12C mass difference above the monoisotopic peak.
    annot = pt.parse("PEPTIDEK[Formula:[13C6]C-6]")
    mono = {f.position: f.mz for f in annot.fragment(ion_types=["y"], isotopes=[0])}
    plus1 = {f.position: f.mz for f in annot.fragment(ion_types=["y"], isotopes=[1])}
    assert plus1.keys() == mono.keys()
    for pos, mz in mono.items():
        assert plus1[pos] - mz == pytest.approx(1.0033548378, abs=1e-6)


def test_isotope_offset_on_global_isotope_label():
    # With <13C> every carbon is already heavy; isotopes=1 is still the M+1 peak.
    annot = pt.parse("<13C>PEPTIDE")
    assert annot.mass(isotopes=1) - annot.mass() == pytest.approx(1.0033548378, abs=1e-6)
    assert pt.parse("K[Formula:[13C6]C-6]").mass(isotopes=2) - pt.parse("K[Formula:[13C6]C-6]").mass() == pytest.approx(2 * 1.0033548378, abs=1e-6)

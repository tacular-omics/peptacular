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


# --------------------------------------------------------------------------- cross-link separator


def test_parse_chimeric_rejects_crosslink_separator():
    # "//" joins cross-linked peptides into one ion; "+" joins separate (chimeric) ions.
    # Cross-links are unsupported, so reading "//" as "+" would silently change meaning.
    with pytest.raises(pt.UnsupportedOperationError, match="//"):
        pt.parse_chimeric("EMEVTK[XLMOD:02001#XL1]SESPEK//EMEVTK[#XL1]SESPEK")
    assert len(pt.parse_chimeric("PEPTIDE+PEPTIDE")) == 2


# --------------------------------------------------------------------------- satellite ions (mzPAF 1.0.1)

# mzPAF 1.0.1 section 4.4.3: d = Σn-1(AA) + offset, v/w = Σc-1(AA) + offset. The residue
# whose side chain is cleaved is not in the residue sum; the offset is its remnant.
_SAT_RESIDUE_TYPES = {
    "d": {"V": ["d-valine"], "I": ["da-isoleucine", "db-isoleucine"], "T": ["da-threonine", "db-threonine"], "G": [], "A": [], "P": []},
    "w": {"V": ["w-valine"], "I": ["wa-isoleucine", "wb-isoleucine"], "T": ["wa-threonine", "wb-threonine"], "G": [], "A": [], "P": []},
    "v": {},
}


def _residue_sum(seq: str) -> float:
    # Residue masses from pyteomics (independent of peptacular/tacular).
    from pyteomics.mass import std_aa_mass

    return sum(std_aa_mass[aa] for aa in seq)


def _expected_satellites(seq: str, family: str) -> dict[tuple[str, int], float]:
    from tacular import FRAGMENT_ION_LOOKUP

    proton = 1.007276466621
    n = len(seq)
    out = {}
    for i in range(1, n + 1):
        if family == "d":
            cleaved, kept = seq[i - 1], seq[: i - 1]
        else:
            cleaved, kept = seq[n - i], seq[n - i + 1 :]
        for ion in _SAT_RESIDUE_TYPES[family].get(cleaved, [family]):
            mz = _residue_sum(kept) + FRAGMENT_ION_LOOKUP[ion].monoisotopic_mass + proton
            out[(ion, i)] = mz
    return out


@pytest.mark.parametrize("seq", ["SAMPLER", "VTIVTI", "EDVKITLS"])
@pytest.mark.parametrize("family", ["d", "v", "w"])
def test_satellite_ions_use_residue_sum_without_cleaved_residue(seq, family):
    got = {(f.ion_type.value, f.position): f.mz for f in pt.fragment(seq, ion_types=[family], charges=[1])}
    expected = _expected_satellites(seq, family)
    assert got.keys() == expected.keys()
    for key, mz in expected.items():
        assert got[key] == pytest.approx(mz, abs=1e-6), key


# SAMPLVER: d5 cleaves L, d1 cleaves S, w3 (VER) cleaves V, v3 cleaves V, d6 cleaves V.
@pytest.mark.parametrize(
    ("ion", "position", "expected_type"),
    [("d", 5, "d"), ("d", 1, "d"), ("d", 6, "d-valine"), ("w", 3, "w-valine"), ("v", 3, "v")],
)
def test_frag_satellite_single_position(ion, position, expected_type):
    seq = "SAMPLVER"
    f = pt.parse(seq).frag(ion_type=ion, charge=1, position=position)
    assert f.ion_type.value == expected_type
    assert f.mz == pytest.approx(_expected_satellites(seq, ion[0])[(expected_type, position)], abs=1e-6)


def test_satellite_offsets_follow_mzpaf_when_tacular_has_them():
    from tacular import FRAGMENT_ION_LOOKUP

    if dict(FRAGMENT_ION_LOOKUP["d"].dict_composition or {}) != {"C": 2, "H": 4, "N": 1}:
        pytest.skip("installed tacular predates the mzPAF 1.0.1 d/v/w offsets")
    # d3 of SAMPLER (cleaves M): S + A residues + C2H4N + H+.
    f = pt.parse("SAMPLER").frag(ion_type="d", charge=1, position=3)
    assert f.mz == pytest.approx(87.032028 + 71.037114 + 42.034374 + 1.007276, abs=1e-5)


# --------------------------------------------------------------------------- mzPAF z / c variants


# mzPAF 1.0.1 "z" is z-dot (sum + H2O - NH2, pyteomics "z-dot"). Biemann z, z+H and c-H
# have no letter of their own, so they are written with a hydrogen delta.
@pytest.mark.parametrize(
    ("ion", "label", "delta_h"),
    [("z.", "z3{IDE}", 0), ("z", "z3{IDE}-H", -1), ("z+H", "z3{IDE}+H", 1), ("c-H", "c3{PEP}-H", -1)],
)
def test_mzpaf_z_and_c_variants(ion, label, delta_h):
    from pyteomics.mass import calculate_mass, nist_mass

    f = pt.parse("PEPIDE").frag(ion_type=ion, charge=1, position=3)
    assert f.to_mzpaf() == label
    series = "z-dot" if label[0] == "z" else "c"
    seq = "IDE" if label[0] == "z" else "PEP"
    # the label's meaning (series + H delta) is the fragment's own m/z
    assert f.mz == pytest.approx(calculate_mass(sequence=seq, ion_type=series, charge=1) + delta_h * nist_mass["H"][0][0], abs=1e-6)


def test_mzpaf_z_variant_delta_precedes_neutral_losses():
    f = pt.parse("PEPIDE").frag(ion_type="z", charge=2, position=3, deltas={"H2O": -1})
    assert f.to_mzpaf() == "z3{IDE}-H-H2O^2"

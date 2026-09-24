"""peptacular against independent references (``tests/reference/reference_data.json``).

The fixture is produced by ``tests/reference/generate_reference.py`` from pyteomics,
Biopython, ExPASy ProtScale, the peptideweb pKa table and the ProForma 2.0 specification;
see that file's header. Tolerance is 1e-6 Da unless a case says why it differs.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import peptacular as pt
from peptacular.property.data import PROPERTY_SCALES, pk_cterminal, pk_nterminal, pk_sidechain

REF = json.loads((Path(__file__).parent / "reference" / "reference_data.json").read_text())
TOL = 1e-6
# Average masses: pyteomics averages NIST abundances, peptacular uses IUPAC standard atomic
# weights; they differ in the fourth decimal for large peptides.
AVG_TOL = 1e-3
# Named mods are weighed from the CV's listed mass (6 decimals), not the formula, so each
# one can differ from the composition by up to 5e-7 Da.
MOD_TOL = 2e-6


# --------------------------------------------------------------------------- masses


@pytest.mark.parametrize("seq", sorted(REF["masses"]["unmodified"]))
def test_unmodified_mass_and_mz(seq):
    ref = REF["masses"]["unmodified"][seq]
    assert pt.mass(seq) == pytest.approx(ref["mono"], abs=TOL)
    assert pt.mass(seq, monoisotopic=False) == pytest.approx(ref["average"], abs=AVG_TOL)
    for z, mz in ref["mz"].items():
        assert pt.mz(seq, charge=int(z)) == pytest.approx(mz, abs=TOL)
        assert pt.mz(f"{seq}/{z}") == pytest.approx(mz, abs=TOL)


@pytest.mark.parametrize("proforma", sorted(REF["masses"]["modified"]))
def test_modified_mass_and_mz(proforma):
    ref = REF["masses"]["modified"][proforma]
    assert pt.mass(proforma) == pytest.approx(ref["mono"], abs=TOL)
    for z, mz in ref["mz"].items():
        assert pt.mz(proforma, charge=int(z)) == pytest.approx(mz, abs=TOL)


# --------------------------------------------------------------------------- fragments


@pytest.mark.parametrize("seq", sorted(REF["fragments"]["terminal"]))
def test_terminal_fragments(seq):
    for key, series in REF["fragments"]["terminal"][seq].items():
        ion, z = key.split("^")
        got = {f.position: f.mz for f in pt.fragment(seq, ion_types=[ion], charges=[int(z)])}
        assert sorted(got) == list(range(1, len(seq) + 1)), key
        for i, mz in enumerate(series, start=1):
            assert got[i] == pytest.approx(mz, abs=TOL), f"{key} {i}"


@pytest.mark.parametrize("seq", sorted(REF["fragments"]["terminal"]))
def test_fast_fragment_matches_reference(seq):
    for key, series in REF["fragments"]["terminal"][seq].items():
        ion, z = key.split("^")
        if ion not in "abcxyz":
            # fast_fragment covers the six main series and rejects the rest.
            with pytest.raises(pt.UnsupportedOperationError):
                pt.parse(seq).fast_fragment(ion_types=[ion], charges=[int(z)])
            continue
        fast = pt.parse(seq).fast_fragment(ion_types=[ion], charges=[int(z)])[(pt.IonType(ion), int(z))]
        assert list(fast) == pytest.approx(series, abs=TOL), key


@pytest.mark.parametrize("seq", sorted(REF["fragments"]["immonium"]))
def test_immonium(seq):
    got = {f.position: f.mz for f in pt.fragment(seq, ion_types=["i"], charges=[1])}
    for i, mz in enumerate(REF["fragments"]["immonium"][seq], start=1):
        assert got[i] == pytest.approx(mz, abs=TOL)


@pytest.mark.parametrize("seq", sorted(REF["fragments"]["internal"]))
def test_internal(seq):
    for ion, rows in REF["fragments"]["internal"][seq].items():
        got = {tuple(f.position): f.mz for f in pt.fragment(seq, ion_types=[ion], charges=[1])}
        for start, end, mz in rows:
            if (start, end) in got:
                assert got[(start, end)] == pytest.approx(mz, abs=TOL), f"{ion} {start}-{end}"
        # Internal ions exclude both termini, so a peptide of n residues has (n-1)(n-2)/2.
        assert len(got) == (len(seq) - 1) * (len(seq) - 2) // 2, ion


@pytest.mark.parametrize("seq", sorted(REF["fragments"]["neutral_loss"]))
def test_neutral_losses(seq):
    for key, series in REF["fragments"]["neutral_loss"][seq].items():
        ion, loss = key.split("-")
        lossy = [f for f in pt.fragment(seq, ion_types=[ion], charges=[1], neutral_deltas=[loss]) if f.deltas]
        sites = "STED" if loss == "H2O" else "RKNQ"
        assert bool(lossy) == any(aa in sites for aa in seq), key
        for f in lossy:
            assert f.mz == pytest.approx(series[f.position - 1], abs=TOL), f"{key} {f.position}"


@pytest.mark.parametrize("proforma", sorted(REF["fragments"]["modified"]))
def test_modified_fragments(proforma):
    for key, series in REF["fragments"]["modified"][proforma].items():
        ion, z = key.split("^")
        got = {f.position: f.mz for f in pt.fragment(proforma, ion_types=[ion], charges=[int(z)])}
        for i, mz in enumerate(series, start=1):
            assert got[i] == pytest.approx(mz, abs=MOD_TOL), f"{key} {i}"


# --------------------------------------------------------------------------- isotopes


@pytest.mark.parametrize("seq", sorted(REF["isotopes"]))
def test_isotope_distribution(seq):
    ref = REF["isotopes"][seq]
    got = pt.isotopic_distribution(seq, max_isotopes=len(ref), min_abundance_threshold=0.0)
    top = max(g.abundance for g in got)
    by_offset = {g.neutron_count: g for g in got}
    for offset, rel, mass in ref:
        g = by_offset[offset]
        assert g.abundance / top == pytest.approx(rel, abs=1e-6), offset
        assert g.mass == pytest.approx(mass, abs=TOL), offset


# --------------------------------------------------------------------------- digestion


@pytest.mark.parametrize("protease", sorted(REF["digestion"]["proteases"]))
def test_digestion(protease):
    protein = REF["digestion"]["protein"]
    for case, peptides in REF["digestion"]["proteases"][protease]["peptides"].items():
        mc = int(case[2])
        semi = case.endswith("_semi")
        got = {p for p, _ in pt.digest(protein, protease, missed_cleavages=mc, semi=semi)}
        assert got == set(peptides), case


# --------------------------------------------------------------------------- charge and pI


def test_pka_tables_match_source():
    for aa, (cterm, nterm, side) in REF["charge_pi"]["pka_table"].items():
        assert pk_cterminal[aa] == cterm
        assert pk_nterminal[aa] == nterm
        assert (pk_sidechain.get(aa) or None) == side


@pytest.mark.parametrize("seq", sorted(REF["charge_pi"]["peptides"]))
def test_charge_and_pi(seq):
    ref = REF["charge_pi"]["peptides"][seq]
    for ph, charge in ref["charge"].items():
        assert pt.charge_at_ph(seq, float(ph)) == pytest.approx(charge, abs=1e-9)
    # peptacular bisects [0, 14] until the bracket is under 0.001 pH units.
    assert pt.pi(seq) == pytest.approx(ref["pi"], abs=1e-3)


# --------------------------------------------------------------------------- property scales


@pytest.mark.parametrize("scale_id", sorted(REF["scales"]))
def test_property_scale_matches_published_table(scale_id):
    ref = REF["scales"][scale_id]
    table = PROPERTY_SCALES[scale_id]
    for aa, value in ref["values"].items():
        assert table[aa] == pytest.approx(value, abs=1e-9), f"{aa} ({ref['source']})"


# --------------------------------------------------------------------------- ProForma 2.0 spec examples

# Spec examples peptacular does not accept as one annotation, and why. Everything else
# in the spec must parse, round-trip and (where a reference exists) weigh the same.
NOT_PARSED = {
    # Cross-linked or branched peptidoforms joined by "//" (unsupported).
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK": pt.UnsupportedOperationError,
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[#XL1]SESPEK": pt.UnsupportedOperationError,
    "ETFGD[MOD:00093#BRANCH]LEVK//GGAGDSTK[#BRANCH]": pt.UnsupportedOperationError,
    "FVNQHLC[MOD:00034#XL1]GSHLVEALYLVC[MOD:00034#XL2]GERGFFYTPK//GIVEQC[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENYC[#XL2]N": (pt.UnsupportedOperationError),
    # Chimeric: parse with parse_chimeric (checked below).
    "EMEVEESPEK/2+ELVISLIVER/3": pt.UnsupportedOperationError,
    "EMEVEESPEK+ELVISLIVER": pt.UnsupportedOperationError,
    # Nested ranges and a group-labelled range: listed in the spec as not yet supported.
    "P(RT(ESFRMS)[+19.0523]IS)[+19.0523]K": pt.ProFormaFormatError,
    "PROT([#g1]EOC[Carbamidomethyl]FORMS)[+19.0523#g1]ISK": pt.ProFormaFormatError,
    # "(...)[Oxidation]^2": ^n is only defined for labile/unknown-position mods.
    "MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxidation]^2[half cystine]^2": pt.ProFormaFormatError,
    # Ion adduct forms in "/z[...]" notation.
    "EMEVEESPEK/2[+2Na+,+H+]": pt.ProFormaFormatError,
    "EMEVEESPEK/1[+2Na+,-H+]": pt.ProFormaFormatError,
    "EMEVEESPEK/-2[2I-]": pt.ProFormaFormatError,
    "EMEVEESPEK/-1[+e-]": pt.ProFormaFormatError,
}
# The spec gives these as incorrect usage: a CV prefix followed by an accession number.
SPEC_INVALID = {
    "EM[M:00719]EVEES[M:00046]PEK": pt.UnknownModificationError,
    "EM[U:35]EVEES[U:56]PEK": pt.UnknownModificationError,
    "EM[R:AA0581]EVEES[R:AA0037]PEK": pt.UnknownModificationError,
    "EM[R: Methionine sulfone]EVEES[O-phospho-L-serine]PEK": pt.UnknownModificationError,
}
# Parse, but a mod has no mass in the bundled vocabularies (or the spec text has a
# stray space from the PDF: "R: L-methionine sulfone").
NO_MASS = {
    "EM[R: L-methionine sulfone]EVEES[O-phospho-L-serine]PEK",
    "EM[RESID:AA0581]EVEES[RESID:AA0037]PEK",
    "LEIK[N6-(L-asparagyl)-L-lysine#XL1]KIPHDN[#XL1]",
    "{Glycan:Hex10HexNAc4}YPVLN[MOD:00006]VTMPN[MOD:00006]NSNGKFDK",
    "SEQUEN[Lipid:OleicAcid]CE",
    "SEQUEN[Formula:CH3(CH2)4CH3]CE",
}
# Mass tolerances wider than 1e-6 Da, and why.
MASS_TOL = {
    # GNO glycans: pyteomics weighs the GlyTouCan composition with its own monosaccharide
    # table, which differs from the formula masses by up to 7e-5 Da per residue.
    "GNO:": 1e-3,
    "G:G59626AS": 1e-3,
    # half cystine (MOD:00798) is weighed from PSI-MOD's listed DiffMono, rounded to 1e-6.
    "MOD:00798": 5e-6,
    "half cystine": 5e-6,
}

SPEC = REF["proforma_spec"]


def _tol(proforma: str) -> float:
    return max([t for key, t in MASS_TOL.items() if key in proforma], default=TOL)


@pytest.mark.parametrize("row", SPEC, ids=[r["proforma"][:50] for r in SPEC])
def test_proforma_spec_example(row):
    s = row["proforma"]
    if s in NOT_PARSED or s in SPEC_INVALID:
        expected = NOT_PARSED.get(s) or SPEC_INVALID[s]
        with pytest.raises(expected):
            pt.parse(s).mass()
        return
    annot = pt.parse(s)
    assert pt.parse(annot.serialize()) == annot
    assert annot.serialize() == pt.parse(annot.serialize()).serialize()
    if s in NO_MASS:
        return
    if row["reference_mass"] is None:
        pytest.fail(f"no reference mass: {row['pyteomics_error']}")
    assert annot.mass(charge=0) == pytest.approx(row["reference_mass"], abs=_tol(s)), row["source"]


def test_spec_chimeric_examples():
    for s in ("EMEVEESPEK/2+ELVISLIVER/3", "EMEVEESPEK+ELVISLIVER"):
        parts = list(pt.parse_chimeric(s))
        assert [p.sequence for p in parts] == ["EMEVEESPEK", "ELVISLIVER"]
        assert pt.serialize_chimeric(parts) == s


def test_spec_lenient_unknown_position_after_nterm():
    # The spec lists "[Acetyl]-[Phospho]^2?EM..." as incorrect (unknown-position mods
    # must come first). peptacular accepts it; recorded here so a change is deliberate.
    annot = pt.parse("[Acetyl]-[Phospho]^2?EM[Oxidation]EVTSESPEK")
    assert annot.mass() == pytest.approx(pt.parse("[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK").mass(), abs=TOL)

"""A bare X (unknown residue, no mass-bearing modification) has no mass and must raise.

``X[+100]`` and other mass-bearing modifications give the residue its mass (ProForma), so those
keep working. Digestion is not a mass calculation and is unaffected.
"""

import pytest

import peptacular as pt
from peptacular.diagnostics import CompositionError, PeptacularError

BARE_X = ["PEPXTIDE", "XPEPTIDE", "PEPTIDEX", "X", "PEPX[INFO:unknown]TIDE", "PEPX[+100]TIDEX"]
ERROR = "Mass not available for amino acid: X"


@pytest.mark.parametrize("seq", BARE_X)
@pytest.mark.parametrize("monoisotopic", [True, False])
def test_bare_x_mass_raises(seq: str, monoisotopic: bool) -> None:
    with pytest.raises(PeptacularError, match=ERROR):
        pt.mass(seq, monoisotopic=monoisotopic)


@pytest.mark.parametrize("seq", BARE_X)
def test_bare_x_mass_with_composition_raises(seq: str) -> None:
    with pytest.raises(CompositionError, match=ERROR):
        pt.mass(seq, calculate_with_composition=True)


@pytest.mark.parametrize("seq", BARE_X)
def test_bare_x_composition_raises(seq: str) -> None:
    with pytest.raises(CompositionError, match=ERROR):
        pt.comp(seq)


@pytest.mark.parametrize("seq", BARE_X)
def test_bare_x_fragment_raises(seq: str) -> None:
    with pytest.raises(PeptacularError, match=ERROR):
        pt.fragment(seq)
    with pytest.raises(PeptacularError, match=ERROR):
        pt.fragment(seq, calculate_with_composition=True)


@pytest.mark.parametrize("seq", BARE_X)
def test_bare_x_fast_fragment_raises(seq: str) -> None:
    with pytest.raises(PeptacularError, match=ERROR):
        pt.parse(seq).fast_fragment()


@pytest.mark.parametrize("seq", BARE_X)
def test_bare_x_isotopic_distribution_raises(seq: str) -> None:
    with pytest.raises(PeptacularError, match=ERROR):
        pt.parse(seq).isotopic_distribution()


def test_bare_x_mz_raises() -> None:
    with pytest.raises(PeptacularError, match=ERROR):
        pt.mz("PEPXTIDE", charge=2)


def test_x_with_mass_mod() -> None:
    assert pt.mass("PEPX[+100]TIDE") == pytest.approx(pt.mass("PEPTIDE") + 100.0)
    assert pt.mass("PEPX[+100]TIDE", monoisotopic=False) == pytest.approx(pt.mass("PEPTIDE", monoisotopic=False) + 100.0)
    frags = pt.fragment("PEPX[+100]TIDE")
    assert frags
    assert pt.parse("PEPX[+100]TIDE").fast_fragment()


def test_x_with_formula_mod_equals_glycine() -> None:
    # C2H3NO is the glycine residue, so X[Formula:C2H3NO] weighs the same as G.
    assert pt.mass("PEPX[Formula:C2H3NO]TIDE") == pytest.approx(pt.mass("PEPGTIDE"))
    assert pt.comp("PEPX[Formula:C2H3NO]TIDE") == pt.comp("PEPGTIDE")
    assert pt.mass("PEPX[Formula:C2H3NO]TIDE", calculate_with_composition=True) == pytest.approx(pt.mass("PEPGTIDE"))
    assert [f.mz for f in pt.fragment("PEPX[Formula:C2H3NO]TIDE")] == pytest.approx([f.mz for f in pt.fragment("PEPGTIDE")])


def test_x_with_named_mod() -> None:
    assert pt.mass("PEPX[Phospho]TIDE") == pytest.approx(pt.mass("PEPT[Phospho]IDE"))


def test_x_with_static_mass_mod() -> None:
    assert pt.mass("<[+5]@X>PEPXTIDE") == pytest.approx(pt.mass("PEPTIDE") + 5.0)


def test_condensed_ambiguity_still_has_mass() -> None:
    annot = pt.parse("PEP(?TIDE)[Phospho]")
    assert annot.condense_ambiguity_to_xnotation().mass() == pytest.approx(annot.mass(), abs=1e-4)


def test_digest_x_sequence_unchanged() -> None:
    peptides = [peptide for peptide, _ in pt.digest("PEPXKTIDXRAA", "trypsin")]
    assert peptides == ["PEPXK", "TIDXR", "AA"]


def test_bare_x_fragment_arrays_raises() -> None:
    pytest.importorskip("numpy")
    with pytest.raises(PeptacularError, match=ERROR):
        pt.fragment_arrays(["PEPXTIDE"])


@pytest.mark.parametrize(
    ("seq", "expected"),
    [("PEP(X)[+10]K", 479.2536), ("(XX)[+10]K", 156.1055)],
)
def test_x_in_interval_with_mass_mod(seq: str, expected: float) -> None:
    assert pt.parse(seq).mass() == pytest.approx(expected, abs=1e-3)


def test_x_in_interval_with_info_mod_raises() -> None:
    with pytest.raises(PeptacularError, match=ERROR):
        pt.parse("PEP(X)[INFO:unknown]K").mass()

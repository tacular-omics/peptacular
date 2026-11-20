"""A bare X (unknown residue, no modification) contributes zero mass and no atoms.

ProForma 2.1 Final, section 7.1 (Ambiguous amino acids):

    "X: Any amino acid (see also section 7.3 (Specifying a gap of known mass), for the use
    of X). We note that the character X itself is assigned zero mass in this notation, and
    it SHOULD be followed by a modification to provide the residue's total mass."

Section 7.3 (Specifying a gap of known mass) says the same thing from the other side: the
gap notation works "meaning that the mass of X is actually zero".

Note that "SHOULD be followed by a modification" is a recommendation, not a requirement, so a
bare X is a legal sequence whose residue weighs nothing. B and Z are different: section 7.1
gives them two possible masses each and asks implementations to expand them combinatorially,
so peptacular still refuses to weigh those.

The absolute expected values in this file come from pyteomics 5.0.1, computed outside
peptacular with X added to ``pyteomics.mass.std_aa_comp`` as an empty composition, which is
what the specification's zero mass means in atoms:

    >>> comp = dict(pyteomics.mass.std_aa_comp); comp["X"] = pyteomics.mass.Composition({})
    >>> pyteomics.mass.calculate_mass(sequence="PEPXTIDE", aa_comp=comp)   # 799.3599634
    >>> pyteomics.mass.calculate_mass(sequence="RTAAXWT", aa_comp=comp)    # 704.3605728
"""

import pytest

import peptacular as pt
from peptacular.diagnostics import CompositionError, PeptacularError

BARE_X = ["PEPXTIDE", "XPEPTIDE", "PEPTIDEX", "PEPX[INFO:unknown]TIDE"]

# pyteomics 5.0.1 with X as an empty composition; see the module docstring.
PYTEOMICS_MONO = {
    "PEPTIDE": 799.359964,
    "PEPXTIDE": 799.359964,
    "RTAAXWT": 704.360573,
    "XPEPTIDEX": 799.359964,
    "PEPXXTIDE": 799.359964,
}
PYTEOMICS_AVERAGE = {"PEPTIDE": 799.823877, "PEPXTIDE": 799.823877}


@pytest.mark.parametrize("seq", sorted(PYTEOMICS_MONO))
def test_monoisotopic_mass_matches_pyteomics(seq: str) -> None:
    """The expected value is pyteomics 5.0.1's, not ours. Fails on 5.0.1 through 5.2.0."""
    assert pt.mass(seq) == pytest.approx(PYTEOMICS_MONO[seq], abs=1e-5)


@pytest.mark.parametrize("seq", sorted(PYTEOMICS_AVERAGE))
def test_average_mass_matches_pyteomics(seq: str) -> None:
    assert pt.mass(seq, monoisotopic=False) == pytest.approx(PYTEOMICS_AVERAGE[seq], abs=1e-4)


def test_composition_matches_pyteomics() -> None:
    # pyteomics.mass.Composition(sequence="PEPXTIDE", aa_comp=comp) is
    # {'H': 53, 'C': 34, 'O': 15, 'N': 7}, the same as for PEPTIDE: X adds no atoms.
    assert pt.chem_formula(pt.comp("PEPXTIDE")) == "C34H53N7O15"
    assert pt.comp("PEPXTIDE") == pt.comp("PEPTIDE")


def test_the_specifications_own_gap_example() -> None:
    """``RTAAX[+367.0537]WT`` is printed in section 7.3; the 367.0537 is the spec's number."""
    # 704.360573 (pyteomics, X = 0) + 367.0537 (printed in the specification) = 1071.414273.
    assert pt.mass("RTAAX[+367.0537]WT") == pytest.approx(1071.414273, abs=1e-5)
    # The other example in 7.3, a sequence tag whose termini are gaps of known mass.
    assert pt.mass("X[+379.147]UEX[+365.113]") == pytest.approx(pt.mass("UE") + 379.147 + 365.113, abs=1e-5)


def test_a_bare_x_weighs_nothing_by_hand() -> None:
    """Computed by hand from the residue masses, independently of the pyteomics numbers."""
    # A peptide is the sum of its residues plus one water. X adds a residue position that
    # contributes nothing, so PEPXTIDE weighs exactly what PEPTIDE weighs, and an X on its own
    # weighs one water.
    assert pt.mass("X") == pytest.approx(pt.chem_mass("H2O"), abs=1e-9)
    assert pt.mass("XXXXX") == pytest.approx(pt.chem_mass("H2O"), abs=1e-9)
    assert pt.mass("PEPXTIDE") == pytest.approx(pt.mass("PEPTIDE"), abs=1e-9)
    assert pt.comp("XXXXX") == pt.chem_comp("H2O")


@pytest.mark.parametrize("seq", BARE_X)
def test_every_mass_entry_point_accepts_a_bare_x(seq: str) -> None:
    """mass, mz, comp, fragment, fast_fragment and isotopic_distribution, i.e. what #46 broke."""
    assert pt.mass(seq) == pytest.approx(799.359964, abs=1e-5)
    assert pt.mz(seq + "/2") == pytest.approx(pt.mz("PEPTIDE/2"), abs=1e-9)
    assert pt.comp(seq) == pt.comp("PEPTIDE")
    # The ion masses are PEPTIDE's, plus the duplicates the zero-mass position adds.
    assert set(round(f.mz, 6) for f in pt.fragment("PEPTIDE")) <= set(round(f.mz, 6) for f in pt.fragment(seq))
    assert len(pt.fast_fragment(seq)) == len(pt.fast_fragment("PEPTIDE"))
    assert [round(d.mass, 6) for d in pt.isotopic_distribution(seq, max_isotopes=3)] == [
        round(d.mass, 6) for d in pt.isotopic_distribution("PEPTIDE", max_isotopes=3)
    ]


def test_fragment_arrays_accepts_a_bare_x() -> None:
    pytest.importorskip("numpy")
    import numpy as np

    arrays = pt.fragment_arrays(["PEPXTIDE"])
    plain = pt.fragment_arrays(["PEPTIDE"])
    # One more ion than PEPTIDE, because the X is a cleavage position of its own, and the ion
    # that ends on the X has the same mass as the one before it.
    assert len(arrays["mz"]) == len(plain["mz"]) + 2
    assert np.isfinite(arrays["mz"]).all()


def test_the_x_adds_a_fragment_position_of_equal_mass() -> None:
    # b3 is PEP and b4 is PEPX, so the two weigh the same.
    by_position = {f.position: f.mz for f in pt.fragment("PEPXTIDE", ion_types="b", charges=1)}
    assert by_position[4] == pytest.approx(by_position[3], abs=1e-9)


@pytest.mark.parametrize("seq", ["PEPX[+100]TIDE", "PEPX[Formula:C2H3NO]TIDE", "PEPX[Phospho]TIDE"])
def test_a_modified_x_is_unchanged(seq: str) -> None:
    """A modification on the X still provides the residue's mass, as section 7.1 recommends."""
    assert pt.mass(seq) > pt.mass("PEPTIDE")


def test_x_with_mass_mod() -> None:
    assert pt.mass("PEPX[+100]TIDE") == pytest.approx(pt.mass("PEPTIDE") + 100.0)
    assert pt.mass("PEPX[+100]TIDE", monoisotopic=False) == pytest.approx(pt.mass("PEPTIDE", monoisotopic=False) + 100.0)


def test_x_with_formula_mod_equals_glycine() -> None:
    # Glycine's residue composition is C2H3NO, so X carrying that formula weighs what G weighs.
    assert pt.mass("PEPX[Formula:C2H3NO]TIDE") == pytest.approx(pt.mass("PEPGTIDE"))
    assert pt.comp("PEPX[Formula:C2H3NO]TIDE") == pt.comp("PEPGTIDE")
    assert pt.mass("PEPX[Formula:C2H3NO]TIDE", calculate_with_composition=True) == pytest.approx(pt.mass("PEPGTIDE"))
    assert [f.mz for f in pt.fragment("PEPX[Formula:C2H3NO]TIDE")] == pytest.approx([f.mz for f in pt.fragment("PEPGTIDE")])


def test_x_with_named_mod() -> None:
    assert pt.mass("PEPX[Phospho]TIDE") == pytest.approx(pt.mass("PEPT[Phospho]IDE"))


def test_x_with_static_mass_mod() -> None:
    assert pt.mass("<[+5]@X>PEPXTIDE") == pytest.approx(pt.mass("PEPTIDE") + 5.0)


def test_condensed_ambiguity_still_has_mass() -> None:
    annot = pt.parse("PEP(?MT)IDE")
    assert annot.condense_ambiguity_to_xnotation().mass() == pytest.approx(annot.mass(), abs=1e-4)


def test_digest_x_sequence_unchanged() -> None:
    assert [s for s, _ in pt.digest("PEPXTIDEK", enzyme="trypsin")] == ["PEPXTIDEK"]


@pytest.mark.parametrize(("seq", "expected"), [("PEP(X)[+10]K", 479.2536), ("(XX)[+10]K", 156.1055)])
def test_x_in_interval_with_mass_mod(seq: str, expected: float) -> None:
    assert pt.parse(seq).mass() == pytest.approx(expected, abs=1e-3)


def test_x_in_interval_with_info_mod_does_not_raise() -> None:
    # An INFO modification carries no mass, so the X inside the interval simply weighs nothing.
    assert pt.parse("PEP(X)[INFO:unknown]K").mass() == pytest.approx(pt.mass("PEPK"), abs=1e-9)


@pytest.mark.parametrize("residue", ["B", "Z"])
@pytest.mark.parametrize("monoisotopic", [True, False])
def test_b_and_z_still_refuse_to_be_weighed(residue: str, monoisotopic: bool) -> None:
    """Section 7.1 gives B and Z two possible masses each, so neither has one mass to report."""
    with pytest.raises(PeptacularError, match=f"Mass not available for amino acid: {residue}"):
        pt.mass(f"PEP{residue}TIDE", monoisotopic=monoisotopic)
    with pytest.raises(CompositionError, match=f"Composition not available for amino acid: {residue}"):
        pt.comp(f"PEP{residue}TIDE")

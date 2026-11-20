"""A bare X (unknown residue, no mass-bearing modification) has no mass and must raise.

This is a deliberate departure from ProForma 2.1 Final, section 7.1 (Ambiguous amino acids):

    "X: Any amino acid (see also section 7.3 (Specifying a gap of known mass), for the use
    of X). We note that the character X itself is assigned zero mass in this notation, and
    it SHOULD be followed by a modification to provide the residue's total mass."

peptacular follows the recommendation in that sentence rather than the zero in front of it: a
sequence weighed with a zero X comes out about one residue light with nothing to say so, so
every mass path refuses a bare X and the error names the two ways to write the mass out.
``X[+100]`` gives the residue a mass and ``X[+0]`` asks for zero mass explicitly, which is the
notation's own value said out loud. Digestion is not a mass calculation and is unaffected.

pyteomics 5.0.1, run outside peptacular on 2026-11-20, refuses the same input and supplies the
absolute expected values in this file:

    >>> pyteomics.mass.fast_mass("PEPXTIDE")              # PyteomicsError: No mass data for residue: X
    >>> pyteomics.mass.calculate_mass(sequence="PEPXTIDE")  # PyteomicsError: Unknown label: X
    >>> pyteomics.mass.calculate_mass(sequence="PEPTIDE")   # 799.35996402671
    >>> pyteomics.mass.calculate_mass(sequence="PEPTIDE", average=True)  # 799.8238768361077
    >>> pyteomics.mass.Composition(sequence="PEPTIDE")    # {'H': 53, 'C': 34, 'O': 15, 'N': 7}
    >>> pyteomics.mass.calculate_mass(sequence="RTAAWT")  # 704.3605731634
"""

import pytest

import peptacular as pt
from peptacular.diagnostics import CompositionError, PeptacularError

BARE_X = ["PEPXTIDE", "XPEPTIDE", "PEPTIDEX", "X", "PEPX[INFO:unknown]TIDE", "PEPX[+100]TIDEX"]
ERROR = "X has no mass of its own"

# pyteomics 5.0.1 values; see the module docstring. X[+0] is a zero-mass residue, so these are
# the values of the same sequence without the X.
PYTEOMICS_PEPTIDE_MONO = 799.359964
PYTEOMICS_PEPTIDE_AVERAGE = 799.823877
PYTEOMICS_RTAAWT_MONO = 704.360573


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


def test_bare_x_fragment_arrays_raises() -> None:
    pytest.importorskip("numpy")
    with pytest.raises(PeptacularError, match=ERROR):
        pt.fragment_arrays(["PEPXTIDE"])


def test_the_a_ion_series_of_a_leading_bare_x_raises_on_both_paths() -> None:
    """An a ion asks for the b ion minus CO, which is where a zero-mass X used to go negative.

    Both mass paths refuse the sequence itself, so neither reports an impossible ion and neither
    drops one silently.
    """
    with pytest.raises(PeptacularError, match=ERROR):
        pt.fragment("XPEPTIDE", ion_types=["a"])
    with pytest.raises(PeptacularError, match=ERROR):
        pt.fragment("XPEPTIDE", ion_types=["a"], calculate_with_composition=True)
    # The a ion series of a modified X is defined, so the refusal is about the bare X only.
    assert len(pt.fragment("X[+100]PEPTIDE", ion_types=["a"])) == 8


def test_the_error_says_what_to_do() -> None:
    """The message names both spellings, because the fix is to write the mass out."""
    with pytest.raises(PeptacularError) as exc:
        pt.mass("PEPXTIDE")
    assert str(exc.value) == "X has no mass of its own. Attach a modification, for example X[+100], or X[+0] for zero mass."


def test_x_with_zero_mass_mod_is_the_notations_value_written_out() -> None:
    """``X[+0]`` is the zero of section 7.1, asked for explicitly. pyteomics supplies the values."""
    assert pt.mass("PEPX[+0]TIDE") == pytest.approx(PYTEOMICS_PEPTIDE_MONO, abs=1e-5)
    assert pt.mass("PEPX[+0]TIDE", monoisotopic=False) == pytest.approx(PYTEOMICS_PEPTIDE_AVERAGE, abs=1e-4)
    assert pt.mass("PEPX[+0]TIDE", calculate_with_composition=True) == pytest.approx(PYTEOMICS_PEPTIDE_MONO, abs=1e-5)
    # pyteomics.mass.Composition(sequence="PEPTIDE") is {'H': 53, 'C': 34, 'O': 15, 'N': 7}.
    assert pt.chem_formula(pt.comp("PEPX[+0]TIDE")) == "C34H53N7O15"


def test_every_mass_entry_point_accepts_a_zero_mass_x() -> None:
    """mass, mz, comp, fragment, fast_fragment, fragment_arrays and isotopic_distribution."""
    pytest.importorskip("numpy")
    seq = "PEPX[+0]TIDE"
    assert pt.mass(seq) == pytest.approx(PYTEOMICS_PEPTIDE_MONO, abs=1e-5)
    assert pt.mz(seq + "/2") == pytest.approx(pt.mz("PEPTIDE/2"), abs=1e-9)
    assert pt.comp(seq) == pt.comp("PEPTIDE")
    # The X is a cleavage position of its own, so the ion that ends on it repeats the mass of the
    # ion before it and PEPTIDE's ions are a subset of this sequence's.
    assert set(round(f.mz, 6) for f in pt.fragment("PEPTIDE")) <= set(round(f.mz, 6) for f in pt.fragment(seq))
    assert len(pt.fast_fragment(seq)) == len(pt.fast_fragment("PEPTIDE"))
    assert [round(d.mass, 6) for d in pt.isotopic_distribution(seq, max_isotopes=3)] == [
        round(d.mass, 6) for d in pt.isotopic_distribution("PEPTIDE", max_isotopes=3)
    ]
    assert len(pt.fragment_arrays([seq])["mz"]) == len(pt.fragment_arrays(["PEPTIDE"])["mz"]) + 2


def test_the_specifications_own_gap_example() -> None:
    """``RTAAX[+367.0537]WT`` is printed in section 7.3; the 367.0537 is the specification's."""
    # 704.360573 (pyteomics RTAAWT) + 367.0537 (printed in section 7.3) = 1071.414273.
    assert pt.mass("RTAAX[+367.0537]WT") == pytest.approx(PYTEOMICS_RTAAWT_MONO + 367.0537, abs=1e-5)
    # The other example in 7.3, a sequence tag whose termini are gaps of known mass.
    assert pt.mass("X[+379.147]UEX[+365.113]") == pytest.approx(pt.mass("UE") + 379.147 + 365.113, abs=1e-5)


def test_x_with_mass_mod() -> None:
    assert pt.mass("PEPX[+100]TIDE") == pytest.approx(PYTEOMICS_PEPTIDE_MONO + 100.0, abs=1e-5)
    assert pt.mass("PEPX[+100]TIDE", monoisotopic=False) == pytest.approx(PYTEOMICS_PEPTIDE_AVERAGE + 100.0, abs=1e-4)
    assert pt.fragment("PEPX[+100]TIDE")
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
    assert pt.mass("<[+5]@X>PEPXTIDE") == pytest.approx(PYTEOMICS_PEPTIDE_MONO + 5.0, abs=1e-5)


def test_unknown_position_modification_does_not_cover_a_bare_x() -> None:
    """A modification of unknown position belongs to no residue, so it leaves the X bare."""
    with pytest.raises(PeptacularError, match=ERROR):
        pt.mass("[+100]?PEPXTIDE")


def test_condensed_ambiguity_still_has_mass() -> None:
    annot = pt.parse("PEP(?TIDE)[Phospho]")
    assert annot.condense_ambiguity_to_xnotation().mass() == pytest.approx(annot.mass(), abs=1e-4)


def test_digest_x_sequence_unchanged() -> None:
    peptides = [peptide for peptide, _ in pt.digest("PEPXKTIDXRAA", "trypsin")]
    assert peptides == ["PEPXK", "TIDXR", "AA"]


@pytest.mark.parametrize(("seq", "expected"), [("PEP(X)[+10]K", 479.2536), ("(XX)[+10]K", 156.1055)])
def test_x_in_interval_with_mass_mod(seq: str, expected: float) -> None:
    assert pt.parse(seq).mass() == pytest.approx(expected, abs=1e-3)


def test_x_in_interval_with_info_mod_raises() -> None:
    with pytest.raises(PeptacularError, match=ERROR):
        pt.parse("PEP(X)[INFO:unknown]K").mass()


@pytest.mark.parametrize("residue", ["B", "Z"])
@pytest.mark.parametrize("monoisotopic", [True, False])
def test_b_and_z_refuse_for_their_own_reason(residue: str, monoisotopic: bool) -> None:
    """Section 7.1 gives B and Z two possible masses each, so neither has one mass to report."""
    with pytest.raises(PeptacularError, match=f"Mass not available for amino acid: {residue}"):
        pt.mass(f"PEP{residue}TIDE", monoisotopic=monoisotopic)
    with pytest.raises(CompositionError, match=f"Composition not available for amino acid: {residue}"):
        pt.comp(f"PEP{residue}TIDE")

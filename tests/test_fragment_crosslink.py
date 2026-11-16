"""fragment() refuses a cross-linked or branched peptidoform instead of answering as if it were linear.

A label on two residues (``#XL1``, ``#BRANCH``) bonds those residues to each other, so one
backbone cleavage between them does not separate the two halves. peptacular does not model that
bond (ProForma 2.1 section 9.2, PROFORMA_COMPLIANCE.md), and the interchain form (``//``) is
already refused when parsed.
"""

import math

import pytest

import peptacular as pt

# Cross-linked and branched peptidoforms: the label is on two residues.
BONDED = [
    "EVTK[X:DSS#XL1]LEK[#XL1]SEFD",
    "EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD",
    "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR",
    "ETFGD[MOD:00093#BRANCH]LEVKGGAGDSTK[#BRANCH]",
]

# Not bonded: one dead-end linker (ProForma 2.1 section 9.2.1), an ordinary modification, and a
# localization group, which names one modification on one of several residues rather than a bond.
UNBONDED = [
    "EMEVTK[X:DSS#XL1]SESPEK",
    "EVTKLEK[XLMOD:02001]SEFD",
    "PEPT[Phospho#g1]IDEK[#g1]",
]


class TestFragmentRefusesCrossLinks:
    @pytest.mark.parametrize("seq", BONDED)
    @pytest.mark.parametrize("ion_type", ["b", "y", "a", "c", "x", "z", "p"])
    def test_fragment_raises_for_every_ion_type(self, seq: str, ion_type: str):
        with pytest.raises(pt.UnsupportedOperationError, match="cross-linked or branched"):
            pt.fragment(seq, ion_types=(ion_type,), charges=(1,))

    @pytest.mark.parametrize("seq", BONDED)
    def test_every_fragment_entry_point_raises(self, seq: str):
        annot = pt.parse(seq)
        with pytest.raises(pt.UnsupportedOperationError, match="cross-linked or branched"):
            annot.fragment(ion_types=("b", "y"), charges=(1,))
        with pytest.raises(pt.UnsupportedOperationError, match="cross-linked or branched"):
            annot.fast_fragment(ion_types=("b", "y"), charges=(1,))
        with pytest.raises(pt.UnsupportedOperationError, match="cross-linked or branched"):
            annot.frag("b", 1, position=4)
        with pytest.raises(pt.UnsupportedOperationError, match="cross-linked or branched"):
            pt.fast_fragment(seq, ion_types=("b",), charges=(1,))

    @pytest.mark.parametrize("seq", BONDED)
    def test_fragment_arrays_raises(self, seq: str):
        pytest.importorskip("numpy")
        with pytest.raises(pt.UnsupportedOperationError, match="cross-linked or branched"):
            pt.fragment_arrays([seq])

    def test_error_names_both_labels_of_a_two_crosslink_peptide(self):
        with pytest.raises(pt.UnsupportedOperationError, match=r"#XL1, #XL2"):
            pt.fragment("EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR", ion_types=("b",), charges=(1,))

    def test_the_impossible_ions_are_no_longer_returned(self):
        # b4, b5 and b6 of this peptide used to come back as 596.329, 709.413 and 838.456: the
        # cleavage sites lie between K4 and K7, which the cross-linker joins.
        with pytest.raises(pt.UnsupportedOperationError):
            pt.fragment("EVTK[X:DSS#XL1]LEK[#XL1]SEFD", ion_types=("b",), charges=(1,))

    def test_diagnose_reports_the_refusal(self):
        diagnostic = pt.diagnose("EVTK[X:DSS#XL1]LEK[#XL1]SEFD", "fragment")
        assert diagnostic is not None
        assert diagnostic.code == "unsupported_operation"
        assert "cross-linked or branched" in diagnostic.message

    @pytest.mark.parametrize("seq", BONDED)
    def test_the_intact_molecule_still_weighs(self, seq: str):
        annot = pt.parse(seq)
        assert math.isfinite(annot.mass())
        assert math.isfinite(annot.neutral_mass())
        assert math.isfinite(annot.mz(charge=2))


class TestFragmentKeepsWorkingWithoutABond:
    @pytest.mark.parametrize("seq", UNBONDED)
    def test_fragment_still_returns_ions(self, seq: str):
        fragments = pt.fragment(seq, ion_types=("b", "y"), charges=(1,))
        assert fragments
        assert all(math.isfinite(fragment.mz) for fragment in fragments)

    def test_dead_end_linker_weighs_as_an_ordinary_modification(self):
        # One #XL1 label and no partner: a hydrolyzed linker, bonded to nothing.
        dead_end = pt.fragment("EMEVTK[X:DSS#XL1]SESPEK", ion_types=("b",), charges=(1,))
        plain = pt.fragment("EMEVTK[X:DSS]SESPEK", ion_types=("b",), charges=(1,))
        assert [fragment.mz for fragment in dead_end] == [pytest.approx(fragment.mz) for fragment in plain]

    @pytest.mark.parametrize("seq", UNBONDED)
    def test_fast_fragment_still_returns_ions(self, seq: str):
        assert pt.fast_fragment(seq, ion_types=("b", "y"), charges=(1,))

"""Charge state / charge carrier notation — ProForma 2.1 section 11.5 compliance."""

import pytest

import peptacular as pt


class TestSpecChargeExamples:
    """Every literal example from ProForma 2.1 section 11.5 must parse, round-trip, and be mass-computable."""

    SPEC_EXAMPLES = [
        "PEPTIDE/2",
        "PEPTIDE/[Na:z+1]",
        "PEPTIDE/[Na:z+1,H:z+1]",
        "PEPTIDE/[[15N1]H4:z+1]",
        "PEPTIDE/[Na:z+1^2]",
        "PEPT[Formula:Zn:z+2]IDE/2",
        "PEPT[Formula:Zn:z+2]IDE/[Na:z+1^2]",
    ]

    @pytest.mark.parametrize("seq", SPEC_EXAMPLES)
    def test_round_trip(self, seq):
        assert pt.parse(seq).serialize() == seq

    @pytest.mark.parametrize("seq", SPEC_EXAMPLES)
    def test_mass_computable(self, seq):
        assert pt.parse(seq).mass() > 0
        assert pt.parse(seq).mz() > 0


class TestBareChargeCarriers:
    """A charge carrier is a bare charged formula and its mass must be correct."""

    def test_sodium_adduct_mass(self):
        # [M+Na]+ = neutral + Na - electron
        base = pt.mass("PEPTIDE")
        na = pt.parse("PEPTIDE/[Na:z+1]")
        # mz for singly charged == the charged mass
        assert na.mz() == pytest.approx(base + 22.98977 - 0.000549, abs=1e-3)

    def test_occurrence_specifier_doubles_charge(self):
        assert pt.parse("PEPTIDE/[H:z+1^2]").charge_state == 2

    def test_bare_formula_carrier(self):
        assert pt.parse("PEPTIDE/[C2H6:z+2]").mz() > 0


class TestPrefixedCarrierRejected:
    """A CV/type prefix ('Formula:', 'Glycan:') is not valid in a charge carrier."""

    @pytest.mark.parametrize("seq", ["PEPTIDE/[Formula:C2H6:z+2]", "PEPTIDE/[Glycan:HexNAc:z+1]"])
    def test_prefixed_carrier_rejected_on_mass(self, seq):
        with pytest.raises(ValueError, match="Invalid charge carrier"):
            pt.parse(seq).mass()

    def test_prefixed_carrier_rejected_on_validate(self):
        with pytest.raises(ValueError, match="Invalid charge carrier"):
            pt.parse("PEPTIDE/[Formula:C2H6:z+2]", validate=True)

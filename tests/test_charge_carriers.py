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


class TestSetChargeInputs:
    """Regression tests for set_charge input handling (pre-release sweep)."""

    def test_mod_wrapped_carrier_matches_bare_carrier(self):
        # A Mod[GlobalChargeCarrier] must serialize its wrapped carrier, not the dataclass
        # repr; the result must equal passing the bare GlobalChargeCarrier.
        from peptacular.annotation.mod import Mod
        from peptacular.proforma_components.comps import GlobalChargeCarrier

        gcc = GlobalChargeCarrier.charged_proton(2)
        bare = pt.parse("PEPTIDE").set_charge(gcc, inplace=False)
        wrapped = pt.parse("PEPTIDE").set_charge(Mod(gcc, 1), inplace=False)
        assert wrapped._charge == bare._charge == ["H:z+1^2"]
        assert wrapped.serialize() == bare.serialize()
        assert wrapped.mass() == pytest.approx(bare.mass())

    def test_charge_zero_clears_to_none(self):
        # A charge of 0 is neutral -> no charge component (equal to an unset peptide).
        a = pt.parse("PEPTIDE")
        assert a.set_charge(0, inplace=False)._charge is None
        assert a.set_charge(0, inplace=False) == a

    def test_bool_charge_rejected(self):
        with pytest.raises(ValueError, match="Unsupported charge type"):
            pt.parse("PEPTIDE").set_charge(True, inplace=False)

    def test_duplicate_adduct_list_not_collapsed(self):
        # Two identical adducts must not be deduped into one: charge_state, serialization
        # and mass must all reflect both carriers (regression: {str: 1} dict collapsed them).
        a = pt.parse("PEPTIDE").set_charge(["Na:z+1", "Na:z+1"], inplace=False)
        assert a.charge_state == 2
        assert a.serialize() == "PEPTIDE/[Na:z+1,Na:z+1]"
        # equivalent to the occurrence-specifier form and mass-consistent with it
        assert pt.parse(a.serialize()).charge_state == 2
        assert a.mass() == pytest.approx(pt.parse("PEPTIDE/[Na:z+1^2]").mass())

    def test_mods_input_preserves_occurrence_count(self):
        # A Mods carrying a carrier with count 2 must keep both carriers, not drop the count.
        from peptacular.annotation.mod import Mods
        from peptacular.constants import ModType

        a = pt.parse("PEPTIDE").set_charge(Mods(mod_type=ModType.CHARGE, _mods={"Na:z+1": 2}), inplace=False)
        assert a._charge == ["Na:z+1", "Na:z+1"]
        assert a.charge_state == 2

    def test_mod_count_zero_clears_to_none(self):
        # A Mod-wrapped carrier with count 0 is neutral: clear to None, never serialize 'PEPTIDE/[]'.
        from peptacular.annotation.mod import Mod
        from peptacular.proforma_components.comps import GlobalChargeCarrier

        a = pt.parse("PEPTIDE").set_charge(Mod(GlobalChargeCarrier.charged_proton(2), 0), inplace=False)
        assert a._charge is None
        assert a.serialize() == "PEPTIDE"


class TestChargeCarrierMzPaf:
    """mzPAF serialization of a charge carrier must render the sign correctly."""

    @pytest.mark.parametrize(
        "charge,expected",
        [(1, "M+H"), (2, "M+2H"), (-1, "M-H"), (-2, "M-2H")],
    )
    def test_proton_carrier_sign(self, charge, expected):
        # A negative charge (negative occurance) previously produced malformed 'M+-2H'.
        from peptacular.proforma_components.comps import GlobalChargeCarrier

        assert GlobalChargeCarrier.charged_proton(charge).to_mz_paf() == expected

    def test_deprotonation_formula_sign(self):
        from peptacular.proforma_components.comps import GlobalChargeCarrier

        assert GlobalChargeCarrier.from_string("H-1:z-1^2").to_mz_paf() == "M-2H"

    def test_negative_occurrence_roundtrips_internally(self):
        # peptacular represents a -1 charge proton carrier as 'H:z+1^-1'; it must re-parse.
        from peptacular.proforma_components.comps import GlobalChargeCarrier

        assert GlobalChargeCarrier.from_string("H:z+1^-1").occurance == -1

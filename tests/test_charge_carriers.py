"""Charge state / charge carrier notation — ProForma 2.1 section 11.5 compliance."""

import pytest

import peptacular as pt
from peptacular.diagnostics import InvalidAdjustmentError


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


class TestChimericChargeCarriers:
    """parse_chimeric must accept the same charge notation as parse (regression)."""

    @pytest.mark.parametrize("chains", [["PEPTIDE", "ELVIS/[Na:z+1]"], ["PEPTIDE/[H:z+1^2]", "ELVIS"], ["PEPTIDE/[Na:z+1,H:z+1]", "ELVIS/[Cl:z-1]"]])
    def test_carriers_in_chimeric(self, chains):
        seq = "+".join(chains)
        parts = list(pt.parse_chimeric(seq))
        assert pt.serialize_chimeric(parts) == seq
        for part, single in zip(parts, chains, strict=True):
            assert part == pt.parse(single)
            assert part.mass() == pytest.approx(pt.parse(single).mass())


ELECTRON = 0.000548579909065
CARRIER_GRID = ["H:z+1", "Na:z+1", "K:z+1", "[15N1]H4:z+1", "H-1:z-1", "Cl:z-1"]


class TestNegativeCarrierSerialization:
    """A negative carrier occurrence is written as the negated carrier so the string re-parses (regression).

    ``charged_proton(-2)`` is stored as ``H:z+1^-2``, which ProForma rejects; it is
    serialized as the equal-mass ``H-1:z-1^2``. A bare ``/-2`` is left as written.
    """

    @pytest.mark.parametrize("z", [-3, -2, -1, 1, 2, 3])
    def test_set_charge_from_fragment_adducts_reparses(self, z):
        frag = pt.parse("PEPTIDEK").frag("p", z)
        annot = pt.parse("PEPTIDEK").set_charge(frag.charge_adducts)
        text = annot.serialize()
        reparsed = pt.parse(text)
        assert reparsed.serialize() == text
        assert reparsed.mass() == pytest.approx(pt.parse(f"PEPTIDEK/{z}").mass(), abs=1e-9)

    @pytest.mark.parametrize("carrier", CARRIER_GRID)
    @pytest.mark.parametrize("occ", [-3, -2, -1, 1, 2, 3])
    def test_carrier_grid_parse_serialize(self, carrier, occ):
        from peptacular.proforma_components.comps import GlobalChargeCarrier

        gcc = GlobalChargeCarrier.from_string(f"{carrier}^{occ}")
        text = f"PEPTIDEK/[{gcc}]"
        annot = pt.parse(text)
        assert annot.serialize() == text
        assert GlobalChargeCarrier.from_string(str(gcc)).to_mz_paf() == gcc.to_mz_paf()
        assert GlobalChargeCarrier.from_string(str(gcc)).get_mass() == pytest.approx(gcc.get_mass(), abs=1e-9)
        try:
            delta = annot.mass() - pt.parse("PEPTIDEK").mass()
        except InvalidAdjustmentError:
            # Removing atoms the peptide lacks (e.g. K-1:z-1^3) is a typed error, not a mass.
            assert any(fe.occurance * occ < 0 for fe in gcc.charged_formula.formula)
        else:
            # occurrence x (carrier formula mass - its charge x electron mass)
            expected = gcc.get_mass() - gcc.occurance * gcc.charged_formula.charge * ELECTRON
            assert delta == pytest.approx(expected, abs=1e-6)

    def test_bare_negative_charge_unchanged(self):
        assert pt.parse("PEPTIDEK/-2").serialize() == "PEPTIDEK/-2"

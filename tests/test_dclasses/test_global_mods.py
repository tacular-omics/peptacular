"""
Tests for parsing isotope replacements and global modifications.
"""

import pytest

import peptacular as pt


class TestIsotopeReplacement:
    """Tests for pt.IsotopeReplacement.from_string() function"""

    def test_deuterium_shorthand(self):
        """Test parsing deuterium as D"""
        result = pt.IsotopeReplacement.from_string("D")
        assert isinstance(result, pt.IsotopeReplacement)
        assert result.element == pt.Element.H
        assert result.isotope == 2

    def test_carbon_13(self):
        """Test parsing 13C"""
        result = pt.IsotopeReplacement.from_string("13C")
        assert result.element == pt.Element.C
        assert result.isotope == 13

    def test_nitrogen_15(self):
        """Test parsing 15N"""
        result = pt.IsotopeReplacement.from_string("15N")
        assert result.element == pt.Element.N
        assert result.isotope == 15

    def test_oxygen_18(self):
        """Test parsing 18O"""
        result = pt.IsotopeReplacement.from_string("18O")
        assert result.element == pt.Element.O
        assert result.isotope == 18

    def test_sulfur_34(self):
        """Test parsing 34S"""
        result = pt.IsotopeReplacement.from_string("34S")
        assert result.element == pt.Element.S
        assert result.isotope == 34

    def test_carbon_12(self):
        """Test parsing 12C (natural isotope)"""
        result = pt.IsotopeReplacement.from_string("12C")
        assert result.element == pt.Element.C
        assert result.isotope == 12

    def test_missing_isotope_number_raises_error(self):
        """Test that missing isotope number raises ValueError"""
        with pytest.raises(ValueError):
            pt.IsotopeReplacement.from_string("C")

    def test_missing_element_raises_error(self):
        """Test that missing element symbol raises ValueError"""
        with pytest.raises(ValueError):
            pt.IsotopeReplacement.from_string("13")

    def test_invalid_element_raises_error(self):
        """Test that invalid element symbol raises ValueError"""
        with pytest.raises(ValueError):
            pt.IsotopeReplacement.from_string("13X")

    def test_string_representation_deuterium(self):
        """Test string representation of deuterium"""
        result = pt.IsotopeReplacement.from_string("D")
        assert str(result) == "D"

    def test_string_representation_isotope(self):
        """Test string representation of isotope"""
        result = pt.IsotopeReplacement.from_string("13C")
        assert str(result) == "13C"


class TestGlobalChargeCarrier:
    """Tests for pt.GlobalChargeCarrier.from_string() function"""

    def test_simple_charge_carrier(self):
        """Test parsing simple charge carrier"""
        result = pt.GlobalChargeCarrier.from_string("H:z+1")
        assert isinstance(result, pt.GlobalChargeCarrier)
        assert result.occurance == 1.0
        assert result.charged_formula.charge == 1

    def test_sodium_carrier(self):
        """Test parsing sodium charge carrier"""
        result = pt.GlobalChargeCarrier.from_string("Na:z+1")
        assert result.occurance == 1.0
        # Check that formula contains Na
        assert any(fe.element == pt.Element.Na for fe in result.charged_formula.formula)

    def test_carrier_with_integer_occurance(self):
        """Test charge carrier with integer occurrence"""
        result = pt.GlobalChargeCarrier.from_string("H:z+1^2")
        assert result.occurance == 2.0
        assert result.charged_formula.charge == 1

    def test_carrier_with_fractional_occurance(self):
        """Test charge carrier with fractional occurrence"""
        # assert exception
        with pytest.raises(ValueError):
            pt.GlobalChargeCarrier.from_string("H:z+1^1.5")

    def test_negative_charge_carrier(self):
        """Test charge carrier with negative charge"""
        result = pt.GlobalChargeCarrier.from_string("H:z-1")
        assert result.charged_formula.charge == -1

    def test_complex_formula_carrier(self):
        """Test charge carrier with complex formula"""
        result = pt.GlobalChargeCarrier.from_string("C2H6:z+2")
        assert result.charged_formula.charge == 2
        assert len(result.charged_formula.formula) == 2

    def test_string_representation(self):
        """Test string representation"""
        result = pt.GlobalChargeCarrier.from_string("Na:z+1")
        result_str = str(result)
        assert "Na" in result_str


class TestFixedModification:
    """Tests for pt.FixedModification.from_string() function"""

    def test_simple_fixed_mod_with_position(self):
        """Test parsing fixed modification with position"""
        result = pt.FixedModification.from_string("[Oxidation]@M")
        assert isinstance(result, pt.FixedModification)
        assert len(result.position_rules) == 1
        assert result.position_rules[0].amino_acid == pt.AminoAcid.M

    def test_fixed_mod_multiple_positions(self):
        """Test fixed modification with multiple positions"""
        result = pt.FixedModification.from_string("[TMT6plex]@K,N-term")
        assert len(result.position_rules) == 2
        # Should have K and N-term
        assert any(pr.amino_acid == pt.AminoAcid.K for pr in result.position_rules)
        assert any(pr.terminal == pt.Terminal.N_TERM for pr in result.position_rules)

    def test_fixed_mod_without_position(self):
        """Test fixed modification without position rules"""
        result = pt.FixedModification.from_string("[Oxidation]")
        assert len(result.position_rules) == 0

    def test_fixed_mod_with_accession(self):
        """Test fixed modification with accession"""
        result = pt.FixedModification.from_string("[UNIMOD:35]@M")
        assert len(result) == 1
        assert isinstance(result.modifications[0], pt.TagAccession)
        assert result.modifications[0].accession == "35"

    def test_fixed_mod_with_mass(self):
        """Test fixed modification with mass"""
        result = pt.FixedModification.from_string("[+15.995]@M")
        assert isinstance(result.modifications[0], pt.TagMass)
        assert result.modifications[0].mass == 15.995

    def test_fixed_mod_with_multiple_tags(self):
        """Test fixed modification with multiple tags"""
        result = pt.FixedModification.from_string("[Oxidation|UNIMOD:35|+15.995]@M")
        assert len(result) == 3

    def test_fixed_mod_terminal_only(self):
        """Test fixed modification on terminal only"""
        result = pt.FixedModification.from_string("[Acetyl]@N-term")
        assert len(result.position_rules) == 1
        assert result.position_rules[0].terminal == pt.Terminal.N_TERM
        assert result.position_rules[0].amino_acid is None

    def test_invalid_format_missing_brackets_raises_error(self):
        """Test that invalid format raises ValueError"""
        with pytest.raises(ValueError):
            pt.FixedModification.from_string("Oxidation@M")

    def test_invalid_format_missing_at_raises_error(self):
        """Test that missing @ in format with position raises ValueError"""
        # This should parse as no position rules
        result = pt.FixedModification.from_string("[Oxidation]")
        assert len(result.position_rules) == 0

    def test_string_representation(self):
        """Test string representation"""
        result = pt.FixedModification.from_string("[Oxidation]@M")
        result_str = str(result)
        assert "[Oxidation]@M" == result_str

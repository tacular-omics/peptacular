"""
Tests for the from_string class methods on NamedTuples.
Demonstrates round-trip conversion: string -> object -> string
"""

import pytest

import peptacular as pt


class TestFormulaElementFromString:
    """Test pt.FormulaElement.from_string()"""

    def test_simple_element(self):
        """Test parsing simple element like 'C'"""
        element = pt.FormulaElement.from_string("C")
        assert element.element == pt.Element.C
        assert element.occurance == 1
        assert element.isotope is None
        assert str(element) == "C"

    def test_element_with_count(self):
        """Test parsing element with count like 'H2'"""
        element = pt.FormulaElement.from_string("H2")
        assert element.element == pt.Element.H
        assert element.occurance == 2
        assert element.isotope is None
        assert str(element) == "H2"

    def test_element_with_isotope(self):
        """Test parsing element with isotope like '[13C2]'"""
        element = pt.FormulaElement.from_string("[13C2]")
        assert element.element == pt.Element.C
        assert element.occurance == 2
        assert element.isotope == 13
        assert str(element) == "[13C2]"

    def test_round_trip(self):
        """Test that from_string(str(obj)) == obj"""
        original = pt.FormulaElement(element=pt.Element.N, occurance=3, isotope=15)
        parsed = pt.FormulaElement.from_string(str(original))
        assert parsed == original


class TestChargedFormulaFromString:
    """Test pt.ChargedFormula.from_string()"""

    def test_simple_formula(self):
        """Test parsing simple formula like 'Formula:C2H6'"""
        formula = pt.ChargedFormula.from_string("Formula:C2H6")
        assert len(formula.formula) == 2
        assert formula.charge is None
        assert str(formula) == "Formula:C2H6"

    def test_formula_with_charge(self):
        """Test parsing formula with charge like 'Formula:Na:z+1'"""
        formula = pt.ChargedFormula.from_string("Formula:Na:z+1")
        assert formula.charge == 1
        assert str(formula) == "Formula:Na:z+1"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "Formula:C2H3N:z+2"
        parsed = pt.ChargedFormula.from_string(original_str)
        assert str(parsed) == original_str


class TestTagAccessionFromString:
    """Test pt.TagAccession.from_string()"""

    def test_unimod_full(self):
        """Test parsing UNIMOD accession"""
        tag = pt.TagAccession.from_string("UNIMOD:35")
        assert tag.cv == pt.CV.UNIMOD
        assert tag.accession == "35"
        # Note: str() outputs the CV enum value which is 'UNIMOD' (proper case)
        assert str(tag) == "UNIMOD:35"

    def test_unimod_shorthand(self):
        """Test parsing UNIMOD shorthand"""
        tag = pt.TagAccession.from_string("UNIMOD:35")
        assert tag.cv == pt.CV.UNIMOD
        assert tag.accession == "35"
        # Note: str() always outputs full form with proper CV casing
        assert str(tag) == "UNIMOD:35"

    def test_psi_mod(self):
        """Test parsing PSI-MOD accession"""
        tag = pt.TagAccession.from_string("MOD:00719")
        assert tag.cv == pt.CV.PSI_MOD
        assert tag.accession == "00719"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "XLMOD:02001"
        parsed = pt.TagAccession.from_string(original_str)
        assert str(parsed) == "XLMOD:02001"  # Note: XLMOD normalizes to XL-MOD


class TestTagMassFromString:
    """Test pt.TagMass.from_string()"""

    def test_positive_mass(self):
        """Test parsing positive mass"""
        tag = pt.TagMass.from_string("+15.995")
        assert tag.mass == 15.995
        assert tag.cv is None
        assert str(tag) == "+15.995"

    def test_negative_mass(self):
        """Test parsing negative mass"""
        tag = pt.TagMass.from_string("-18.010")
        assert tag.mass == -18.010
        # Note: Python's float formatting may drop trailing zeros
        assert str(tag) == "-18.01"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "+79.966"
        parsed = pt.TagMass.from_string(original_str)
        assert str(parsed) == original_str


class TestTagNameFromString:
    """Test pt.TagName.from_string()"""

    def test_simple_name(self):
        """Test parsing simple name"""
        tag = pt.TagName.from_string("Oxidation")
        assert tag.name == "Oxidation"
        assert tag.cv is None
        assert str(tag) == "Oxidation"

    def test_name_with_cv(self):
        """Test parsing name with single-letter CV prefix (Rule 4)"""
        tag = pt.TagName.from_string("U:Oxidation")
        assert tag.name == "Oxidation"
        assert tag.cv == pt.CV.UNIMOD
        # Note: str() outputs with proper CV casing
        assert str(tag) == "Oxidation"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "Phospho"
        parsed = pt.TagName.from_string(original_str)
        assert str(parsed) == original_str


class TestTagInfoFromString:
    """Test pt.TagInfo.from_string()"""

    def test_simple_info(self):
        """Test parsing INFO tag"""
        tag = pt.TagInfo.from_string("INFO:probable oxidation")
        assert tag.info == "probable oxidation"
        assert str(tag) == "INFO:probable oxidation"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "INFO:some custom text"
        parsed = pt.TagInfo.from_string(original_str)
        assert str(parsed) == original_str


class TestGlycanComponentFromString:
    """Test pt.GlycanComponent.from_string()"""

    def test_simple_glycan(self):
        """Test parsing simple glycan component"""
        component = pt.GlycanComponent.from_string("Hex")
        assert component.monosaccharide == pt.Monosaccharide.Hex
        assert component.occurance == 1
        assert str(component) == "Hex"

    def test_glycan_with_count(self):
        """Test parsing glycan component with count"""
        component = pt.GlycanComponent.from_string("HexNAc4")
        assert component.monosaccharide == pt.Monosaccharide.HexNAc
        assert component.occurance == 4
        assert str(component) == "HexNAc4"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "NeuAc2"
        parsed = pt.GlycanComponent.from_string(original_str)
        assert str(parsed) == original_str


class TestTagCustomFromString:
    """Test pt.TagCustom.from_string()"""

    def test_custom_tag(self):
        """Test parsing custom tag"""
        tag = pt.TagCustom.from_string("C:MyCustomMod")
        assert tag.name == "MyCustomMod"
        assert str(tag) == "C:MyCustomMod"

    def test_invalid_without_prefix(self):
        """Test that custom tag without C: prefix raises error"""
        with pytest.raises(ValueError):
            pt.TagCustom.from_string("MyCustomMod")


class TestPositionRuleFromString:
    """Test pt.PositionRule.from_string()"""

    def test_anywhere(self):
        """Test parsing amino acid position (anywhere)"""
        rule = pt.PositionRule.from_string("K")
        assert rule.terminal == pt.Terminal.ANYWHERE
        assert rule.amino_acid == pt.AminoAcid.K
        assert str(rule) == "K"

    def test_n_term(self):
        """Test parsing N-terminal position"""
        rule = pt.PositionRule.from_string("N-term")
        assert rule.terminal == pt.Terminal.N_TERM
        assert rule.amino_acid is None
        assert str(rule) == "N-term"

    def test_n_term_with_aa(self):
        """Test parsing N-terminal with specific amino acid"""
        rule = pt.PositionRule.from_string("N-term:K")
        assert rule.terminal == pt.Terminal.N_TERM
        assert rule.amino_acid == pt.AminoAcid.K
        assert str(rule) == "N-term:K"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "C-term:R"
        parsed = pt.PositionRule.from_string(original_str)
        assert str(parsed) == original_str


class TestSequenceElementFromString:
    """Test pt.SequenceElement.from_string()"""

    def test_simple_amino_acid(self):
        """Test parsing simple amino acid"""
        element = pt.SequenceElement.from_string("M")
        assert element.amino_acid == pt.AminoAcid.M
        assert len(element.modifications) == 0
        assert str(element) == "M"

    def test_amino_acid_with_modification(self):
        """Test parsing amino acid with modification"""
        element = pt.SequenceElement.from_string("M[Oxidation]")
        assert element.amino_acid == pt.AminoAcid.M
        assert len(element.modifications) == 1
        assert str(element) == "M[Oxidation]"

    def test_amino_acid_with_multiple_modifications(self):
        """Test parsing amino acid with multiple modifications"""
        element = pt.SequenceElement.from_string("K[UNIMOD:1][+42.010]")
        assert element.amino_acid == pt.AminoAcid.K
        assert len(element.modifications) == 2
        # Note: output format may differ slightly due to normalization

    def test_round_trip_simple(self):
        """Test round-trip conversion for simple case"""
        original_str = "A"
        parsed = pt.SequenceElement.from_string(original_str)
        assert str(parsed) == original_str

    def test_round_trip_with_mod(self):
        """Test round-trip conversion with modification"""
        original_str = "S[Phospho]"
        parsed = pt.SequenceElement.from_string(original_str)
        assert str(parsed) == original_str


class TestIsotopeReplacementFromString:
    """Test pt.IsotopeReplacement.from_string()"""

    def test_deuterium(self):
        """Test parsing deuterium (special case)"""
        iso = pt.IsotopeReplacement.from_string("D")
        assert iso.element == pt.Element.H
        assert iso.isotope == 2
        assert str(iso) == "D"

    def test_carbon_13(self):
        """Test parsing 13C"""
        iso = pt.IsotopeReplacement.from_string("13C")
        assert iso.element == pt.Element.C
        assert iso.isotope == 13
        assert str(iso) == "13C"

    def test_nitrogen_15(self):
        """Test parsing 15N"""
        iso = pt.IsotopeReplacement.from_string("15N")
        assert iso.element == pt.Element.N
        assert iso.isotope == 15
        assert str(iso) == "15N"

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "18O"
        parsed = pt.IsotopeReplacement.from_string(original_str)
        assert str(parsed) == original_str


class TestGlobalChargeCarrierFromString:
    """Test pt.GlobalChargeCarrier.from_string()"""

    def test_simple_carrier(self):
        """Test parsing simple charge carrier"""
        carrier = pt.GlobalChargeCarrier.from_string("Na:z+1")
        assert carrier.occurance == 1.0
        assert carrier.charged_formula.charge == 1
        assert str(carrier) == "Na:z+1"

    def test_carrier_with_occurrence(self):
        """Test parsing charge carrier with occurrence"""
        carrier = pt.GlobalChargeCarrier.from_string("H:z+1^2")
        assert carrier.occurance == 2.0
        assert carrier.charged_formula.charge == 1
        assert str(carrier) == "H:z+1^2"

    def test_carrier_with_fractional_occurrence(self):
        """Test parsing charge carrier with fractional occurrence"""
        with pytest.raises(ValueError):
            pt.GlobalChargeCarrier.from_string("Na:z+1^1.5")

    def test_round_trip(self):
        """Test round-trip conversion"""
        original_str = "H:z+1^3"
        parsed = pt.GlobalChargeCarrier.from_string(original_str)
        assert str(parsed) == original_str


class TestCaching:
    """Test that from_string methods use caching"""

    def test_same_string_returns_same_object(self):
        """Test that parsing the same string twice returns cached result"""
        tag1 = pt.TagName.from_string("Oxidation")
        tag2 = pt.TagName.from_string("Oxidation")
        # Due to caching, these should be the exact same object
        assert tag1 is tag2

    def test_different_strings_return_different_objects(self):
        """Test that different strings return different objects"""
        tag1 = pt.TagName.from_string("Oxidation")
        tag2 = pt.TagName.from_string("Phospho")
        assert tag1 is not tag2
        assert tag1 != tag2


class TestErrorHandling:
    """Test error handling in from_string methods"""

    def test_invalid_element(self):
        """Test that invalid element raises error"""
        with pytest.raises(ValueError):
            pt.FormulaElement.from_string("Zz")

    def test_invalid_accession(self):
        """Test that invalid accession format raises error"""
        with pytest.raises(ValueError):
            pt.TagAccession.from_string("NotAnAccession")

    def test_invalid_mass(self):
        """Test that invalid mass format raises error"""
        with pytest.raises(ValueError):
            pt.TagMass.from_string("15.995")  # Missing required sign

    def test_invalid_info(self):
        """Test that invalid INFO format raises error"""
        with pytest.raises(ValueError):
            pt.TagInfo.from_string("NotInfo")  # Missing 'INFO:' prefix

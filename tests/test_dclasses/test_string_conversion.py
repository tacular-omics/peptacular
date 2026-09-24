"""
Tests for converting dataclass objects back to ProForma strings.
"""

import tacular

import peptacular as pt


class TestTagStringConversion:
    """Tests for converting tag objects to strings"""

    def test_tag_name_to_string(self):
        """Test TagName string conversion"""
        tag = pt.TagName(name="Oxidation")
        assert str(tag) == "Oxidation"

    def test_tag_name_with_cv_to_string(self):
        """Test TagName with CV string conversion"""
        tag = pt.TagName(name="Oxidation", cv=pt.CV.UNIMOD)
        assert str(tag) == "Oxidation"

    def test_tag_accession_to_string(self):
        """Test TagAccession string conversion"""
        tag = pt.TagAccession(accession="35", cv=pt.CV.UNIMOD)
        assert str(tag) == "UNIMOD:35"

    def test_tag_accession_psi_mod(self):
        """Test TagAccession PSI-MOD string conversion"""
        tag = pt.TagAccession(accession="00719", cv=pt.CV.PSI_MOD)
        assert str(tag) == "MOD:00719"

    def test_tag_mass_positive(self):
        """Test TagMass positive value string conversion"""
        tag = pt.TagMass(mass_str="15.995")
        assert str(tag) == "+15.995"

    def test_tag_mass_negative(self):
        """Test TagMass negative value string conversion"""
        tag = pt.TagMass(mass_str="-18.010")
        assert str(tag) == "-18.01"

    def test_tag_mass_with_cv(self):
        """Test TagMass with CV string conversion"""
        tag = pt.TagMass(mass_str="15.995", cv=pt.CV.UNIMOD)
        assert str(tag) == "U:+15.995"

    def test_tag_info_to_string(self):
        """Test TagInfo string conversion"""
        tag = pt.TagInfo(info="custom annotation")
        assert str(tag) == "INFO:custom annotation"


class TestFormulaStringConversion:
    """Tests for converting formula objects to strings"""

    def test_formula_element_simple(self):
        """Test FormulaElement string conversion"""
        elem = pt.FormulaElement(element=tacular.Element.C, occurance=2)
        assert str(elem) == "C2"

    def test_formula_element_single_count(self):
        """Test FormulaElement with count of 1"""
        elem = pt.FormulaElement(element=tacular.Element.O, occurance=1)
        assert str(elem) == "O"

    def test_formula_element_with_isotope(self):
        """Test FormulaElement with isotope"""
        elem = pt.FormulaElement(element=tacular.Element.C, occurance=2, isotope=13)
        assert str(elem) == "[13C2]"

    def test_charged_formula_simple(self):
        """Test ChargedFormula string conversion"""
        formula = pt.ChargedFormula(
            formula=(
                pt.FormulaElement(element=tacular.Element.C, occurance=2),
                pt.FormulaElement(element=tacular.Element.H, occurance=6),
            )
        )
        assert str(formula) == "Formula:C2H6"

    def test_charged_formula_with_charge(self):
        """Test ChargedFormula with charge"""
        formula = pt.ChargedFormula(
            formula=(
                pt.FormulaElement(element=tacular.Element.C, occurance=2),
                pt.FormulaElement(element=tacular.Element.H, occurance=6),
            ),
            charge=2,
        )
        assert str(formula) == "Formula:C2H6:z+2"

    def test_charged_formula_negative_charge(self):
        """Test ChargedFormula with negative charge"""
        formula = pt.ChargedFormula(formula=(pt.FormulaElement(element=tacular.Element.O, occurance=1),), charge=-1)
        assert str(formula) == "Formula:O:z-1"


class TestGlycanStringConversion:
    """Tests for converting glycan objects to strings"""

    def test_glycan_component_simple(self):
        """Test GlycanComponent string conversion"""
        glycan = pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.Hex, occurance=1)
        assert str(glycan) == "Hex"

    def test_glycan_component_with_count(self):
        """Test GlycanComponent with count"""
        glycan = pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.Hex, occurance=5)
        assert str(glycan) == "Hex5"

    def test_glycan_tuple_to_string(self):
        """Test tuple of GlycanComponents"""
        glycan_tuple = (
            pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.Hex, occurance=5),
            pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.HexNAc, occurance=4),
        )
        glycan_str = "Glycan:" + "".join(str(g) for g in glycan_tuple)
        assert glycan_str == "Glycan:Hex5HexNAc4"


class TestPositionRuleStringConversion:
    """Tests for converting position rule objects to strings"""

    def test_position_rule_anywhere(self):
        """Test PositionRule with ANYWHERE terminal"""
        rule = pt.PositionRule(terminal=pt.Terminal.ANYWHERE, amino_acid=tacular.AminoAcid.M)
        assert str(rule) == "M"

    def test_position_rule_n_term(self):
        """Test PositionRule with N-term"""
        rule = pt.PositionRule(terminal=pt.Terminal.N_TERM)
        assert str(rule) == "N-term"

    def test_position_rule_n_term_with_aa(self):
        """Test PositionRule with N-term and amino acid"""
        rule = pt.PositionRule(terminal=pt.Terminal.N_TERM, amino_acid=tacular.AminoAcid.K)
        assert str(rule) == "N-term:K"


class TestModificationStringConversion:
    """Tests for converting modification objects to strings"""

    def test_modification_ambiguous_secondary(self):
        """Test ModificationAmbiguousSecondary string conversion"""
        mod = pt.ModificationAmbiguousSecondary(label="g1")
        assert str(mod) == "#g1"

    def test_modification_ambiguous_secondary_with_score(self):
        """Test ModificationAmbiguousSecondary with score"""
        mod = pt.ModificationAmbiguousSecondary(label="g1", score=0.95)
        assert str(mod) == "#g1(0.95)"

    def test_modification_ambiguous_primary_simple(self):
        """Test ModificationAmbiguousPrimary string conversion"""
        tags = pt.ModificationTags.from_string("Phospho")
        mod = pt.ModificationAmbiguousPrimary(label="g1", tags=tags)
        assert str(mod) == "Phospho#g1"

    def test_modification_cross_linker_primary(self):
        """Test ModificationCrossLinker primary definition"""
        tags = pt.ModificationTags.from_string("XLMOD:02001")
        assert str(tags) == "XLMOD:02001"
        mod = pt.ModificationCrossLinker(label="1", tags=tags)
        assert str(mod) == "XLMOD:02001#XL1"

    def test_modification_cross_linker_secondary(self):
        """Test ModificationCrossLinker secondary reference"""
        mod = pt.ModificationCrossLinker(label="1")
        assert str(mod) == "#XL1"


class TestRoundTrip:
    """Tests for round-trip conversion: parse -> string -> parse"""

    def test_round_trip_simple_peptide(self):
        """Test round-trip: parse modification and convert back to string"""
        original = "Oxidation"
        parsed = pt.ModificationTags.from_string(original).tags[0]
        assert str(parsed) == original

    def test_round_trip_accession(self):
        """Test round-trip for accession"""
        original = "Unimod:35"
        parsed = pt.ModificationTags.from_string(original).tags[0]
        result_str = str(parsed)
        # Note: case might differ (UNIMOD vs Unimod)
        assert result_str.upper() == original.upper()

    def test_round_trip_mass(self):
        """Test round-trip for mass"""
        original = "+15.995"
        parsed = pt.ModificationTags.from_string(original).tags[0]
        assert str(parsed) == original

    def test_round_trip_formula(self):
        """Test round-trip for formula"""
        original = "Formula:C2H6"
        parsed = pt.ModificationTags.from_string(original).tags[0]
        assert str(parsed) == original

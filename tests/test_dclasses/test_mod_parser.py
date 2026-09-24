"""
Tests for modification parsing functionality.

Tests the conversion of ProForma modification strings into structured tag objects.
"""

import pytest
import tacular

import peptacular as pt


class TestTagName:
    """Tests for parsing named modifications"""

    def test_simple_name(self):
        """Test parsing a simple modification name"""
        result = pt.ModificationTags.from_string("Oxidation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Oxidation"
        assert result.cv is None

    def test_various_names(self):
        """Test parsing various modification names"""
        names = ["Phospho", "Acetyl", "Amidated", "Methylation"]
        for name in names:
            result = pt.ModificationTags.from_string(name).tags[0]
            assert isinstance(result, pt.TagName)
            assert result.name == name

    def test_case_sensitive_name(self):
        """Test that modification names are case-sensitive"""
        result1 = pt.ModificationTags.from_string("Phospho").tags[0]
        assert isinstance(result1, pt.TagName)
        assert result1.name == "Phospho"

        result2 = pt.ModificationTags.from_string("phospho").tags[0]
        assert isinstance(result2, pt.TagName)
        assert result2.name == "phospho"


class TestTagAccession:
    """Tests for parsing CV accession modifications"""

    def test_unimod_accession(self):
        """Test parsing UNIMOD accession"""
        result = pt.ModificationTags.from_string("UNIMOD:35").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "35"
        assert result.cv == pt.CV.UNIMOD

    def test_unimod_case_insensitive(self):
        """Test that UNIMOD is case-insensitive"""
        result = pt.ModificationTags.from_string("unimod:35").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.cv == pt.CV.UNIMOD

    def test_psi_mod_accession(self):
        """Test parsing PSI-MOD accession"""
        result = pt.ModificationTags.from_string("MOD:00719").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "00719"
        assert result.cv == pt.CV.PSI_MOD

    def test_resid_accession(self):
        """Test parsing RESID accession"""
        result = pt.ModificationTags.from_string("RESID:AA0037").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "AA0037"
        assert result.cv == pt.CV.RESID

    def test_gnome_accession(self):
        """Test parsing GNO accession (use GNO:)"""
        result = pt.ModificationTags.from_string("GNO:G12345").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "G12345"
        assert result.cv == pt.CV.GNOME

    def test_xlmod_accession(self):
        """Test parsing XL-MOD accession"""
        result = pt.ModificationTags.from_string("XLMOD:02001").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "02001"
        assert result.cv == pt.CV.XL_MOD

    def test_xlmod_short_form(self):
        """Test parsing XLMOD (without hyphen) accession"""
        result = pt.ModificationTags.from_string("XLMOD:02001").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "02001"
        assert result.cv == pt.CV.XL_MOD


class TestShorthandForNamedMods:
    """Tests for single-letter CV prefixes with named modifications (Rule 4)"""

    def test_unimod_shorthand(self):
        """Test parsing UNIMOD shorthand U:Oxidation"""
        result = pt.ModificationTags.from_string("U:Oxidation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Oxidation"
        assert result.cv == pt.CV.UNIMOD

    def test_unimod_shorthand_lowercase(self):
        """Test parsing UNIMOD shorthand u:Phospho (case-insensitive)"""
        result = pt.ModificationTags.from_string("u:Phospho").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Phospho"
        assert result.cv == pt.CV.UNIMOD

    def test_psi_mod_shorthand(self):
        """Test parsing PSI-MOD shorthand M:Phospho"""
        result = pt.ModificationTags.from_string("M:Phospho").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Phospho"
        assert result.cv == pt.CV.PSI_MOD

    def test_psi_mod_shorthand_lowercase(self):
        """Test parsing PSI-MOD shorthand m:Acetyl (case-insensitive)"""
        result = pt.ModificationTags.from_string("m:Acetyl").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Acetyl"
        assert result.cv == pt.CV.PSI_MOD

    def test_resid_shorthand(self):
        """Test parsing RESID shorthand R:Methylation"""
        result = pt.ModificationTags.from_string("R:Methylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Methylation"
        assert result.cv == pt.CV.RESID

    def test_resid_shorthand_lowercase(self):
        """Test parsing RESID shorthand r:Biotinylation (case-insensitive)"""
        result = pt.ModificationTags.from_string("r:Biotinylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Biotinylation"
        assert result.cv == pt.CV.RESID

    def test_gnome_shorthand(self):
        """Test parsing GNO shorthand G:Glycosylation"""
        result = pt.ModificationTags.from_string("G:Glycosylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Glycosylation"
        assert result.cv == pt.CV.GNOME

    def test_gnome_shorthand_lowercase(self):
        """Test parsing GNO shorthand g:Mannosylation (case-insensitive)"""
        result = pt.ModificationTags.from_string("g:Mannosylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Mannosylation"
        assert result.cv == pt.CV.GNOME

    def test_xlmod_shorthand(self):
        """Test parsing XLMOD shorthand X:DSS"""
        result = pt.ModificationTags.from_string("X:DSS").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "DSS"
        assert result.cv == pt.CV.XL_MOD

    def test_xlmod_shorthand_lowercase(self):
        """Test parsing XLMOD shorthand x:BS3 (case-insensitive)"""
        result = pt.ModificationTags.from_string("x:BS3").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "BS3"
        assert result.cv == pt.CV.XL_MOD


class TestTagMass:
    """Tests for parsing mass delta modifications"""

    def test_positive_mass(self):
        """Test parsing positive mass delta"""
        result = pt.ModificationTags.from_string("+15.995").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(15.995)  # type: ignore
        assert result.cv is None

    def test_negative_mass(self):
        """Test parsing negative mass delta"""
        result = pt.ModificationTags.from_string("-18.010").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(-18.010)  # type: ignore

    def test_mass_without_sign_invalid(self):
        """Test that mass without sign is treated as a name (not valid ProForma)"""
        result = pt.ModificationTags.from_string("42.010").tags[0]
        # Without a sign, it's not recognized as a mass, so it becomes a pt.TagName
        assert isinstance(result, pt.TagName)
        assert result.name == "42.010"

    def test_integer_mass(self):
        """Test parsing integer mass with required sign"""
        result = pt.ModificationTags.from_string("+16").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == 16.0


class TestChargedFormula:
    """Tests for parsing formula modifications"""

    def test_simple_formula(self):
        """Test parsing simple formula"""
        result = pt.ModificationTags.from_string("Formula:C2H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == tacular.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[1].element == tacular.Element.H
        assert result.formula[1].occurance == 6
        assert result.charge is None

    def test_single_element(self):
        """Test parsing single element formula"""
        result = pt.ModificationTags.from_string("Formula:O").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 1
        assert result.formula[0].element == tacular.Element.O
        assert result.formula[0].occurance == 1

    def test_formula_with_negative_count(self):
        """Test parsing formula with negative element counts"""
        result = pt.ModificationTags.from_string("Formula:H-2O-1").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == tacular.Element.H
        assert result.formula[0].occurance == -2
        assert result.formula[1].element == tacular.Element.O
        assert result.formula[1].occurance == -1

    def test_formula_with_charge(self):
        """Test parsing formula with charge state"""
        result = pt.ModificationTags.from_string("Formula:C2H6:z+2").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.charge == 2

    def test_formula_with_negative_charge(self):
        """Test parsing formula with negative charge"""
        result = pt.ModificationTags.from_string("Formula:O:z-1").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert result.charge == -1

    def test_complex_formula(self):
        """Test parsing complex formula"""
        result = pt.ModificationTags.from_string("Formula:C10H15N3O6S").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 5
        # Check carbon
        assert result.formula[0].element == tacular.Element.C
        assert result.formula[0].occurance == 10

    def test_formula_with_spaces(self):
        """Test parsing formula with spaces between element pairs (ProForma Rule 1)"""
        result = pt.ModificationTags.from_string("Formula:C12 H20 O2").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 3
        assert result.formula[0].element == tacular.Element.C
        assert result.formula[0].occurance == 12
        assert result.formula[1].element == tacular.Element.H
        assert result.formula[1].occurance == 20
        assert result.formula[2].element == tacular.Element.O
        assert result.formula[2].occurance == 2

    def test_formula_isotope_prefix_notation(self):
        """Test parsing formula with isotope prefix notation [13C2] (ProForma Rule 3)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == tacular.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        assert result.formula[1].element == tacular.Element.H
        assert result.formula[1].occurance == 6

    def test_formula_multiple_isotopes(self):
        """Test parsing formula with multiple isotope specifications [13C2][12C-2]H2N"""
        result = pt.ModificationTags.from_string("Formula:[13C2][12C-2]H2N").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 4
        # First carbon: 13C with count 2
        assert result.formula[0].element == tacular.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        # Second carbon: 12C with count -2
        assert result.formula[1].element == tacular.Element.C
        assert result.formula[1].occurance == -2
        assert result.formula[1].isotope == 12
        # Hydrogen
        assert result.formula[2].element == tacular.Element.H
        assert result.formula[2].occurance == 2
        # Nitrogen
        assert result.formula[3].element == tacular.Element.N
        assert result.formula[3].occurance == 1

    def test_formula_isotope_replacement_example(self):
        """Test ProForma spec example: [13C2]C-2H2N (2 12C replaced by 2 13C)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]C-2H2N").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 4
        # 13C with count 2
        assert result.formula[0].element == tacular.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        # Natural C with count -2
        assert result.formula[1].element == tacular.Element.C
        assert result.formula[1].occurance == -2
        assert result.formula[1].isotope is None

    def test_formula_zero_cardinality_not_allowed(self):
        """Test that zero cardinality raises error (ProForma Rule 2)"""
        with pytest.raises(ValueError):
            pt.ModificationTags.from_string("Formula:C0H2")


class TestGlycanComposition:
    """Tests for parsing glycan composition modifications"""

    def test_simple_glycan(self):
        """Test parsing simple glycan composition"""
        result = pt.ModificationTags.from_string("Glycan:Hex").tags[0]
        assert isinstance(result, pt.GlycanTag)
        assert len(result) == 1
        res: pt.GlycanComponent = result[0]
        assert res.monosaccharide == tacular.Monosaccharide.Hex
        assert res.occurance == 1

    def test_glycan_with_count(self):
        """Test parsing glycan with count"""
        result = pt.ModificationTags.from_string("Glycan:Hex5").tags[0]
        assert isinstance(result, pt.GlycanTag)
        assert len(result) == 1
        res: pt.GlycanComponent = result[0]
        assert res.monosaccharide == tacular.Monosaccharide.Hex
        assert res.occurance == 5

    def test_complex_glycan_composition(self):
        """Test parsing complex glycan composition"""
        result = pt.ModificationTags.from_string("Glycan:Hex5HexNAc4").tags[0]
        assert isinstance(result, pt.GlycanTag)
        assert len(result) == 2
        res1: pt.GlycanComponent = result[0]
        res2: pt.GlycanComponent = result[1]

        assert res1.monosaccharide == tacular.Monosaccharide.Hex
        assert res1.occurance == 5
        assert res2.monosaccharide == tacular.Monosaccharide.HexNAc
        assert res2.occurance == 4

    def test_various_monosaccharides(self):
        """Test parsing various monosaccharide types"""
        monosaccharides = ["Fuc", "Hep", "NeuGc", "dHex"]
        for mono in monosaccharides:
            result = pt.ModificationTags.from_string(f"Glycan:{mono}").tags[0]
            assert isinstance(result, pt.GlycanTag)
            assert len(result) == 1


class TestTagInfo:
    """Tests for parsing INFO tag modifications"""

    def test_simple_info(self):
        """Test parsing simple INFO tag"""
        result = pt.ModificationTags.from_string("INFO:custom annotation").tags[0]
        assert isinstance(result, pt.TagInfo)
        assert result.info == "custom annotation"

    def test_info_with_special_chars(self):
        """Test parsing INFO tag with special characters"""
        result = pt.ModificationTags.from_string("INFO:annotation-with_special.chars").tags[0]
        assert isinstance(result, pt.TagInfo)
        assert result.info == "annotation-with_special.chars"


class TestEdgeCases:
    """Tests for edge cases and error conditions"""

    def test_empty_string(self):
        """Test parsing empty string raises error"""
        with pytest.raises(ValueError):
            pt.ModificationTags.from_string("")

    def test_whitespace_only(self):
        """Test parsing whitespace-only string raises error"""
        with pytest.raises(ValueError):
            pt.ModificationTags.from_string("   ")

    def test_unknown_cv_prefix_treated_as_name(self):
        """Test that unknown CV prefix is treated as a name"""
        result = pt.ModificationTags.from_string("UNKNOWN:123").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "UNKNOWN:123"

    def test_invalid_formula_element(self):
        """Test parsing formula with invalid element"""
        with pytest.raises(ValueError):
            pt.ModificationTags.from_string("Formula:Zz")

    def test_invalid_glycan_composition(self):
        """Test parsing invalid glycan composition"""
        with pytest.raises(ValueError):
            pt.ModificationTags.from_string("Glycan:InvalidMono")


class TestStringConversion:
    """Tests for converting dataclass objects back to ProForma strings"""

    def test_tag_name_to_string(self):
        """Test pt.TagName string conversion"""
        tag = pt.TagName(name="Oxidation")
        assert str(tag) == "Oxidation"

    def test_tag_name_with_cv_to_string(self):
        """Test pt.TagName with CV string conversion"""
        tag = pt.TagName(name="Oxidation", cv=pt.CV.UNIMOD)
        assert str(tag) == "Oxidation"

    def test_tag_accession_to_string(self):
        """Test pt.TagAccession string conversion"""
        tag = pt.TagAccession(accession="35", cv=pt.CV.UNIMOD)
        assert str(tag) == "UNIMOD:35"

    def test_tag_accession_psi_mod(self):
        """Test pt.TagAccession PSI-MOD string conversion"""
        tag = pt.TagAccession(accession="00719", cv=pt.CV.PSI_MOD)
        assert str(tag) == "MOD:00719"

    def test_tag_mass_positive(self):
        """Test pt.TagMass positive value string conversion"""
        tag = pt.TagMass(mass_str="15.995")
        assert str(tag) == "+15.995"

    def test_tag_mass_negative(self):
        """Test pt.TagMass negative value string conversion"""
        tag = pt.TagMass(mass_str="-18.010")
        assert str(tag) == "-18.01"

    def test_tag_mass_with_cv(self):
        """Test pt.TagMass with CV string conversion"""
        tag = pt.TagMass(mass_str="15.995", cv=pt.CV.UNIMOD)
        assert str(tag) == "U:+15.995"

    def test_tag_info_to_string(self):
        """Test pt.TagInfo string conversion"""
        tag = pt.TagInfo(info="custom annotation")
        assert str(tag) == "INFO:custom annotation"

    def test_formula_element_simple(self):
        """Test pt.FormulaElement string conversion"""
        elem = pt.FormulaElement(element=tacular.Element.C, occurance=2)
        assert str(elem) == "C2"

    def test_formula_element_single_count(self):
        """Test pt.FormulaElement with count of 1"""
        elem = pt.FormulaElement(element=tacular.Element.O, occurance=1)
        assert str(elem) == "O"

    def test_formula_element_with_isotope(self):
        """Test pt.FormulaElement with isotope"""
        elem = pt.FormulaElement(element=tacular.Element.C, occurance=2, isotope=13)
        assert str(elem) == "[13C2]"

    def test_charged_formula_simple(self):
        """Test pt.ChargedFormula string conversion"""
        formula = pt.ChargedFormula(
            formula=(
                pt.FormulaElement(element=tacular.Element.C, occurance=2),
                pt.FormulaElement(element=tacular.Element.H, occurance=6),
            )
        )
        assert str(formula) == "Formula:C2H6"

    def test_charged_formula_with_charge(self):
        """Test pt.ChargedFormula with charge"""
        formula = pt.ChargedFormula(
            formula=(
                pt.FormulaElement(element=tacular.Element.C, occurance=2),
                pt.FormulaElement(element=tacular.Element.H, occurance=6),
            ),
            charge=2,
        )
        assert str(formula) == "Formula:C2H6:z+2"

    def test_charged_formula_negative_charge(self):
        """Test pt.ChargedFormula with negative charge"""
        formula = pt.ChargedFormula(formula=(pt.FormulaElement(element=tacular.Element.O, occurance=1),), charge=-1)
        assert str(formula) == "Formula:O:z-1"

    def test_glycan_component_simple(self):
        """Test pt.GlycanComponent string conversion"""
        glycan = pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.Hex, occurance=1)
        assert str(glycan) == "Hex"

    def test_glycan_component_with_count(self):
        """Test pt.GlycanComponent with count"""
        glycan = pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.Hex, occurance=5)
        assert str(glycan) == "Hex5"

    def test_glycan_tuple_to_string(self):
        """Test tuple of pt.GlycanComponents"""
        glycan_tuple = (
            pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.Hex, occurance=5),
            pt.GlycanComponent(monosaccharide=tacular.Monosaccharide.HexNAc, occurance=4),
        )
        glycan_str = "Glycan:" + "".join(str(g) for g in glycan_tuple)
        assert glycan_str == "Glycan:Hex5HexNAc4"

    def test_position_rule_anywhere(self):
        """Test pt.PositionRule with ANYWHERE terminal"""

        rule = pt.PositionRule(terminal=pt.Terminal.ANYWHERE, amino_acid=tacular.AminoAcid.M)
        assert str(rule) == "M"

    def test_position_rule_n_term(self):
        """Test pt.PositionRule with N-term"""

        rule = pt.PositionRule(terminal=pt.Terminal.N_TERM)
        assert str(rule) == "N-term"

    def test_position_rule_n_term_with_aa(self):
        """Test pt.PositionRule with N-term and amino acid"""
        rule = pt.PositionRule(terminal=pt.Terminal.N_TERM, amino_acid=tacular.AminoAcid.K)
        assert str(rule) == "N-term:K"

    def test_sequence_element_simple(self):
        """Test pt.SequenceElement string conversion"""
        elem = pt.SequenceElement(amino_acid=tacular.AminoAcid.M)
        assert str(elem) == "M"

    def test_modification_ambiguous_secondary(self):
        """Test pt.ModificationAmbiguousSecondary string conversion"""
        mod = pt.ModificationAmbiguousSecondary(label="g1")
        assert str(mod) == "#g1"

    def test_modification_ambiguous_secondary_with_score(self):
        """Test pt.ModificationAmbiguousSecondary with score"""
        mod = pt.ModificationAmbiguousSecondary(label="g1", score=0.95)
        assert str(mod) == "#g1(0.95)"

    def test_modification_cross_linker_secondary(self):
        """Test pt.ModificationCrossLinker secondary reference"""
        mod = pt.ModificationCrossLinker(label="1")
        assert str(mod) == "#XL1"

    def test_isotope_replacement_deuterium(self):
        """Test IsotopeReplacement for deuterium"""
        iso = pt.IsotopeReplacement(element=tacular.Element.H, isotope=2)
        assert str(iso) == "D"

    def test_isotope_replacement_c13(self):
        """Test IsotopeReplacement for 13C"""
        iso = pt.IsotopeReplacement(element=tacular.Element.C, isotope=13)
        assert str(iso) == "13C"

    def test_round_trip_simple_peptide(self):
        """Test round-trip: parse modification and convert back to string"""
        original = "Oxidation"
        parsed = pt.ModificationTags.from_string(original)
        assert str(parsed) == original

    def test_round_trip_accession(self):
        """Test round-trip for accession"""
        original = "Unimod:35"
        parsed = pt.ModificationTags.from_string(original)
        result_str = str(parsed)
        # Note: case might differ (UNIMOD vs Unimod)
        assert result_str.upper() == original.upper()

    def test_round_trip_mass(self):
        """Test round-trip for mass"""
        original = "+15.995"
        parsed = pt.ModificationTags.from_string(original)
        assert str(parsed) == original

    def test_round_trip_formula(self):
        """Test round-trip for formula"""
        original = "Formula:C2H6"
        parsed = pt.ModificationTags.from_string(original)
        assert str(parsed) == original

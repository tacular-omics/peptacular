"""
Tests for parsing modification strings with multiple tags (pipe-separated).
"""

import pytest
import tacular

import peptacular as pt


class TestParseModificationTags:
    """Tests for ModificationTags.from_string function"""

    def test_returns_tuple(self):
        """Test that ModificationTags.from_string returns a tuple"""
        result = pt.ModificationTags.from_string("Oxidation")
        assert isinstance(result, pt.ModificationTags)
        assert len(result) == 1

    def test_contains_correct_tag(self):
        """Test that tuple contains the correct tag"""
        result = pt.ModificationTags.from_string("UNIMOD:35")
        assert len(result) == 1
        first_tag = result[0]
        assert isinstance(first_tag, pt.TagAccession)
        assert first_tag.accession == "35"


class TestMultipleTags:
    """Tests for parsing modification strings with multiple tags (pipe-separated)"""

    def test_two_tags_name_and_accession(self):
        """Test parsing two tags: name and accession"""
        result = pt.ModificationTags.from_string("Oxidation|UNIMOD:35")
        assert isinstance(result, pt.ModificationTags)
        assert len(result) == 2
        first_tag = result[0]
        second_tag = result[1]
        assert isinstance(first_tag, pt.TagName)
        assert first_tag.name == "Oxidation"
        assert isinstance(second_tag, pt.TagAccession)
        assert second_tag.accession == "35"
        assert second_tag.cv == pt.CV.UNIMOD

    def test_two_tags_accession_and_mass(self):
        """Test parsing accession and mass"""
        result = pt.ModificationTags.from_string("UNIMOD:35|+15.995")
        assert len(result) == 2
        first_tag = result[0]
        second_tag = result[1]
        assert isinstance(first_tag, pt.TagAccession)
        assert first_tag.accession == "35"
        assert isinstance(second_tag, pt.TagMass)
        assert second_tag.mass == pytest.approx(15.995)

    def test_three_tags_name_accession_mass(self):
        """Test parsing three tags: name, accession, and mass"""
        result = pt.ModificationTags.from_string("Oxidation|UNIMOD:35|+15.995")
        assert len(result) == 3
        assert isinstance(result[0], pt.TagName)
        assert result[0].name == "Oxidation"
        assert isinstance(result[1], pt.TagAccession)
        assert result[1].accession == "35"
        assert isinstance(result[2], pt.TagMass)
        assert result[2].mass == pytest.approx(15.995)

    def test_four_tags_complete_oxidation(self):
        """Test parsing all four ways to represent oxidation"""
        result = pt.ModificationTags.from_string("Oxidation|UNIMOD:35|+15.995|Formula:O")
        assert len(result) == 4
        assert isinstance(result[0], pt.TagName)
        assert isinstance(result[1], pt.TagAccession)
        assert isinstance(result[2], pt.TagMass)
        assert isinstance(result[3], pt.ChargedFormula)
        assert result[3].formula[0].element == tacular.Element.O

    def test_name_and_info(self):
        """Test parsing name with INFO tag"""
        result = pt.ModificationTags.from_string("Phospho|INFO:probable")
        assert len(result) == 2
        assert isinstance(result[0], pt.TagName)
        assert result[0].name == "Phospho"
        assert isinstance(result[1], pt.TagInfo)
        assert result[1].info == "probable"

    def test_accession_and_info(self):
        """Test parsing accession with INFO tag"""
        result = pt.ModificationTags.from_string("MOD:00046|INFO:confident")
        assert len(result) == 2
        assert isinstance(result[0], pt.TagAccession)
        assert result[0].cv == pt.CV.PSI_MOD
        assert isinstance(result[1], pt.TagInfo)
        assert result[1].info == "confident"

    def test_shorthand_with_mass(self):
        """Test parsing shorthand named mod with mass (Rule 4 + Rule 2)"""
        result = pt.ModificationTags.from_string("U:Oxidation|+15.995")
        assert len(result) == 2
        assert isinstance(result[0], pt.TagName)
        assert result[0].cv == pt.CV.UNIMOD
        assert result[0].name == "Oxidation"
        assert isinstance(result[1], pt.TagMass)
        assert result[1].mass == pytest.approx(15.995)

    def test_formula_with_mass(self):
        """Test parsing formula with mass"""
        result = pt.ModificationTags.from_string("Formula:H-2O-1|+79.966")
        assert len(result) == 2
        assert isinstance(result[0], pt.ChargedFormula)
        assert isinstance(result[1], pt.TagMass)

    def test_whitespace_handling(self):
        """Test that whitespace around pipe is handled correctly"""
        result = pt.ModificationTags.from_string("Oxidation | UNIMOD:35 | +15.995")
        assert len(result) == 3
        assert isinstance(result[0], pt.TagName)
        assert result[0].name == "Oxidation"
        assert isinstance(result[1], pt.TagAccession)
        assert result[1].accession == "35"
        assert isinstance(result[2], pt.TagMass)
        assert result[2].mass == pytest.approx(15.995)

    def test_empty_parts_ignored(self):
        """Test that empty parts from consecutive pipes throw error"""
        with pytest.raises(ValueError):
            _ = pt.ModificationTags.from_string("Oxidation||UNIMOD:35")

    def test_complex_real_world_example(self):
        """Test complex real-world example with multiple tags"""
        result = pt.ModificationTags.from_string("Phosphorylation|MOD:00046|M:Phospho|+79.966|Formula:H P O3|INFO:high confidence")
        assert len(result) == 6
        assert isinstance(result[0], pt.TagName)
        assert result[0].name == "Phosphorylation"
        assert isinstance(result[1], pt.TagAccession)
        assert result[1].cv == pt.CV.PSI_MOD
        assert isinstance(result[2], pt.TagName)
        assert result[2].cv == pt.CV.PSI_MOD  # M: is shorthand for named mods
        assert result[2].name == "Phospho"
        assert isinstance(result[3], pt.TagMass)
        assert isinstance(result[4], pt.ChargedFormula)
        assert isinstance(result[5], pt.TagInfo)

    def test_glycan_with_info(self):
        """Test glycan composition with INFO tag"""
        result = pt.ModificationTags.from_string("Glycan:Hex5HexNAc4|INFO:N-glycan")
        assert len(result) == 2
        # First tag should be GlycanTag
        assert isinstance(result[0], pt.GlycanTag)
        assert len(result[0].components) == 2
        assert all(isinstance(g, pt.GlycanComponent) for g in result[0].components)
        assert isinstance(result[1], pt.TagInfo)

    def test_multiple_accessions_different_cvs(self):
        """Test multiple accessions from different CVs"""
        result = pt.ModificationTags.from_string("UNIMOD:35|MOD:00719|RESID:AA0037")
        assert len(result) == 3
        assert isinstance(result[0], pt.TagAccession)
        assert result[0].cv == pt.CV.UNIMOD
        assert isinstance(result[1], pt.TagAccession)
        assert result[1].cv == pt.CV.PSI_MOD
        assert isinstance(result[2], pt.TagAccession)
        assert result[2].cv == pt.CV.RESID

    def test_mixed_shorthand_and_full(self):
        """Test mixing shorthand named mods and full accessions"""
        result = pt.ModificationTags.from_string("U:Oxidation|MOD:00719|M:Phospho")
        assert len(result) == 3
        assert isinstance(result[0], pt.TagName)
        assert result[0].cv == pt.CV.UNIMOD
        assert isinstance(result[1], pt.TagAccession)
        assert result[1].cv == pt.CV.PSI_MOD
        assert isinstance(result[2], pt.TagName)
        assert result[2].cv == pt.CV.PSI_MOD

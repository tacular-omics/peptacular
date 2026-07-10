"""
Tests for parsing formula modification tags.
"""

import pytest

import peptacular as pt


class TestChargedFormula:
    """Tests for parsing formula modifications"""

    def test_simple_formula(self):
        """Test parsing simple formula"""
        result = pt.ModificationTags.from_string("Formula:C2H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 6
        assert result.charge is None

    def test_single_element(self):
        """Test parsing single element formula"""
        result = pt.ModificationTags.from_string("Formula:O").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 1
        assert result.formula[0].element == pt.Element.O
        assert result.formula[0].occurance == 1

    def test_formula_with_negative_count(self):
        """Test parsing formula with negative element counts"""
        result = pt.ModificationTags.from_string("Formula:H-2O-1").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.H
        assert result.formula[0].occurance == -2
        assert result.formula[1].element == pt.Element.O
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

    def test_formula_with_isotope(self):
        """Test parsing formula with isotope specification [13C2] (count inside bracket)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 6

    def test_complex_formula(self):
        """Test parsing complex formula"""
        result = pt.ModificationTags.from_string("Formula:C10H15N3O6S").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 5
        # Check carbon
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 10

    def test_formula_with_spaces(self):
        """Test parsing formula with spaces between element pairs (ProForma Rule 1)"""
        result = pt.ModificationTags.from_string("Formula:C12 H20 O2").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 3
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 12
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 20
        assert result.formula[2].element == pt.Element.O
        assert result.formula[2].occurance == 2

    def test_formula_isotope_prefix_notation(self):
        """Test parsing formula with isotope prefix notation [13C2] (ProForma Rule 3)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 6

    def test_formula_isotope_single_count(self):
        """Test parsing isotope with default count of 1: [13C]"""
        result = pt.ModificationTags.from_string("Formula:[13C]H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 1
        assert result.formula[0].isotope == 13

    def test_formula_multiple_isotopes(self):
        """Test parsing formula with multiple isotope specifications [13C2][12C-2]H2N"""
        result = pt.ModificationTags.from_string("Formula:[13C2][12C-2]H2N").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 4
        # First carbon: 13C with count 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        # Second carbon: 12C with count -2
        assert result.formula[1].element == pt.Element.C
        assert result.formula[1].occurance == -2
        assert result.formula[1].isotope == 12
        # Hydrogen
        assert result.formula[2].element == pt.Element.H
        assert result.formula[2].occurance == 2
        # Nitrogen
        assert result.formula[3].element == pt.Element.N
        assert result.formula[3].occurance == 1

    def test_formula_isotope_replacement_example(self):
        """Test ProForma spec example: [13C2]C-2H2N (2 12C replaced by 2 13C)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]C-2H2N").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 4
        # 13C with count 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        # Natural C with count -2
        assert result.formula[1].element == pt.Element.C
        assert result.formula[1].occurance == -2
        assert result.formula[1].isotope is None

    def test_formula_zero_cardinality_not_allowed(self):
        """Test that zero cardinality raises error (ProForma Rule 2)"""
        with pytest.raises(ValueError):
            pt.ModificationTags.from_string("Formula:C0H2").tags[0]

    def test_formula_deuterium_shorthand_keeps_isotope(self):
        """Deuterium 'D' in a formula must carry isotope H-2, not collapse to protium."""
        result = pt.ModificationTags.from_string("Formula:CD3").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        # C then D3 -> hydrogen with isotope 2, count 3
        h = next(fe for fe in result.formula if fe.element == pt.Element.H)
        assert h.isotope == 2
        assert h.occurance == 3
        # Mass must equal the explicit [2H3] isotope form, not natural H3
        assert abs(pt.mass("PEPT[Formula:CD3]IDE") - pt.mass("PEPT[Formula:C[2H3]]IDE")) < 1e-9
        assert pt.mass("PEPT[Formula:CD3]IDE") != pytest.approx(pt.mass("PEPT[Formula:CH3]IDE"))
        # Round-trips back to the D shorthand
        assert pt.parse("PEPT[Formula:CD3]IDE").serialize() == "PEPT[Formula:CD3]IDE"

    def test_formula_tritium_shorthand_keeps_isotope(self):
        """Tritium 'T' in a formula must carry isotope H-3."""
        result = pt.ModificationTags.from_string("Formula:CT3").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        h = next(fe for fe in result.formula if fe.element == pt.Element.H)
        assert h.isotope == 3
        assert h.occurance == 3
        assert abs(pt.mass("PEPT[Formula:CT3]IDE") - pt.mass("PEPT[Formula:C[3H3]]IDE")) < 1e-9
        assert pt.parse("PEPT[Formula:CT3]IDE").serialize() == "PEPT[Formula:CT3]IDE"

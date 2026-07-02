"""
Tests for chem.py module functions.
"""

import pytest
from tacular import ELEMENT_LOOKUP, ElementInfo

import peptacular as pt


class TestParseFormula:
    """Test parse_formula function."""

    def test_parse_simple_formula(self):
        """Parse a simple chemical formula string."""
        result = pt.parse_formula("H2O")
        assert result[ELEMENT_LOOKUP["H"]] == 2
        assert result[ELEMENT_LOOKUP["O"]] == 1

    def test_parse_complex_formula(self):
        """Parse a complex formula with multiple elements."""
        result = pt.parse_formula("C6H12O6")
        assert result[ELEMENT_LOOKUP["C"]] == 6
        assert result[ELEMENT_LOOKUP["H"]] == 12
        assert result[ELEMENT_LOOKUP["O"]] == 6

    def test_parse_formula_from_mapping(self):
        """Parse formula from a mapping of elements to counts."""
        formula_map = {"C": 6, "H": 12, "O": 6}
        result = pt.parse_formula(formula_map)
        assert result[ELEMENT_LOOKUP["C"]] == 6
        assert result[ELEMENT_LOOKUP["H"]] == 12
        assert result[ELEMENT_LOOKUP["O"]] == 6

    def test_parse_formula_from_element_info_mapping(self):
        """Parse formula from a mapping with ElementInfo objects."""
        formula_map = {
            ELEMENT_LOOKUP["C"]: 6,
            ELEMENT_LOOKUP["H"]: 12,
            ELEMENT_LOOKUP["O"]: 6,
        }
        result = pt.parse_formula(formula_map)
        assert result[ELEMENT_LOOKUP["C"]] == 6
        assert result[ELEMENT_LOOKUP["H"]] == 12
        assert result[ELEMENT_LOOKUP["O"]] == 6

    def test_parse_formula_list_sequential(self):
        """Parse multiple formulas sequentially."""
        formulas = ["H2O", "CO2", "NH3"]
        results = pt.parse_formula(formulas, n_workers=None, method="sequential")
        assert results[0][ELEMENT_LOOKUP["H"]] == 2
        assert results[0][ELEMENT_LOOKUP["O"]] == 1
        assert results[1][ELEMENT_LOOKUP["C"]] == 1
        assert results[1][ELEMENT_LOOKUP["O"]] == 2
        assert results[2][ELEMENT_LOOKUP["N"]] == 1
        assert results[2][ELEMENT_LOOKUP["H"]] == 3

    def test_parse_formula_list_parallel(self):
        """Parse multiple formulas in parallel."""
        formulas = ["H2O", "CO2", "NH3", "CH4"]
        results = pt.parse_formula(formulas, n_workers=2, method="thread")
        assert len(results) == 4
        assert results[0][ELEMENT_LOOKUP["H"]] == 2
        assert results[3][ELEMENT_LOOKUP["C"]] == 1
        assert results[3][ELEMENT_LOOKUP["H"]] == 4

    def test_parse_formula_with_separator(self):
        """Parse formula with custom separator."""
        result = pt.parse_formula("C 6 H 12 O 6", sep=" ")
        assert result[ELEMENT_LOOKUP["C"]] == 6
        assert result[ELEMENT_LOOKUP["H"]] == 12
        assert result[ELEMENT_LOOKUP["O"]] == 6

    def test_parse_formula_empty(self):
        """Parse an empty-like formula."""
        result = pt.parse_formula({})
        assert len(result) == 0


class TestChemComp:
    """Test chem_comp function."""

    def test_chem_comp_simple(self):
        """Get composition from simple formula."""
        result = pt.chem_comp("H2O")
        assert result[ELEMENT_LOOKUP["H"]] == 2
        assert result[ELEMENT_LOOKUP["O"]] == 1

    def test_chem_comp_from_mapping(self):
        """Get composition from mapping."""
        formula_map = {"N": 2, "H": 6}
        result = pt.chem_comp(formula_map)
        assert result[ELEMENT_LOOKUP["N"]] == 2
        assert result[ELEMENT_LOOKUP["H"]] == 6

    def test_chem_comp_list(self):
        """Get compositions from list of formulas."""
        formulas = ["H2O", "CO2", "NH3"]
        results = pt.chem_comp(formulas, n_workers=None, method="sequential")
        assert len(results) == 3
        assert results[0][ELEMENT_LOOKUP["H"]] == 2
        assert results[1][ELEMENT_LOOKUP["O"]] == 2
        assert results[2][ELEMENT_LOOKUP["N"]] == 1

    def test_chem_comp_parallel(self):
        """Get compositions in parallel."""
        formulas = ["CH4", "C2H6", "C3H8", "C4H10"]
        results = pt.chem_comp(formulas, n_workers=2, method="thread")
        assert len(results) == 4
        assert results[0][ELEMENT_LOOKUP["C"]] == 1
        assert results[1][ELEMENT_LOOKUP["C"]] == 2
        assert results[2][ELEMENT_LOOKUP["C"]] == 3
        assert results[3][ELEMENT_LOOKUP["C"]] == 4


class TestChemMass:
    """Test chem_mass function."""

    def test_chem_mass_water(self):
        """Calculate mass of water."""
        mass = pt.chem_mass("H2O", monoisotopic=True)
        # H: 1.007825, O: 15.994915
        expected = 2 * 1.007825 + 15.994915
        assert abs(mass - expected) < 0.001

    def test_chem_mass_glucose(self):
        """Calculate mass of glucose."""
        mass = pt.chem_mass("C6H12O6", monoisotopic=True)
        assert abs(mass - 180.063) < 0.01

    def test_chem_mass_from_mapping(self):
        """Calculate mass from element mapping."""
        formula_map = {"C": 1, "H": 4}  # Methane
        mass = pt.chem_mass(formula_map, monoisotopic=True)
        assert abs(mass - 16.031) < 0.01

    def test_chem_mass_average(self):
        """Calculate average mass instead of monoisotopic."""
        mass_mono = pt.chem_mass("H2O", monoisotopic=True)
        mass_avg = pt.chem_mass("H2O", monoisotopic=False)
        # Average mass should be slightly higher
        assert mass_avg > mass_mono

    def test_chem_mass_list(self):
        """Calculate masses of multiple formulas."""
        formulas = ["H2O", "CO2", "NH3"]
        masses = pt.chem_mass(formulas, monoisotopic=True, n_workers=None, method="sequential")
        assert len(masses) == 3
        assert all(isinstance(m, float) for m in masses)
        assert abs(masses[0] - 18.015) < 0.01  # H2O
        assert abs(masses[1] - 43.990) < 0.01  # CO2
        assert abs(masses[2] - 17.027) < 0.01  # NH3

    def test_chem_mass_parallel(self):
        """Calculate masses in parallel."""
        formulas = ["CH4", "C2H6", "C3H8", "C4H10"]
        masses = pt.chem_mass(formulas, monoisotopic=True, n_workers=2, method="thread")
        assert len(masses) == 4
        # Each successive alkane should be heavier
        assert masses[0] < masses[1] < masses[2] < masses[3]


class TestChemFormula:
    """Test chem_formula function."""

    def test_chem_formula_from_mapping(self):
        """Generate formula string from element mapping."""
        comp_map = {"C": 6, "H": 12, "O": 6}
        result = pt.chem_formula(comp_map, hill_order=True)
        # Hill order: C, H, then alphabetical
        assert result == "C6H12O6"

    def test_chem_formula_hill_order(self):
        """Test Hill order (C first, then H, then alphabetical)."""
        comp_map = {"O": 2, "H": 4, "C": 2}
        result = pt.chem_formula(comp_map, hill_order=True)
        assert result == "C2H4O2"

    def test_chem_formula_no_hill_order(self):
        """Test without Hill order."""
        comp_map = {"O": 2, "H": 4, "C": 2}
        result = pt.chem_formula(comp_map, hill_order=False)
        # When Hill order is disabled, still gets some ordering
        assert "C2" in result and "H4" in result and "O2" in result

    def test_chem_formula_with_separator(self):
        """Generate formula with custom separator."""
        comp_map = {"C": 6, "H": 12, "O": 6}
        result = pt.chem_formula(comp_map, hill_order=True, sep=" ")
        assert " " in result
        assert "C" in result and "H" in result and "O" in result

    def test_chem_formula_with_prefix(self):
        """Generate formula with Formula: prefix."""
        comp_map = {"C": 6, "H": 12, "O": 6}
        result = pt.chem_formula(comp_map, include_formula_prefix=True)
        assert result.startswith("Formula:")

    def test_chem_formula_list(self):
        """Generate formula strings from list of compositions."""
        comps = [{"H": 2, "O": 1}, {"C": 1, "O": 2}, {"N": 1, "H": 3}]
        results = pt.chem_formula(comps, hill_order=True, n_workers=None, method="sequential")
        assert len(results) == 3
        assert results[0] == "H2O"
        assert results[1] == "CO2"
        assert results[2] == "H3N"

    def test_chem_formula_parallel(self):
        """Generate formula strings in parallel."""
        comps = [
            {"C": 1, "H": 4},
            {"C": 2, "H": 6},
            {"C": 3, "H": 8},
            {"C": 4, "H": 10},
        ]
        results = pt.chem_formula(comps, n_workers=2, method="thread")
        assert len(results) == 4
        assert results[0] == "CH4"
        assert results[1] == "C2H6"
        assert results[2] == "C3H8"
        assert results[3] == "C4H10"

    def test_chem_formula_single_element(self):
        """Generate formula for single element."""
        comp_map = {"N": 2}
        result = pt.chem_formula(comp_map)
        assert result == "N2"


class TestChemIntegration:
    """Integration tests combining multiple chem functions."""

    def test_round_trip_parse_and_generate(self):
        """Parse a formula and generate it back."""
        original = "C6H12O6"
        comp = pt.parse_formula(original)
        regenerated = pt.chem_formula(comp, hill_order=True)
        assert regenerated == original

    def test_mass_from_parsed_formula(self):
        """Parse formula and calculate its mass."""
        formula = "C6H12O6"
        comp = pt.parse_formula(formula)
        mass = pt.chem_mass(comp, monoisotopic=True)
        assert abs(mass - 180.063) < 0.01

    def test_workflow_multiple_formulas(self):
        """Test complete workflow with multiple formulas."""
        formulas = ["H2O", "CO2", "NH3"]

        # Parse all formulas
        comps = pt.parse_formula(formulas, n_workers=None, method="sequential")
        assert len(comps) == 3

        # Calculate masses
        masses = pt.chem_mass(formulas, monoisotopic=True, n_workers=None, method="sequential")
        assert len(masses) == 3

        # Generate formula strings from compositions
        regenerated = pt.chem_formula(comps, hill_order=True, n_workers=None, method="sequential")
        assert regenerated[0] == "H2O"
        assert regenerated[1] == "CO2"
        assert regenerated[2] == "H3N"

    def test_composition_arithmetic(self):
        """Test that compositions can be combined."""
        comp1 = pt.parse_formula("H2O")
        comp2 = pt.parse_formula("CO2")

        # Add compositions together
        combined = comp1 + comp2
        formula = pt.chem_formula(combined, hill_order=True)
        assert formula == "CH2O3"

        # Calculate combined mass
        mass = pt.chem_mass(combined, monoisotopic=True)
        expected = pt.chem_mass("H2O") + pt.chem_mass("CO2")
        assert abs(mass - expected) < 0.001

    def test_element_info_consistency(self):
        """Test that ElementInfo objects are used consistently."""
        comp = pt.parse_formula("C6H12O6")

        # All keys should be ElementInfo instances
        assert all(isinstance(k, ElementInfo) for k in comp.keys())

        # Should be able to look up by string
        assert comp[ELEMENT_LOOKUP["C"]] == 6


class TestChemEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_composition(self):
        """Handle empty composition."""
        comp = {}
        formula = pt.chem_formula(comp)
        assert formula == ""

    def test_single_count_element(self):
        """Test element with count of 1."""
        comp = {"O": 1}
        formula = pt.chem_formula(comp)
        # Count of 1 might be shown as "O" or "O1" depending on implementation
        assert "O" in formula

    def test_large_element_counts(self):
        """Handle large element counts."""
        comp = {"C": 100, "H": 200}
        formula = pt.chem_formula(comp)
        assert "C100" in formula
        assert "H200" in formula
        mass = pt.chem_mass(comp)
        assert mass > 0

    def test_chem_formula_single_string(self):
        """A single formula string must not be treated as a batch of characters.

        Regression: chem_formula used ``isinstance(comp, Sequence)`` without excluding
        ``str``, so ``chem_formula("H2O")`` iterated the string character-by-character
        and raised on the digit ``"2"``.
        """
        assert pt.chem_formula("H2O") == "H2O"
        assert pt.chem_formula("C6H12O6") == "C6H12O6"

    def test_chem_formula_string_matches_dict(self):
        """String and equivalent dict inputs produce the same formula."""
        assert pt.chem_formula("H2O") == pt.chem_formula({"H": 2, "O": 1})

    def test_chem_formula_list_still_batches(self):
        """A list of formulas still returns a list of results."""
        assert pt.chem_formula(["H2O", "CO2"]) == ["H2O", "CO2"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

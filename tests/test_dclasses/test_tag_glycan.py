"""
Tests for parsing glycan composition modification tags.
"""

import peptacular as pt


class TestGlycanComposition:
    """Tests for parsing glycan composition modifications"""

    def test_simple_glycan(self):
        """Test parsing simple glycan composition"""
        result = pt.ModificationTags.from_string("Glycan:Hex").tags[0]
        # assert that component is GlycanTag
        assert isinstance(result, pt.GlycanTag)
        res: pt.GlycanComponent = result.components[0]
        assert res.monosaccharide == pt.Monosaccharide.Hex
        assert res.occurance == 1

    def test_glycan_with_count(self):
        """Test parsing glycan with count"""
        result = pt.ModificationTags.from_string("Glycan:Hex5").tags[0]
        assert isinstance(result, pt.GlycanTag)
        res: pt.GlycanComponent = result.components[0]
        assert res.monosaccharide == pt.Monosaccharide.Hex
        assert res.occurance == 5

    def test_complex_glycan_composition(self):
        """Test parsing complex glycan composition"""
        result = pt.ModificationTags.from_string("Glycan:Hex5HexNAc4").tags[0]
        assert isinstance(result, pt.GlycanTag)
        assert len(result) == 2
        res1: pt.GlycanComponent = result.components[0]
        res2: pt.GlycanComponent = result.components[1]

        assert res1.monosaccharide == pt.Monosaccharide.Hex
        assert res1.occurance == 5
        assert res2.monosaccharide == pt.Monosaccharide.HexNAc
        assert res2.occurance == 4

    def test_various_monosaccharides(self):
        """Test parsing various monosaccharide types"""
        monosaccharides = ["Fuc", "Hep", "NeuGc", "dHex"]
        for mono in monosaccharides:
            result = pt.ModificationTags.from_string(f"Glycan:{mono}").tags[0]
            assert isinstance(result, pt.GlycanTag)
            assert len(result) == 1


class TestParseGlycan:
    """Tests for parse_glycan() function"""

    def test_simple_glycan(self):
        """Test parsing simple glycan with prefix"""
        from peptacular.proforma_components.parsers import parse_glycan

        result = parse_glycan("Glycan:Hex5")
        assert len(result) == 1
        assert result[0].monosaccharide == pt.Monosaccharide.Hex
        assert result[0].occurance == 5

    def test_complex_glycan(self):
        """Test parsing complex glycan composition"""
        from peptacular.proforma_components.parsers import parse_glycan

        result = parse_glycan("Glycan:Hex5HexNAc4NeuAc2")
        assert isinstance(result, tuple)
        assert len(result) == 3
        assert result[0].monosaccharide == pt.Monosaccharide.Hex
        assert result[0].occurance == 5
        assert result[1].monosaccharide == pt.Monosaccharide.HexNAc
        assert result[1].occurance == 4
        assert result[2].monosaccharide == pt.Monosaccharide.NeuAc
        assert result[2].occurance == 2

    def test_case_insensitive_prefix(self):
        """Test that Glycan: prefix is case insensitive"""
        from peptacular.proforma_components.parsers import parse_glycan

        result1 = parse_glycan("Glycan:Hex")
        result2 = parse_glycan("glycan:Hex")
        result3 = parse_glycan("GLYCAN:Hex")
        assert result1 == result2 == result3

    def test_missing_prefix_raises_error(self):
        """Test that missing Glycan: prefix raises ValueError"""
        import pytest

        from peptacular.proforma_components.parsers import parse_glycan

        with pytest.raises(ValueError):
            parse_glycan("Hex5HexNAc4")

    def test_empty_string_raises_error(self):
        """Test that empty string raises ValueError"""
        import pytest

        from peptacular.proforma_components.parsers import parse_glycan

        with pytest.raises(ValueError):
            parse_glycan("")

    def test_only_prefix_raises_error(self):
        """Test that only prefix without composition raises ValueError"""
        import pytest

        from peptacular.proforma_components.parsers import parse_glycan

        # This should fail during composition parsing
        with pytest.raises(ValueError):
            parse_glycan("Glycan:")


class TestGlycanWhitespace:
    """ProForma allows optional whitespace between monosaccharide/count tokens."""

    def test_space_separated_components(self):
        from peptacular.proforma_components.parsers import parse_glycan

        spaced = parse_glycan("Glycan:Hex5 HexNAc4")
        compact = parse_glycan("Glycan:Hex5HexNAc4")
        assert spaced == compact

    def test_spaced_and_compact_masses_match(self):
        assert pt.parse("N[Glycan:Hex5 HexNAc4]").mass() == pt.parse("N[Glycan:Hex5HexNAc4]").mass()

    def test_multiple_and_surrounding_whitespace(self):
        from peptacular.proforma_components.parsers import parse_glycan

        result = parse_glycan("Glycan: Hex5  HexNAc4  NeuAc2 ")
        assert [(c.monosaccharide, c.occurance) for c in result] == [
            (pt.Monosaccharide.Hex, 5),
            (pt.Monosaccharide.HexNAc, 4),
            (pt.Monosaccharide.NeuAc, 2),
        ]

    def test_whitespace_between_name_and_count(self):
        from peptacular.proforma_components.parsers import parse_glycan

        assert parse_glycan("Glycan:Hex 5") == parse_glycan("Glycan:Hex5")

    def test_formula_in_glycan_still_rejected(self):
        """A chemical formula embedded in a glycan is not valid ProForma (reference rejects it too)."""
        import pytest

        from peptacular.proforma_components.parsers import parse_glycan

        with pytest.raises(ValueError):
            parse_glycan("Glycan:Hex1Formula:CH2")

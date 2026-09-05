"""
Tests for parsing glycan composition modification tags.
"""

import pytest

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


class TestGlycanFormulaAndMassComponents:
    """ProForma 2.1 §10.2: components not in the named list are given as a molecular
    formula or monoisotopic mass wrapped in curly braces, intermixed with monosaccharides."""

    SPEC_EXAMPLES = [
        "SEQUEN[Glycan:{C8H13N1O5}1Hex2]CE",  # molecular formula
        "SEQUEN[Glycan:{C8H13[15N1]O5}1Hex2]CE",  # isotope-labelled formula
        "SEQUEN[Glycan:{C8H13N1O5Na1:z+1}1Hex2]CE",  # charged formula (level 3)
        "SEQUEN[Glycan:{+203.079}1Hex2]CE",  # bare monoisotopic mass
    ]

    @pytest.mark.parametrize("seq", SPEC_EXAMPLES)
    def test_spec_examples_round_trip(self, seq):
        assert pt.parse(seq).serialize() == seq

    @pytest.mark.parametrize("seq", SPEC_EXAMPLES)
    def test_spec_examples_have_mass(self, seq):
        assert pt.parse(seq).mass() > 0

    def test_formula_component_equals_named_monosaccharide(self):
        # {C8H13N1O5} is exactly a HexNAc.
        formula = pt.parse("N[Glycan:{C8H13N1O5}1]")
        named = pt.parse("N[Glycan:HexNAc1]")
        assert abs(formula.mass() - named.mass()) < 1e-6
        assert dict(formula.comp(charge=0)) == dict(named.comp(charge=0))

    def test_formula_component_composition_resolves(self):
        # A formula component contributes an elemental composition (unlike a bare mass).
        comp = pt.parse("N[Glycan:{C8H13N1O5}2]").comp(charge=0)
        assert sum(comp.values()) > 0

    def test_charged_formula_component_surfaces_charge(self):
        annot = pt.parse("SEQUEN[Glycan:{C8H13N1O5Na1:z+1}1Hex2]CE")
        # The +1 from the charged formula component is reflected in the peptide's composition path.
        assert annot.mass() > 0
        assert annot.serialize() == "SEQUEN[Glycan:{C8H13N1O5Na1:z+1}1Hex2]CE"

    def test_bare_mass_component_adds_its_mass(self):
        base = pt.parse("N[Glycan:Hex2]").mass()
        with_mass = pt.parse("N[Glycan:{+500.0}1Hex2]").mass()
        assert abs((with_mass - base) - 500.0) < 1e-6

    def test_bare_mass_occurrence_multiplies(self):
        base = pt.parse("N").mass()
        assert abs((pt.parse("N[Glycan:{+500.0}2]").mass() - base) - 1000.0) < 1e-6

    def test_integer_mass_component_normalizes_and_round_trips(self):
        component = pt.GlycanComponent(monosaccharide=203, occurance=1)
        assert component.is_mass
        assert component.get_mass() == 203.0
        assert component.serialize() == "{+203.0}"
        assert pt.GlycanComponent.from_string(component.serialize()).get_mass() == 203.0

    def test_repeated_bare_mass_mod_scales_delta_mass(self):
        mods = pt.Mods(pt.ModType.INTERNAL, {"Glycan:{+203.079}": 2})
        composition, delta_mass, charge = mods.get_composition_with_delta_mass_charge()
        assert not composition
        assert delta_mass == pytest.approx(406.158)
        assert charge == 0

    def test_bare_mass_component_has_no_composition(self):
        # Consistent with any bare-mass modification: comp() cannot resolve a pure delta mass.
        with pytest.raises(ValueError):
            pt.parse("N[Glycan:{+203.079}1Hex2]").comp()

    def test_mixed_mass_and_monosaccharide_mass_is_additive(self):
        combined = pt.parse("N[Glycan:{+203.079}1Hex2]").mass()
        parts = pt.parse("N[Glycan:Hex2]").mass() + 203.079
        assert abs(combined - parts) < 1e-6

    @pytest.mark.parametrize("value", ["{C8H13N1O5}2", "{C8H13N1O5Na1:z+1}1", "{+203.079}1", "HexNAc4"])
    def test_component_reserializes_without_formula_prefix(self, value):
        # The real serializer (not the raw-string echo) must not emit a 'Formula:' prefix
        # inside the braces, and must re-parse to the same mass.
        from peptacular.proforma_components.parsers import parse_glycan_component

        gc = parse_glycan_component(value)
        out = gc.serialize()
        assert "Formula:" not in out
        assert abs(parse_glycan_component(out).get_mass() - gc.get_mass()) < 1e-9

    def test_glycan_tag_reserialization_preserves_mass(self):
        gt = pt.GlycanTag.from_string("Glycan:{C8H13N1O5}1Hex2")
        reparsed = pt.GlycanTag.from_string("Glycan:" + gt.serialize().split(":", 1)[1])
        assert abs(reparsed.get_mass() - gt.get_mass()) < 1e-9

    def test_whitespace_before_curly_component_count(self):
        spaced = pt.GlycanTag.from_string("Glycan:{C8H13N1O5} 2")
        compact = pt.GlycanTag.from_string("Glycan:{C8H13N1O5}2")
        assert spaced == compact

    @pytest.mark.parametrize(
        ("value", "message"),
        [
            ("Glycan:{}", "Empty '{}' glycan component"),
            ("Glycan:{+not-a-mass}", "Invalid glycan mass component"),
            ("Glycan:{C8H13N1O5", "Unclosed '{' in glycan composition"),
        ],
    )
    def test_invalid_curly_components_raise_clear_errors(self, value, message):
        with pytest.raises(ValueError, check=lambda e: message in str(e)):
            pt.GlycanTag.from_string(value)

    @pytest.mark.parametrize("value", ["+nan", "-nan", "+inf", "-inf", "1e309"])
    def test_non_finite_or_overflowing_mass_is_rejected(self, value):
        with pytest.raises(ValueError, match="Invalid glycan mass component"):
            pt.GlycanTag.from_string(f"Glycan:{{{value}}}")

    @pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
    def test_programmatic_non_finite_mass_is_rejected(self, value):
        with pytest.raises(ValueError, match="must be finite"):
            pt.GlycanComponent(monosaccharide=value, occurance=1)

"""
Tests for parsing ambiguous modifications and cross-linkers.
"""

import pytest

import peptacular as pt


class TestModificationAmbiguousPrimary:
    """Tests for pt.ModificationAmbiguousPrimary.from_string() function"""

    def test_simple_ambiguous_primary(self):
        """Test parsing simple ambiguous primary modification"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1")
        assert isinstance(result, pt.ModificationAmbiguousPrimary)
        assert result.label == "g1"  # Label includes everything from # part
        assert len(result.tags) == 1  # Phospho was in the label part, not tags

    def test_ambiguous_primary_with_score(self):
        """Test ambiguous primary with score"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1(0.8)")
        assert result.label == "g1"  # Label is extracted before score
        assert result.score == 0.8

    def test_ambiguous_primary_with_position(self):
        """Test ambiguous primary with position constraints"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1|Position:S,T,Y")
        assert result.label == "g1"
        assert result.position is not None
        assert len(result.position) == 3

    def test_ambiguous_primary_with_limit(self):
        """Test ambiguous primary with limit"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1|Limit:2")
        assert result.label == "g1"
        assert result.limit == 2

    def test_ambiguous_primary_with_comkp(self):
        """Test ambiguous primary with CoMKP flag"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1|CoMKP")
        assert result.label == "g1"
        assert result.comkp is True
        assert result.comup is False

    def test_ambiguous_primary_with_comup(self):
        """Test ambiguous primary with CoMUP flag"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1|CoMUP")
        assert result.label == "g1"
        assert result.comup is True
        assert result.comkp is False

    def test_ambiguous_primary_full_spec(self):
        """Test ambiguous primary with all features"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1(0.9)|Position:S,T|Limit:1|CoMKP")
        assert result.label == "g1"
        assert result.score == 0.9
        assert result.position is not None
        assert len(result.position) == 2
        assert result.limit == 1
        assert result.comkp is True

    def test_ambiguous_primary_with_accession(self):
        """Test ambiguous primary with accession tag"""
        result = pt.ModificationAmbiguousPrimary.from_string("UNIMOD:21#g1")
        assert result.label == "g1"  # Label includes the accession part
        assert len(result.tags) == 1

    def test_ambiguous_primary_multiple_tags(self):
        """Test ambiguous primary with multiple tags"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho|UNIMOD:21|+79.966#g1")
        # The last part with # is the label, previous parts are tags
        assert result.label == "g1"  # Last part containing # becomes label
        assert len(result.tags) == 3  # Phospho and UNIMOD:21 are tags

    def test_missing_label_raises_error(self):
        """Test that missing label raises ValueError"""
        with pytest.raises(ValueError):
            pt.ModificationAmbiguousPrimary.from_string("Phospho")

    def test_string_representation(self):
        """Test string representation"""
        result = pt.ModificationAmbiguousPrimary.from_string("Phospho#g1(0.8)")
        result_str = str(result)
        assert "Phospho" in result_str
        assert "#g1" in result_str
        assert "(0.8)" in result_str


class TestModificationAmbiguousSecondary:
    """Tests for pt.ModificationAmbiguousSecondary.from_string() function"""

    def test_simple_ambiguous_secondary(self):
        """Test parsing simple ambiguous secondary modification"""
        result = pt.ModificationAmbiguousSecondary.from_string("#g1")
        assert isinstance(result, pt.ModificationAmbiguousSecondary)
        assert result.label == "g1"
        assert result.score is None

    def test_ambiguous_secondary_with_score(self):
        """Test ambiguous secondary with score"""
        result = pt.ModificationAmbiguousSecondary.from_string("#g1(0.6)")
        assert result.label == "g1"
        assert result.score == 0.6

    def test_missing_hash_raises_error(self):
        """Test that missing # raises ValueError"""
        with pytest.raises(ValueError):
            pt.ModificationAmbiguousSecondary.from_string("g1")

    def test_various_label_names(self):
        """Test various label names"""
        labels = ["g1", "group1", "site_a", "123"]
        for label in labels:
            result = pt.ModificationAmbiguousSecondary.from_string(f"#{label}")
            assert result.label == label

    def test_string_representation(self):
        """Test string representation"""
        result = pt.ModificationAmbiguousSecondary.from_string("#g1")
        assert str(result) == "#g1"

    def test_string_representation_with_score(self):
        """Test string representation with score"""
        result = pt.ModificationAmbiguousSecondary.from_string("#g1(0.5)")
        assert str(result) == "#g1(0.5)"


class TestModificationCrossLinker:
    """Tests for pt.ModificationCrossLinker.from_string() function"""

    def test_cross_linker_primary_with_tags(self):
        """Test parsing cross-linker primary definition"""
        result = pt.ModificationCrossLinker.from_string("XLMOD:02001#XL1")
        assert isinstance(result, pt.ModificationCrossLinker)
        assert result.label == "1"  # XL prefix removed
        assert result.tags is not None
        assert len(result.tags) == 1

    def test_cross_linker_secondary_reference(self):
        """Test parsing cross-linker secondary reference"""
        result = pt.ModificationCrossLinker.from_string("#XL1")
        assert result.label == "1"  # XL prefix removed
        assert result.tags is None or len(result.tags) == 0

    def test_cross_linker_branch(self):
        """Test parsing branch cross-linker"""
        result = pt.ModificationCrossLinker.from_string("MOD:00093#BRANCH")
        assert result.label == "BRANCH"
        assert result.tags is not None

    def test_cross_linker_branch_secondary(self):
        """Test parsing branch secondary reference"""
        result = pt.ModificationCrossLinker.from_string("#BRANCH")
        assert result.label == "BRANCH"

    def test_cross_linker_with_mass(self):
        """Test cross-linker with mass tag"""
        result = pt.ModificationCrossLinker.from_string("+138.068#XL1")
        assert result.label == "1"
        assert result.tags is not None
        assert isinstance(result.tags[0], pt.TagMass)

    def test_cross_linker_multiple_tags(self):
        """Test cross-linker with multiple tags"""
        result = pt.ModificationCrossLinker.from_string("XLMOD:02001|DSS|+138.068#XL1")
        assert result.label == "1"
        assert result.tags is not None
        assert len(result.tags) == 3

    def test_missing_hash_raises_error(self):
        """Test that missing # raises ValueError"""
        with pytest.raises(ValueError):
            pt.ModificationCrossLinker.from_string("XLMOD:02001")

    def test_xl_prefix_removed_from_label(self):
        """Test that XL prefix is removed from label"""
        result = pt.ModificationCrossLinker.from_string("#XL2")
        assert result.label == "2"

        result2 = pt.ModificationCrossLinker.from_string("#XL123")
        assert result2.label == "123"

    def test_non_xl_label_preserved(self):
        """Test that non-XL labels are preserved"""
        result = pt.ModificationCrossLinker.from_string("#BRANCH")
        assert result.label == "BRANCH"

        result2 = pt.ModificationCrossLinker.from_string("#custom")
        assert result2.label == "custom"

    def test_string_representation_primary(self):
        """Test string representation of primary cross-linker"""
        result = pt.ModificationCrossLinker.from_string("XLMOD:02001#XL1")
        result_str = str(result)
        # CV enum might normalize XLMOD to XL-MOD
        assert "XLMOD:02001" in result_str or "XL-MOD:02001" in result_str
        assert "#XL1" in result_str

    def test_string_representation_secondary(self):
        """Test string representation of secondary cross-linker"""
        result = pt.ModificationCrossLinker.from_string("#XL1")
        assert str(result) == "#XL1"


class TestParseModificationDispatch:
    """Tests for the ``parse_modification`` router that picks the modification type.

    Regression: cross-link *definitions* like ``XLMOD:02001#XL1`` (which do not start
    with ``#xl``) were routed to the ambiguous-primary parser, and any modification
    containing a ``|`` was routed to the cross-linker parser.
    """

    @pytest.mark.parametrize(
        "text, expected",
        [
            ("XLMOD:02001#XL1", "ModificationCrossLinker"),
            ("#XL1", "ModificationCrossLinker"),
            ("XLMOD:02001|DSS|+138.068#XL1", "ModificationCrossLinker"),
            ("MOD:00093#BRANCH", "ModificationCrossLinker"),
            ("#BRANCH", "ModificationCrossLinker"),
            ("Oxidation#1", "ModificationAmbiguousPrimary"),
            ("Phospho|UNIMOD:21|+79.966#g1", "ModificationAmbiguousPrimary"),
            ("Phospho#g1|Position:S,T,Y", "ModificationAmbiguousPrimary"),
            ("#g1", "ModificationAmbiguousSecondary"),
            ("Oxidation", "ModificationTags"),
        ],
    )
    def test_dispatch(self, text, expected):
        from peptacular.proforma_components.parsers import parse_modification

        assert type(parse_modification(text)).__name__ == expected

"""Tests for parsing ProForma name annotations"""

import pytest

import peptacular as pt


class TestNameParsing:
    """Test parsing of (>name), (>>name), and (>>>name) annotations"""

    def test_simple_peptidoform_name(self):
        """Test parsing a peptidoform with single > name"""
        parser = pt.annotation.parser.ProFormaParser("(>myPeptide)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, connection = result[0]
        assert parsed.peptide_name == "myPeptide"
        assert parsed.unmod_sequence == "PEPTIDE"
        assert connection is None

    def test_invalid_name_format(self):
        """Test that invalid name formats raise errors"""
        # Missing '>'
        with pytest.raises(ValueError):
            parser = pt.annotation.parser.ProFormaParser("(myPeptide)PEPTIDE")
            list(parser.parse())

        # Too many '>'
        with pytest.raises(ValueError):
            parser = pt.annotation.parser.ProFormaParser("(>>>>myPeptide)PEPTIDE")
            list(parser.parse())

        # out of order names
        with pytest.raises(ValueError):
            parser = pt.annotation.parser.ProFormaParser("(>compound)(>>>ion)(>>peptide)PEPTIDE")
            list(parser.parse())

        # multiple names of same type
        with pytest.raises(ValueError):
            parser = pt.annotation.parser.ProFormaParser("(>peptide1)(>peptide2)PEPTIDE")
            list(parser.parse())

    def test_peptidoform_ion_name(self):
        """Test parsing a peptidoform ion with >> name"""
        parser = pt.annotation.parser.ProFormaParser("(>>myIon)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, connection = result[0]
        assert parsed.ion_name == "myIon"
        assert parsed.unmod_sequence == "PEPTIDE"
        assert connection is None

    def test_compound_peptidoform_name(self):
        """Test parsing a compound peptidoform with >>> name"""
        parser = pt.annotation.parser.ProFormaParser("(>>>myCompound)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, connection = result[0]
        assert parsed.compound_name == "myCompound"
        assert parsed.unmod_sequence == "PEPTIDE"
        assert connection is None

    def test_name_with_modifications(self):
        """Test name combined with modifications"""
        parser = pt.annotation.parser.ProFormaParser("(>named)[Acetyl]-PEM[Oxidation]TIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "named"
        assert parsed.unmod_sequence == "PEMTIDE"
        assert parsed.nterm_mods is not None
        assert "Acetyl" in parsed.nterm_mods
        assert parsed.internal_mods is not None
        assert 2 in parsed.internal_mods
        assert "Oxidation" in parsed.internal_mods[2]

    def test_name_with_labile_mods(self):
        """Test name with labile modifications"""
        parser = pt.annotation.parser.ProFormaParser("(>glyco){Glycan:Hex}PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "glyco"
        assert parsed.labile_mods is not None
        assert "Glycan:Hex" in parsed.labile_mods

    def test_name_with_unknown_mods(self):
        """Test name with unknown position modifications"""
        parser = pt.annotation.parser.ProFormaParser("(>unknown)[Phospho]?PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "unknown"
        assert parsed.unknown_mods is not None
        assert "Phospho" in parsed.unknown_mods

    def test_name_with_global_mods(self):
        """Test name with global modifications"""
        parser = pt.annotation.parser.ProFormaParser("<13C>(>isotope)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "isotope"
        assert parsed.global_mods is not None
        assert "13C" in parsed.global_mods

    def test_name_with_charge(self):
        """Test name with charge state"""
        parser = pt.annotation.parser.ProFormaParser("(>charged)PEPTIDE/2")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "charged"
        assert parsed.charge == 2

    def test_name_with_special_characters(self):
        """Test name with various characters"""
        parser = pt.annotation.parser.ProFormaParser("(>my_peptide-123)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "my_peptide-123"

    def test_name_with_spaces(self):
        """Test name with spaces"""
        parser = pt.annotation.parser.ProFormaParser("(>my peptide name)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "my peptide name"

    def test_name_with_balanced_parentheses(self):
        """A name may contain balanced parentheses; the outer ')' still closes the name."""
        parser = pt.annotation.parser.ProFormaParser("(>my (special) peptide)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "my (special) peptide"
        assert parsed.unmod_sequence == "PEPTIDE"

    def test_name_with_unbalanced_parenthesis(self):
        """A stray unbalanced '(' in a name stays parseable (falls back to the first ')');
        regression: balanced-paren scanning previously rejected these outright."""
        parser = pt.annotation.parser.ProFormaParser("(>a(b)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == "a(b"
        assert parsed.unmod_sequence == "PEPTIDE"

    def test_crosslinked_with_names(self):
        """Test crosslinked peptides with names"""
        parser = pt.annotation.parser.ProFormaParser("(>pep1)PEPTK[#XL1]IDE//(>pep2)SEQK[#XL1]UENCE")

        # Capture data from each segment as we iterate
        segments: list[dict[str, object]] = []
        for p, conn in parser.parse():
            # Copy the data we need from the parser state
            segments.append(
                {
                    "name": p.peptide_name,
                    "sequence": p.unmod_sequence,
                    "connection": conn,
                }
            )

        assert len(segments) == 2

        assert segments[0]["name"] == "pep1"
        assert segments[0]["sequence"] == "PEPTKIDE"
        assert segments[0]["connection"] is True  # Crosslinked

        assert segments[1]["name"] == "pep2"
        assert segments[1]["sequence"] == "SEQKUENCE"
        assert segments[1]["connection"] is None

    def test_chimeric_with_names(self):
        """Test chimeric peptides with names"""
        parser = pt.annotation.parser.ProFormaParser("(>first)PEPTIDE+(>second)SEQUENCE")

        # Capture data from each segment as we iterate
        segments: list[dict[str, object]] = []
        for p, conn in parser.parse():
            # Copy the data we need from the parser state
            segments.append(
                {
                    "name": p.peptide_name,
                    "sequence": p.unmod_sequence,
                    "connection": conn,
                }
            )

        assert len(segments) == 2

        assert segments[0]["name"] == "first"
        assert segments[0]["sequence"] == "PEPTIDE"
        assert segments[0]["connection"] is False  # Chimeric

        assert segments[1]["name"] == "second"
        assert segments[1]["sequence"] == "SEQUENCE"
        assert segments[1]["connection"] is None

    def test_no_name(self):
        """Test that sequences without names work as before"""
        parser = pt.annotation.parser.ProFormaParser("PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name is None
        assert parsed.unmod_sequence == "PEPTIDE"

    def test_annotation_parse_with_name(self):
        """Test ProFormaAnnotation.parse() with name"""
        annot = pt.ProFormaAnnotation.parse("(>myPeptide)PEPTIDE")

        assert annot.peptide_name == "myPeptide"
        assert annot.sequence == "PEPTIDE"

    def test_annotation_parse_without_name(self):
        """Test ProFormaAnnotation.parse() without name"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE")

        assert annot.peptide_name == ""
        assert annot.sequence == "PEPTIDE"

    def test_invalid_arrow_count(self):
        """Test that invalid arrow counts raise errors"""
        # Zero arrows
        with pytest.raises(ValueError):
            parser = pt.annotation.parser.ProFormaParser("()PEPTIDE")
            list(parser.parse())

        # Four arrows
        with pytest.raises(ValueError):
            parser = pt.annotation.parser.ProFormaParser("(>>>>name)PEPTIDE")
            list(parser.parse())

    def test_unclosed_name(self):
        """Test that unclosed name parenthesis raises error"""
        with pytest.raises(ValueError):
            parser = pt.annotation.parser.ProFormaParser("(>myName")
            list(parser.parse())

    def test_empty_name(self):
        """Test that empty names are allowed"""
        parser = pt.annotation.parser.ProFormaParser("(>)PEPTIDE")
        result = list(parser.parse())

        assert len(result) == 1
        parsed, _ = result[0]
        assert parsed.peptide_name == ""
        assert parsed.unmod_sequence == "PEPTIDE"

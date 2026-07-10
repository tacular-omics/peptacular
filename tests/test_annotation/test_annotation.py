"""Tests for the ProFormaAnnotation class"""

import pytest

import peptacular as pt


class TestAnnotationConstruction:
    """Test creating ProFormaAnnotation objects directly"""

    def test_empty_annotation(self):
        """Test creating an empty annotation"""
        annot = pt.ProFormaAnnotation()
        assert annot.sequence == ""
        assert annot.peptide_name == ""
        assert annot._isotope_mods is None
        assert annot._static_mods is None
        assert annot._labile_mods is None
        assert annot._unknown_mods is None
        assert annot._nterm_mods is None
        assert annot._cterm_mods is None
        assert annot._internal_mods is None
        assert annot._intervals is None
        assert annot.charge is None

    def test_annotation_with_sequence(self):
        """Test creating annotation with just a sequence"""
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot.peptide_name == ""

    def test_annotation_with_name(self):
        """Test creating annotation with name"""
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", peptide_name="myPeptide")
        assert annot.sequence == "PEPTIDE"
        assert annot.peptide_name == "myPeptide"

    def test_annotation_with_charge(self):
        """Test creating annotation with charge"""
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE", charge=2)
        assert annot.sequence == "PEPTIDE"
        assert annot.charge == 2


class TestAnnotationParse:
    """Test parsing ProForma strings into annotations"""

    def test_parse_simple_sequence(self):
        """Test parsing simple unmodified sequence"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot.peptide_name == ""
        assert annot._nterm_mods is None
        assert annot._cterm_mods is None
        assert annot._internal_mods is None

    def test_parse_simple_sequence_with_validation(self):
        """Test parsing simple unmodified sequence with validation"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot.peptide_name == ""

    def test_parse_with_name(self):
        """Test parsing sequence with name"""
        annot = pt.ProFormaAnnotation.parse("(>myPeptide)PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot.peptide_name == "myPeptide"

    def test_parse_with_name_with_validation(self):
        """Test parsing sequence with name with validation"""
        annot = pt.ProFormaAnnotation.parse("(>myPeptide)PEPTIDE", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot.peptide_name == "myPeptide"

    def test_parse_with_nterm_mod(self):
        """Test parsing with N-terminal modification"""
        annot = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot._nterm_mods is not None
        assert "Acetyl" in annot._nterm_mods
        assert annot._nterm_mods["Acetyl"] == 1

    def test_parse_with_nterm_mod_with_validation(self):
        """Test parsing with N-terminal modification with validation"""
        annot = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot._nterm_mods is not None
        assert "Acetyl" in annot._nterm_mods

    def test_parse_with_cterm_mod(self):
        """Test parsing with C-terminal modification"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE-[Amidated]")
        assert annot.sequence == "PEPTIDE"
        assert annot._cterm_mods is not None
        assert "Amidated" in annot._cterm_mods
        assert annot._cterm_mods["Amidated"] == 1

    def test_parse_with_cterm_mod_with_validation(self):
        """Test parsing with C-terminal modification with validation"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE-[Amidated]", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot._cterm_mods is not None
        assert "Amidated" in annot._cterm_mods

    def test_parse_with_internal_mod(self):
        """Test parsing with internal modification"""
        annot = pt.ProFormaAnnotation.parse("PEM[Oxidation]TIDE")
        assert annot.sequence == "PEMTIDE"
        assert annot._internal_mods is not None
        assert 2 in annot._internal_mods
        assert "Oxidation" in annot._internal_mods[2]
        assert annot._internal_mods[2]["Oxidation"] == 1

    def test_parse_with_internal_mod_with_validation(self):
        """Test parsing with internal modification with validation"""
        annot = pt.ProFormaAnnotation.parse("PEM[Oxidation]TIDE", validate=True)
        assert annot.sequence == "PEMTIDE"
        assert annot._internal_mods is not None
        assert "Oxidation" in annot._internal_mods[2]

    def test_parse_with_multiple_internal_mods(self):
        """Test parsing with multiple internal modifications"""
        annot = pt.ProFormaAnnotation.parse("PEM[Oxidation]TIS[Phospho]DE")
        assert annot.sequence == "PEMTISDE"
        assert annot._internal_mods is not None
        assert 2 in annot._internal_mods
        assert "Oxidation" in annot._internal_mods[2]
        assert 5 in annot._internal_mods
        assert "Phospho" in annot._internal_mods[5]

    def test_parse_with_multiple_internal_mods_with_validation(self):
        """Test parsing with multiple internal modifications with validation"""
        annot = pt.ProFormaAnnotation.parse("PEM[Oxidation]TIS[Phospho]DE", validate=True)
        assert annot.sequence == "PEMTISDE"
        assert annot._internal_mods is not None
        assert "Oxidation" in annot._internal_mods[2]
        assert "Phospho" in annot._internal_mods[5]

    def test_parse_with_labile_mod(self):
        """Test parsing with labile modification"""
        annot = pt.ProFormaAnnotation.parse("{Glycan:Hex}PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot._labile_mods is not None
        assert "Glycan:Hex" in annot._labile_mods
        assert annot._labile_mods["Glycan:Hex"] == 1

    def test_parse_with_labile_mod_with_validation(self):
        """Test parsing with labile modification with validation"""
        annot = pt.ProFormaAnnotation.parse("{Glycan:Hex}PEPTIDE", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot._labile_mods is not None
        assert "Glycan:Hex" in annot._labile_mods

    def test_parse_with_unknown_mod(self):
        """Test parsing with unknown position modification"""
        annot = pt.ProFormaAnnotation.parse("[Phospho]?PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot._unknown_mods is not None
        assert "Phospho" in annot._unknown_mods
        assert annot._unknown_mods["Phospho"] == 1

    def test_parse_with_unknown_mod_with_validation(self):
        """Test parsing with unknown position modification with validation"""
        annot = pt.ProFormaAnnotation.parse("[Phospho]?PEPTIDE", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot._unknown_mods is not None
        assert "Phospho" in annot._unknown_mods

    def test_parse_with_isotope_mod(self):
        """Test parsing with isotope modification"""
        annot = pt.ProFormaAnnotation.parse("<13C>PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot._isotope_mods is not None
        assert "13C" in annot._isotope_mods
        assert annot._isotope_mods["13C"] == 1

    def test_parse_with_isotope_mod_with_validation(self):
        """Test parsing with isotope modification with validation"""
        annot = pt.ProFormaAnnotation.parse("<13C>PEPTIDE", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot._isotope_mods is not None
        assert "13C" in annot._isotope_mods

    def test_parse_with_static_mod(self):
        """Test parsing with static modification"""
        annot = pt.ProFormaAnnotation.parse("<[Carbamidomethyl@C]>PEPTIDE")
        assert annot.sequence == "PEPTIDE"
        assert annot._static_mods is not None
        assert "[Carbamidomethyl@C]" in annot._static_mods
        assert annot._static_mods["[Carbamidomethyl@C]"] == 1

    def test_parse_with_static_mod_with_validation(self):
        """Test parsing with static modification with validation"""
        annot = pt.ProFormaAnnotation.parse("<[Oxidation]@M>PEPTIDE", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert len(annot.static_mods) == 1
        assert "[Oxidation]@M" in annot.static_mods

    def test_parse_with_charge(self):
        """Test parsing with charge state"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        assert annot.sequence == "PEPTIDE"
        assert annot.charge == 2

    def test_parse_with_charge_with_validation(self):
        """Test parsing with charge state with validation"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE/2", validate=True)
        assert annot.sequence == "PEPTIDE"
        assert annot.charge == 2

    def test_parse_with_charge_adducts(self):
        """Test parsing with charge adducts (adduct only, no charge number)"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE/[H:z+1]")
        assert annot.sequence == "PEPTIDE"
        assert annot.charge_state == 1
        adducts = annot.charge_adducts
        assert len(adducts) == 1
        assert adducts.serialize() == "[H:z+1]"

    def test_parse_with_adduct_only(self):
        """Test parsing with adduct only (no charge number)"""
        annot = pt.ProFormaAnnotation.parse("SEQUEN/[Na:z+1,H:z+1]")
        assert annot.sequence == "SEQUEN"
        assert annot.charge_state == 2
        adducts = annot.charge_adducts
        assert len(adducts) == 2

    def test_parse_with_multiple_adducts(self):
        """Test parsing with multiple adduct modifications - must be comma-separated in single bracket"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE/[H:z+1,Na:z+1]")
        assert annot.sequence == "PEPTIDE"
        assert annot.charge_state == 2
        adducts = annot.charge_adducts
        assert len(adducts) == 2

    def test_parse_with_adduct_multiplier(self):
        """Test parsing with adduct multiplier (^n)"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE/[Na:z+1^2]")
        assert annot.sequence == "PEPTIDE"
        assert annot.charge_state == 2
        adducts = annot.charge_adducts
        assert len(adducts) == 1

    def test_parse_with_adduct_comma_and_multiplier(self):
        """Test parsing with comma-separated adducts and multiplier"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE/[Na:z+1^2,H:z+1]")
        assert annot.sequence == "PEPTIDE"
        assert annot.charge_state == 3
        adducts = annot.charge_adducts
        assert len(adducts) == 2

    def test_parse_with_charge_and_adduct_raises_error(self):
        """Test that having both charge and adduct raises an error"""
        with pytest.raises(ValueError):
            pt.ProFormaAnnotation.parse("PEPTIDE/2[+H]")

    def test_parse_with_bracket_multiplier(self):
        """Test parsing with comma-separated adducts and multiplier"""
        with pytest.raises(ValueError):
            _ = pt.ProFormaAnnotation.parse("PEPTIDE/[Na:z+1,H:z+1]^2")

    def test_parse_with_multi_bracket(self):
        """Test parsing with comma-separated adducts and multiplier"""
        with pytest.raises(ValueError):
            _ = pt.ProFormaAnnotation.parse("PEPTIDE/[Na:z+1][H:z+1]")

    def test_parse_with_invalid_amino_acid_raises_error(self):
        """Test that parsing with invalid amino acid raises error with validation"""
        with pytest.raises(ValueError):
            pt.ProFormaAnnotation.parse("PEPT@IDE", validate=True)

    def test_parse_rejects_chimeric(self):
        """Test that parsing rejects chimeric peptides"""
        with pytest.raises(
            ValueError,
        ):
            pt.ProFormaAnnotation.parse("PEPTIDE+SEQUENCE")

    def test_parse_crosslinked_returns_multi(self):
        """Cross-linked ('//') peptides parse into a MultiProFormaAnnotation."""
        ion = pt.ProFormaAnnotation.parse("PEPTK[XLMOD:02001#XL1]IDE//SEQK[#XL1]")
        assert isinstance(ion, pt.MultiProFormaAnnotation)
        assert len(ion) == 2
        # Round-trips back to the original string.
        assert ion.serialize() == "PEPTK[XLMOD:02001#XL1]IDE//SEQK[#XL1]"

    def test_round_trip_simple_peptide(self):
        """Test parse -> serialize round trip for simple peptide"""
        original = "PEPTIDE"
        annot = pt.ProFormaAnnotation.parse(original)
        result = annot.serialize()
        assert result == original

    def test_round_trip_with_modifications(self):
        """Test parse -> serialize round trip with modifications"""
        original = "[Acetyl]-PEM[Oxidation]TIDE-[Amidated]"
        annot = pt.ProFormaAnnotation.parse(original)
        result = annot.serialize()
        assert result == original

    def test_round_trip_with_name(self):
        """Test parse -> serialize round trip with name"""
        original = "(>myPeptide)PEPTIDE"
        annot = pt.ProFormaAnnotation.parse(original)
        result = annot.serialize()
        assert result == original

    def test_round_trip_with_charge(self):
        """Test parse -> serialize round trip with charge"""
        original = "PEPTIDE/2"
        annot = pt.ProFormaAnnotation.parse(original)
        result = annot.serialize()
        assert result == original

    def test_round_trip_complex(self):
        """Test parse -> serialize round trip with complex annotation"""
        original = "(>>>complex)<13C>[Acetyl]-PEM[Oxidation]TIS[Phospho]DE-[Amidated]/2"
        annot = pt.ProFormaAnnotation.parse(original)
        result = annot.serialize()
        assert result == original

    def test_round_trip_labile(self):
        """Test parse -> serialize round trip with labile modification"""
        original = "{Glycan:Hex}PEPTIDE"
        annot = pt.ProFormaAnnotation.parse(original)
        result = annot.serialize()
        assert result == original

    def test_round_trip_unknown(self):
        """Test parse -> serialize round trip with unknown modification"""
        original = "[Phospho]?PEPTIDE"
        annot = pt.ProFormaAnnotation.parse(original)
        result = annot.serialize()
        assert result == original

    def test_serialize_empty_annotation(self):
        """Test serializing empty annotation"""
        annot = pt.ProFormaAnnotation()
        result = annot.serialize()
        assert result == ""

    def test_serialize_only_name(self):
        """Test serializing annotation with only name (no sequence)"""
        annot = pt.ProFormaAnnotation(peptide_name="test")
        result = annot.serialize()
        assert result == "(>test)"

    def test_sequecne_regions(self):
        """Test getting sequence regions"""
        annot = pt.ProFormaAnnotation.parse("PE(PTID)[Oxidation]E")
        regions: tuple[pt.SequenceRegion, ...] = annot.sequence_regions
        assert len(regions) == 1
        assert str(regions[0]) == "(PTID)[Oxidation]"

        annot = pt.ProFormaAnnotation.parse("PE(?PTID)[Oxidation]E")
        regions: tuple[pt.SequenceRegion, ...] = annot.sequence_regions
        assert len(regions) == 1
        assert str(regions[0]) == "(?PTID)[Oxidation]"

    def test_sequence_elements(self):
        """Test getting sequence elements"""
        annot = pt.ProFormaAnnotation.parse("PEPTIDE")
        elements: tuple[pt.SequenceElement, ...] = annot.sequence_elements
        assert len(elements) == 7
        assert str(elements[0]) == "P"

    def test_sequence_elements_with_regions(self):
        """Test getting sequence elements"""
        annot = pt.ProFormaAnnotation.parse("P(?EPT)[oxidation]IDE")
        elements: tuple[pt.SequenceElement, ...] = annot.sequence_elements
        assert len(elements) == 7
        assert str(elements[0]) == "P"

    def test_sequence_elements_and_regions(self):
        """Test getting sequence elements and regions with modifications"""
        annot = pt.ProFormaAnnotation.parse("P(?E[Phospho]PT)[Oxidation]IDE")
        elems: tuple[pt.SequenceElement | pt.SequenceRegion, ...] = annot.sequence_elements_and_regions
        assert len(elems) == 5
        assert str(elems[0]) == "P"
        assert str(elems[1]) == "(?E[Phospho]PT)[Oxidation]"
        assert str(elems[2]) == "I"


if __name__ == "__main__":
    pytest.main([__file__])

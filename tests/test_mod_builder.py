import unittest
import warnings

import peptacular as pt
from peptacular.annotation.mod_builder import (
    apply_mods,
    apply_static_mods_infront,
    modify,
    ensure_single_static_mod,
    get_mod_index_from_aa,
    get_mod_index_from_regex,
    get_sites,
)


class TestGetModIndexFromAA(unittest.TestCase):
    def test_single_match_beginning(self):
        """Test single match at beginning of peptide."""
        result = get_mod_index_from_aa("PEPTIDE", "P")
        self.assertEqual(result, {0, 2})

    def test_single_match_end(self):
        """Test single match at end of peptide."""
        result = get_mod_index_from_aa("PEPTIDE", "E")
        self.assertEqual(result, {1, 6})

    def test_single_match_middle(self):
        """Test single match in middle of peptide."""
        result = get_mod_index_from_aa("PEPTIDE", "T")
        self.assertEqual(result, {3})

    def test_multiple_amino_acids(self):
        """Test multiple amino acids in mod pattern."""
        result = get_mod_index_from_aa("PEPTIDE", "PE")
        self.assertEqual(result, {0, 1, 2, 6})

    def test_no_match(self):
        """Test when no matches are found."""
        result = get_mod_index_from_aa("PEPTIDE", "X")
        self.assertEqual(result, set())

    def test_empty_peptide(self):
        """Test with empty peptide sequence."""
        result = get_mod_index_from_aa("", "P")
        self.assertEqual(result, set())

    def test_empty_mod_aa(self):
        """Test with empty modification amino acids."""
        result = get_mod_index_from_aa("PEPTIDE", "")
        self.assertEqual(result, set())

    def test_all_positions_match(self):
        """Test when all positions match."""
        result = get_mod_index_from_aa("PPPP", "P")
        self.assertEqual(result, {0, 1, 2, 3})

    def test_single_character_peptide(self):
        """Test with single character peptide."""
        result = get_mod_index_from_aa("P", "P")
        self.assertEqual(result, {0})


class TestGetModIndexFromRegexImplemented(unittest.TestCase):
    def test_single_character_regex_match(self):
        """Test regex matching single character."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?=P)")
        self.assertEqual(result, {0})  # Should match first P only

    def test_regex_match_at_end(self):
        """Test regex matching at end of sequence."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?<=E)$")
        self.assertEqual(result, {7})  # Should match after E at end

    def test_regex_match_at_beginning(self):
        """Test regex matching at beginning of sequence."""
        result = get_mod_index_from_regex("PEPTIDE", r"^(?=P)")
        self.assertEqual(result, {0})  # Should match P at start

    def test_regex_no_match(self):
        """Test regex with no matches."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?=X)")
        self.assertEqual(result, set())

    def test_regex_multiple_matches(self):
        """Test regex that could match multiple positions."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?=E)")
        self.assertEqual(result, {1})  # Should match first E only

    def test_regex_character_class(self):
        """Test regex with character class."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?=[PT])")
        self.assertEqual(result, {0})  # Should match first P

    def test_regex_lookahead_pattern(self):
        """Test regex with lookahead pattern."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?=PE)")
        self.assertEqual(result, {0})  # P followed by E

    def test_regex_lookbehind_pattern(self):
        """Test regex with lookbehind pattern."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?<=PE)")
        self.assertEqual(result, {2})  # After PE

    def test_regex_empty_sequence(self):
        """Test regex with empty sequence."""
        result = get_mod_index_from_regex("", r"(?=P)")
        self.assertEqual(result, set())

    def test_regex_case_sensitive(self):
        """Test regex case sensitivity."""
        result = get_mod_index_from_regex("peptide", r"(?=P)")
        self.assertEqual(result, set())  # No match due to case

    def test_regex_case_insensitive_flag(self):
        """Test regex with case insensitive flag."""
        result = get_mod_index_from_regex("peptide", r"(?i)(?=P)")
        self.assertEqual(result, {0})  # Should match with flag


class TestGetSites(unittest.TestCase):
    def setUp(self):
        self.phospho_mod = "Phospho"
        self.oxidation_mod = "Oxidation"

    def test_empty_mods(self):
        """Test with empty modifications dictionary."""
        result = get_sites("PEPTIDE", {}, is_regex=False)
        self.assertEqual(result, {})

    def test_single_mod_single_aa(self):
        """Test single modification on single amino acid."""
        mods = {"P": [self.phospho_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {0: [self.phospho_mod], 2: [self.phospho_mod]}
        self.assertEqual(result, expected)

    def test_single_mod_multiple_aa(self):
        """Test single modification on multiple amino acids."""
        mods = {"PE": [self.oxidation_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.oxidation_mod],
            1: [self.oxidation_mod],
            2: [self.oxidation_mod],
            6: [self.oxidation_mod],
        }
        self.assertEqual(result, expected)

    def test_multiple_mods_same_aa(self):
        """Test multiple modifications on same amino acid."""
        mods = {"P": [self.phospho_mod, self.oxidation_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.phospho_mod, self.oxidation_mod],
            2: [self.phospho_mod, self.oxidation_mod],
        }
        self.assertEqual(result, expected)

    def test_multiple_mod_types(self):
        """Test multiple different modification types."""
        mods = {"P": [self.phospho_mod], "E": [self.oxidation_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.phospho_mod],
            1: [self.oxidation_mod],
            2: [self.phospho_mod],
            6: [self.oxidation_mod],
        }
        self.assertEqual(result, expected)

    def test_numeric_modifications(self):
        """Test with numeric modification values."""
        mods = {"P": [123.456]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {0: [123.456], 2: [123.456]}
        self.assertEqual(result, expected)

    def test_string_modifications(self):
        """Test with string modification values."""
        mods = {"P": ["CustomMod"]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {0: ["CustomMod"], 2: ["CustomMod"]}
        self.assertEqual(result, expected)

    def test_mixed_modification_types(self):
        """Test with mixed modification value types."""
        mods = {"P": [self.phospho_mod, 123.456, "StringMod"]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.phospho_mod, 123.456, "StringMod"],
            2: [self.phospho_mod, 123.456, "StringMod"],
        }
        self.assertEqual(result, expected)

    def test_no_matching_sites(self):
        """Test when no sites match the modification pattern."""
        mods = {"X": [self.phospho_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        self.assertEqual(result, {})

    def test_empty_mod_list(self):
        """Test with empty modification list."""
        mods = {"P": []}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        self.assertEqual(result, {0: [], 2: []})


class TestEnsureSingleStaticMod(unittest.TestCase):
    def setUp(self):
        self.phospho_mod = "Phospho"
        self.oxidation_mod = "Oxidation"

    def test_single_mod_per_site_no_warning(self):
        """Test no warning when single mod per site."""
        mods = {0: [self.phospho_mod], 1: [self.oxidation_mod]}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ensure_single_static_mod(mods)
            self.assertEqual(len(w), 0)

    def test_multiple_mods_per_site_warning(self):
        """Test warning when multiple mods per site."""
        mods = {0: [self.phospho_mod, self.oxidation_mod]}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ensure_single_static_mod(mods)
            self.assertEqual(len(w), 1)
            self.assertIn("Multiple static modifications", str(w[0].message))

    def test_empty_mods_no_warning(self):
        """Test no warning with empty modifications."""
        mods = {}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ensure_single_static_mod(mods)
            self.assertEqual(len(w), 0)


class TestApplyMods(unittest.TestCase):
    def setUp(self):
        self.annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        self.phospho_mod = "Phospho"
        self.oxidation_mod = "Oxidation"
        self.acetyl_mod = "Acetyl"

    def test_apply_static_mod_infront(self):
        """Test applying only static modifications."""
        static_mods = {"P": [self.phospho_mod]}
        result = apply_static_mods_infront(self.annotation, internal_static=static_mods)
        self.assertEqual(result.serialize(), "<[Phospho]@P>PEPTIDE")

    def test_apply_static_mod_infront2(self):
        """Test applying only static modifications."""
        static_mods = {"PC": [self.phospho_mod]}
        result = apply_static_mods_infront(self.annotation, internal_static=static_mods)
        # When using multi-character patterns, they're applied to each character separately
        self.assertEqual(result.serialize(), "<[Phospho]@P>PEPTIDE")

    def test_apply_static_mod_infront3(self):
        """Test applying only static modifications."""
        static_mods = {"PE": [self.phospho_mod]}
        result = apply_static_mods_infront(self.annotation, internal_static=static_mods)
        # When using multi-character patterns, they're applied to each character separately
        self.assertEqual(result.serialize(), "<[Phospho]@P,E>PEPTIDE")

    def test_apply_internal_mods_only(self):
        """Test applying only internal modifications."""
        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=False)
        self.assertEqual(result.serialize(), "P[Phospho]EP[Phospho]TIDE")

    def test_apply_nterm_mods_only(self):
        """Test applying only N-terminal modifications."""
        nterm_mods = {"P": [self.acetyl_mod]}
        result = apply_mods(self.annotation, nterm=nterm_mods, inplace=False)
        self.assertEqual(result.serialize(), "[Acetyl]-PEPTIDE")

    def test_apply_cterm_mods_only(self):
        """Test applying only C-terminal modifications."""
        cterm_mods = {"E": [self.oxidation_mod]}
        result = apply_mods(self.annotation, cterm=cterm_mods, inplace=False)
        self.assertEqual(result.serialize(), "PEPTIDE-[Oxidation]")

    def test_apply_all_mod_types(self):
        """Test applying all modification types together."""
        nterm_mods = {"P": [self.acetyl_mod]}
        cterm_mods = {"E": [self.oxidation_mod]}
        internal_mods = {"T": [self.phospho_mod]}

        result = apply_mods(
            self.annotation,
            nterm=nterm_mods,
            cterm=cterm_mods,
            internal=internal_mods,
            inplace=False,
        )
        self.assertEqual(result.serialize(), "[Acetyl]-PEPT[Phospho]IDE-[Oxidation]")

    def test_inplace_true_returns_same_object(self):
        """Test inplace=True returns the same annotation object."""
        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=True)
        self.assertIs(result, self.annotation)

    def test_inplace_false_returns_copy(self):
        """Test inplace=False returns a different annotation object."""
        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=False)
        self.assertIsNot(result, self.annotation)
        self.assertEqual(self.annotation.serialize(), "PEPTIDE")  # Original unchanged

    def test_merge_strategy_update(self):
        """Test merge strategy 'update' clears existing mods."""
        # Pre-populate annotation with modifications
        self.annotation.set_internal_mods_at_index(0, {"Oxidation": 1})

        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=True)

        # Should only have phospho mod, oxidation should be cleared
        self.assertEqual(result.serialize(), "P[Oxidation][Phospho]EP[Phospho]TIDE")

    def test_merge_strategy_merge(self):
        """Test merge strategy 'merge' combines with existing mods."""
        # Pre-populate annotation with modifications
        self.annotation.set_internal_mods_at_index(0, {"Oxidation": 1})

        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=True)

        # Should have both oxidation and phospho mods
        if result._internal_mods is None:
            raise AssertionError("Internal mods should not be None")

        mod_dict = result._internal_mods[0]
        self.assertIn("Oxidation", mod_dict)
        self.assertIn("Phospho", mod_dict)

    def test_none_modifications(self):
        """Test with None modification parameters."""
        result = apply_mods(self.annotation, nterm=None, cterm=None, internal=None, inplace=False)
        self.assertEqual(result.serialize(), "PEPTIDE")

    def test_empty_modification_dicts(self):
        """Test with empty modification dictionaries."""
        result = apply_mods(self.annotation, nterm={}, cterm={}, internal={}, inplace=False)
        self.assertEqual(result.serialize(), "PEPTIDE")


class TestBuildMods(unittest.TestCase):
    def setUp(self):
        self.annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        self.phospho_mod = "Phospho"
        self.oxidation_mod = "Oxidation"
        self.acetyl_mod = "Acetyl"

    def test_no_modifications_single_result(self):
        """Test modify with no modifications returns single result."""
        results = list(modify(self.annotation, max_variable_mods=0))
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "PEPTIDE")

    def test_static_modifications_only(self):
        """Test modify with only static modifications."""
        results = list(
            modify(
                self.annotation,
                internal_static={"P": [self.phospho_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "P[Phospho]EP[Phospho]TIDE")

    def test_variable_modifications_only(self):
        """Test modify with only variable modifications."""
        results = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate: no mods, P@0, P@2
        self.assertEqual(len(results), 3)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("P[Phospho]EPTIDE", serialized)
        self.assertIn("PEP[Phospho]TIDE", serialized)

    def test_static_and_variable_modifications(self):
        """Test modify with both static and variable modifications."""
        results = list(
            modify(
                self.annotation,
                internal_static={"T": [self.phospho_mod]},
                internal_variable={"P": [self.oxidation_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate: static only, static + P@0, static + P@2
        self.assertEqual(len(results), 3)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPT[Phospho]IDE", serialized)
        self.assertIn("P[Oxidation]EPT[Phospho]IDE", serialized)
        self.assertIn("PEP[Oxidation]T[Phospho]IDE", serialized)

    def test_static_and_variable_modifications_infront(self):
        """Test modify with both static and variable modifications."""
        results = list(
            modify(
                self.annotation,
                internal_static={"T": [self.phospho_mod]},
                internal_variable={"P": [self.oxidation_mod]},
                max_variable_mods=1,
                use_static_notation=True,
            )
        )
        # Should generate: static only, static + P@0, static + P@2
        self.assertEqual(len(results), 3)
        serialized = [r.serialize() for r in results]
        self.assertIn("<[Phospho]@T>PEPTIDE", serialized)
        self.assertIn("<[Phospho]@T>P[Oxidation]EPTIDE", serialized)
        self.assertIn("<[Phospho]@T>PEP[Oxidation]TIDE", serialized)

    def test_nterm_static_modifications(self):
        """Test N-terminal static modifications."""
        results = list(
            modify(
                self.annotation,
                nterm_static={"P": [self.acetyl_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "[Acetyl]-PEPTIDE")

    def test_cterm_static_modifications(self):
        """Test C-terminal static modifications."""
        results = list(
            modify(
                self.annotation,
                cterm_static={"E": [self.oxidation_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "PEPTIDE-[Oxidation]")

    def test_nterm_variable_modifications(self):
        """Test N-terminal variable modifications."""
        results = list(
            modify(
                self.annotation,
                nterm_variable={"P": [self.acetyl_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results), 2)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("[Acetyl]-PEPTIDE", serialized)

    def test_cterm_variable_modifications(self):
        """Test C-terminal variable modifications."""
        results = list(
            modify(
                self.annotation,
                cterm_variable={"E": [self.oxidation_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results), 2)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("PEPTIDE-[Oxidation]", serialized)

    def test_labile_static_modifications(self):
        """Test labile static modifications."""
        results = list(
            modify(
                self.annotation,
                labile_static={"P": [self.phospho_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        # Labile mods should appear in serialization
        serialized = results[0].serialize()
        self.assertIn("{Phospho}{Phospho}PEPTIDE", serialized)

    def test_labile_variable_modifications(self):
        """Test labile variable modifications."""
        results = list(
            modify(
                self.annotation,
                labile_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate combinations with labile mods
        self.assertEqual(len(results), 3)
        serialized = [r.serialize() for r in results]

        serialized_counter = {s: serialized.count(s) for s in serialized}
        self.assertEqual(serialized_counter.get("PEPTIDE", 0), 1)
        self.assertEqual(serialized_counter.get("{Phospho}PEPTIDE", 0), 2)

    def test_max_variable_mods_limit(self):
        """Test max_variable_mods parameter limits combinations."""
        # With max_variable_mods=1, should get 3 results for 2 possible sites
        results_1 = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results_1), 3)  # no mods, P@0, P@2

        # With max_variable_mods=2, should get 4 results
        results_2 = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=2,
            )
        )
        serialized = [r.serialize() for r in results_2]
        self.assertEqual(len(results_2), 4)  # no mods, P@0, P@2, P@0+P@2
        self.assertTrue("PEPTIDE" in serialized)
        self.assertTrue("P[Phospho]EPTIDE" in serialized)
        self.assertTrue("PEP[Phospho]TIDE" in serialized)
        self.assertTrue("P[Phospho]EP[Phospho]TIDE" in serialized)

    def test_site_conflict_prevention(self):
        """Test that site conflicts are prevented."""
        # Two different mods for the same amino acid should not conflict
        # since they're applied to the same sites
        results = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod, self.oxidation_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate: no mods, phospho@0, phospho@2, oxidation@0, oxidation@2
        self.assertEqual(len(results), 5)
        serializxed = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serializxed)
        self.assertIn("P[Phospho]EPTIDE", serializxed)
        self.assertIn("PEP[Phospho]TIDE", serializxed)
        self.assertIn("P[Oxidation]EPTIDE", serializxed)
        self.assertIn("PEP[Oxidation]TIDE", serializxed)

    def test_multiple_variable_mod_types(self):
        """Test with multiple types of variable modifications."""
        results = list(
            modify(
                self.annotation,
                nterm_variable={"P": [self.acetyl_mod]},
                internal_variable={"T": [self.phospho_mod]},
                cterm_variable={"E": [self.oxidation_mod]},
                max_variable_mods=2,
            )
        )
        # Should generate many combinations
        self.assertEqual(len(results), 7)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("[Acetyl]-PEPTIDE", serialized)
        self.assertIn("PEPT[Phospho]IDE", serialized)
        self.assertIn("PEPTIDE-[Oxidation]", serialized)
        self.assertIn("[Acetyl]-PEPT[Phospho]IDE", serialized)
        self.assertIn("[Acetyl]-PEPTIDE-[Oxidation]", serialized)
        self.assertIn("PEPT[Phospho]IDE-[Oxidation]", serialized)

    def test_use_regex_parameter(self):
        """Test use_regex parameter is passed through."""
        # This mainly tests that the parameter doesn't cause errors
        results = list(
            modify(
                self.annotation,
                internal_static={"P": [self.phospho_mod]},
                max_variable_mods=0,
                use_regex=False,
            )
        )
        self.assertEqual(len(results), 1)

    def test_empty_sequence_annotation(self):
        """Test with empty sequence annotation."""
        empty_annotation = pt.ProFormaAnnotation(sequence="")
        results = list(modify(empty_annotation, max_variable_mods=0))
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "")

    def test_single_amino_acid_sequence(self):
        """Test with single amino acid sequence."""
        single_aa = pt.ProFormaAnnotation(sequence="P")
        results = list(
            modify(
                single_aa,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results), 2)
        serialized = [r.serialize() for r in results]
        self.assertIn("P", serialized)
        self.assertIn("P[Phospho]", serialized)

    def test_condense_peptidoform(self):
        """Test condense_to_peptidoform function."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        condensed = annotation.condense_to_peptidoform()
        self.assertEqual(condensed.serialize(), "PEPTIDE")

    def test_condense_peptidoform_mod(self):
        """Test condense_to_peptidoform function."""
        annotation = pt.ProFormaAnnotation.parse("P[100]EPTIDE")
        condensed = annotation.condense_to_peptidoform()
        self.assertEqual(condensed.serialize(), "[100]?PEPTIDE")

    def test_condense_peptidoform_multi_mod(self):
        """Test condense_to_peptidoform function."""
        annotation = pt.ProFormaAnnotation.parse("[3.14]-P[100]EPTID[10]E-[Oxidation]")
        condensed = annotation.condense_to_peptidoform()
        # Note: modifications may be in different order, and only one ? follows all unknown mods
        self.assertEqual(condensed.serialize(), "[100][10]?[3.14]-PEPTIDE-[Oxidation]")


class TestBuildModsFirstProteoformOnly(unittest.TestCase):
    def setUp(self):
        self.annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        self.phospho_mod = "Phospho"
        self.oxidation_mod = "Oxidation"
        self.acetyl_mod = "Acetyl"
        self.methyl_mod = "Methyl"

    def test_first_proteoform_only_removes_positional_duplicates(self):
        """Test that first_proteoform_only=True removes positional duplicates."""
        # Without first_proteoform_only: should get both P@0 and P@2
        results_all = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=False,
            )
        )
        self.assertEqual(len(results_all), 3)  # no mods, P@0, P@2

        # With first_proteoform_only: should get only one P modification
        results_unique = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )
        self.assertEqual(len(results_unique), 2)  # no mods, P (position-independent)

        serialized = [r.serialize() for r in results_unique]
        self.assertIn("PEPTIDE", serialized)
        # Should have exactly one phosphorylated P (at either position)
        phospho_count = sum(1 for s in serialized if "Phospho" in s)
        self.assertEqual(phospho_count, 1)

    def test_first_proteoform_only_with_two_modifications(self):
        """Test that first_proteoform_only works with max_variable_mods=2."""
        results = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=2,
                unique_peptidoforms=True,
            )
        )
        # Should get: no mods, 1 P, 2 P (but not different positional combinations)
        self.assertEqual(len(results), 3)

        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)

        # Count how many have phospho
        phospho_counts = {}
        for s in serialized:
            count = s.count("Phospho")
            phospho_counts[count] = phospho_counts.get(count, 0) + 1

        self.assertEqual(phospho_counts.get(0, 0), 1)  # no mods
        self.assertEqual(phospho_counts.get(1, 0), 1)  # one P
        self.assertEqual(phospho_counts.get(2, 0), 1)  # two P

    def test_first_proteoform_only_with_different_mod_types(self):
        """Test that first_proteoform_only distinguishes different modification types."""
        results = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod, self.oxidation_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )
        # Should get: no mods, one phospho, one oxidation
        self.assertEqual(len(results), 3)

        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)

        # Should have one of each mod type
        phospho_count = sum(1 for s in serialized if "Phospho" in s)
        oxidation_count = sum(1 for s in serialized if "Oxidation" in s)
        self.assertEqual(phospho_count, 1)
        self.assertEqual(oxidation_count, 1)

    def test_first_proteoform_only_with_mixed_mods(self):
        """Test first_proteoform_only with multiple mod types on different residues."""
        results = list(
            modify(
                self.annotation,
                internal_variable={
                    "P": [self.phospho_mod],
                    "E": [self.oxidation_mod],
                },
                max_variable_mods=2,
                unique_peptidoforms=True,
            )
        )

        serialized = [r.serialize() for r in results]

        self.assertEqual(len(results), 6)

        self.assertIn("PEPTIDE", serialized)

        # Check that we have the right combinations
        has_phospho_only = any("Phospho" in s and "Oxidation" not in s for s in serialized)
        has_oxidation_only = any("Oxidation" in s and "Phospho" not in s for s in serialized)
        has_both = any("Phospho" in s and "Oxidation" in s for s in serialized)

        self.assertTrue(has_phospho_only)
        self.assertTrue(has_oxidation_only)
        self.assertTrue(has_both)

    def test_first_proteoform_only_with_nterm_mods(self):
        """Test first_proteoform_only with N-terminal modifications."""
        results = list(
            modify(
                self.annotation,
                nterm_variable={"P": [self.acetyl_mod]},
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        # Should get: no mods, nterm acetyl, internal phospho
        self.assertEqual(len(results), 3)

        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)

        # Check we have one of each
        acetyl_count = sum(1 for s in serialized if "Acetyl" in s)
        phospho_count = sum(1 for s in serialized if "Phospho" in s)
        self.assertEqual(acetyl_count, 1)
        self.assertEqual(phospho_count, 1)

    def test_first_proteoform_only_with_cterm_mods(self):
        """Test first_proteoform_only with C-terminal modifications."""
        results = list(
            modify(
                self.annotation,
                cterm_variable={"E": [self.oxidation_mod]},
                internal_variable={"E": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        # Should get: no mods, cterm oxidation, internal phospho
        self.assertEqual(len(results), 3)

        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)

    def test_first_proteoform_only_with_static_mods(self):
        """Test that static mods don't affect first_proteoform_only logic."""
        results_without = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        results_with_static = list(
            modify(
                self.annotation,
                internal_static={"T": [self.methyl_mod]},
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        # Should have same number of results (static doesn't affect variable uniqueness)
        self.assertEqual(len(results_without), len(results_with_static))

    def test_first_proteoform_only_with_labile_mods(self):
        """Test first_proteoform_only with labile modifications."""
        results = list(
            modify(
                self.annotation,
                labile_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        # Should get: no mods, one labile phospho (position-independent)
        self.assertEqual(len(results), 2)

        serialized = [r.serialize() for r in results]
        labile_count = sum(1 for s in serialized if "{Phospho}" in s)
        self.assertEqual(labile_count, 1)

    def test_first_proteoform_only_complex_scenario(self):
        """Test first_proteoform_only with complex modification scenario."""
        # Sequence with multiple mod sites
        complex_annotation = pt.ProFormaAnnotation(sequence="MPEPTIPEP")

        results = list(
            modify(
                complex_annotation,
                internal_variable={
                    "M": [self.oxidation_mod],
                    "P": [self.phospho_mod],
                },
                max_variable_mods=2,
                unique_peptidoforms=True,
            )
        )

        # Should get unique combinations by mod type count, not position:
        # - no mods
        # - 1 oxidation
        # - 1 phospho
        # - 2 phospho
        # - 1 oxidation + 1 phospho
        self.assertEqual(len(results), 5)

    def test_first_proteoform_only_false_gives_all_combinations(self):
        """Test that first_proteoform_only=False gives all positional combinations."""
        results = list(
            modify(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=2,
                unique_peptidoforms=False,
            )
        )

        # Should get: no mods, P@0, P@2, P@0+P@2
        self.assertEqual(len(results), 4)

        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("P[Phospho]EPTIDE", serialized)
        self.assertIn("PEP[Phospho]TIDE", serialized)
        self.assertIn("P[Phospho]EP[Phospho]TIDE", serialized)

    def test_first_proteoform_only_with_numeric_mods(self):
        """Test first_proteoform_only with numeric modification values."""
        results = list(
            modify(
                self.annotation,
                internal_variable={"P": [15.99]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        # Should get: no mods, one +15.99 (position-independent)
        self.assertEqual(len(results), 2)

    def test_first_proteoform_only_with_string_mods(self):
        """Test first_proteoform_only with string modification values."""
        results = list(
            modify(
                self.annotation,
                internal_variable={"P": ["CustomMod"]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        # Should get: no mods, one CustomMod (position-independent)
        self.assertEqual(len(results), 2)

    def test_first_proteoform_only_empty_sequence(self):
        """Test first_proteoform_only with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation(sequence="")
        results = list(
            modify(
                empty_annotation,
                max_variable_mods=0,
                unique_peptidoforms=True,
            )
        )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "")

    def test_first_proteoform_only_single_amino_acid(self):
        """Test first_proteoform_only with single amino acid sequence."""
        single_aa = pt.ProFormaAnnotation(sequence="P")
        results = list(
            modify(
                single_aa,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
                unique_peptidoforms=True,
            )
        )

        # Should get: no mods, one phospho (only one position anyway)
        self.assertEqual(len(results), 2)

    def test_first_proteoform_only_triple_mods_same_type(self):
        """Test first_proteoform_only with three mods of the same type."""
        triple_p_annotation = pt.ProFormaAnnotation(sequence="PPPTIDE")

        results = list(
            modify(
                triple_p_annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=3,
                unique_peptidoforms=True,
            )
        )

        # Should get: 0 mods, 1 phospho, 2 phospho, 3 phospho
        self.assertEqual(len(results), 4)

        serialized = [r.serialize() for r in results]
        phospho_count_distribution = {}
        for s in serialized:
            count = s.count("Phospho")
            phospho_count_distribution[count] = phospho_count_distribution.get(count, 0) + 1

        # Each count should appear exactly once
        self.assertEqual(phospho_count_distribution, {0: 1, 1: 1, 2: 1, 3: 1})

    def test_first_proteoform_only_nterm_cterm_combination(self):
        """Test first_proteoform_only with both nterm and cterm modifications."""
        results = list(
            modify(
                self.annotation,
                nterm_variable={"P": [self.acetyl_mod]},
                cterm_variable={"E": [self.oxidation_mod]},
                max_variable_mods=2,
                unique_peptidoforms=True,
            )
        )

        # Should get: no mods, nterm only, cterm only, both
        self.assertEqual(len(results), 4)

        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)

        has_nterm = any("[Acetyl]-" in s and "Oxidation" not in s for s in serialized)
        has_cterm = any("-[Oxidation]" in s and "Acetyl" not in s for s in serialized)
        has_both = any("[Acetyl]-" in s and "-[Oxidation]" in s for s in serialized)

        self.assertTrue(has_nterm)
        self.assertTrue(has_cterm)
        self.assertTrue(has_both)


if __name__ == "__main__":
    unittest.main()

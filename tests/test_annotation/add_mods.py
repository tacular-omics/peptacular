# filepath: /workspaces/peptacular/tests/test_proforma_funcs/add_mods.py
import unittest

import peptacular as pt


class TestAddMods(unittest.TestCase):
    """
    TESTS FOR: add_labile_mods
    """

    def test_add_labile_mods_replace_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        result = annotation.set_labile_mods("HexNAc", inplace=True)

        self.assertIs(result, annotation)  # Should return same object
        self.assertEqual(annotation.serialize(), "{HexNAc}PEPTIDE")

    def test_add_labile_mods_replace_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        result = annotation.set_labile_mods("HexNAc", inplace=False)

        self.assertIsNot(result, annotation)  # Should return different object
        self.assertEqual(result.serialize(), "{HexNAc}PEPTIDE")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "{Glycan}PEPTIDE")

    def test_add_labile_mods_append_to_existing(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        annotation.append_labile_mod("HexNAc", inplace=True)

        self.assertEqual(annotation.serialize(), "{Glycan}{HexNAc}PEPTIDE")

    def test_add_labile_mods_append_to_empty(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.append_labile_mod("Glycan", inplace=True)

        self.assertEqual(annotation.serialize(), "{Glycan}PEPTIDE")

    def test_add_labile_mods_list_input(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.set_labile_mods(["Glycan", "HexNAc"], inplace=True)

        self.assertEqual(annotation.serialize(), "{Glycan}{HexNAc}PEPTIDE")

    def test_add_labile_mods_mod_object_input(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mod = pt.Mod("Glycan", 1)
        annotation.set_labile_mods(mod, inplace=True)

        self.assertEqual(annotation.serialize(), "{Glycan}PEPTIDE")

    def test_add_labile_mods_none_input_no_append(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        annotation.set_labile_mods(None, inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: add_unknown_mods
    """

    def test_add_unknown_mods_replace_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        result = annotation.set_unknown_mods("Mystery", inplace=True)

        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Mystery]?PEPTIDE")

    def test_add_unknown_mods_replace_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        result = annotation.set_unknown_mods("Mystery", inplace=False)

        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Mystery]?PEPTIDE")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "[Unknown]?PEPTIDE")

    def test_add_unknown_mods_append_to_existing(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        annotation.append_unknown_mod("Mystery", inplace=True)

        self.assertEqual(annotation.serialize(), "[Unknown][Mystery]?PEPTIDE")

    def test_add_unknown_mods_append_to_empty(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.append_unknown_mod("Unknown", inplace=True)

        self.assertEqual(annotation.serialize(), "[Unknown]?PEPTIDE")

    def test_add_unknown_mods_list_input(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.set_unknown_mods(["Unknown1", "Unknown2"], inplace=True)

        self.assertEqual(annotation.serialize(), "[Unknown1][Unknown2]?PEPTIDE")

    def test_add_unknown_mods_none_input_no_append(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        annotation.set_unknown_mods(None, inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: add_nterm_mods
    """

    def test_add_nterm_mods_replace_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        result = annotation.set_nterm_mods("Formyl", inplace=True)

        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Formyl]-PEPTIDE")

    def test_add_nterm_mods_replace_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        result = annotation.set_nterm_mods("Formyl", inplace=False)

        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Formyl]-PEPTIDE")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_nterm_mods_append_to_existing(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        annotation.append_nterm_mod("Formyl", inplace=True)

        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE")

    def test_add_nterm_mods_append_to_empty(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.append_nterm_mod("Acetyl", inplace=True)

        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_nterm_mods_list_input(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.set_nterm_mods(["Acetyl", "Formyl"], inplace=True)

        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE")

    def test_add_nterm_mods_none_input_no_append(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        annotation.set_nterm_mods(None, inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: add_cterm_mods
    """

    def test_add_cterm_mods_replace_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        result = annotation.set_cterm_mods("Hydroxyl", inplace=True)

        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Hydroxyl]")

    def test_add_cterm_mods_replace_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        result = annotation.set_cterm_mods("Hydroxyl", inplace=False)

        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "PEPTIDE-[Hydroxyl]")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide]")

    def test_add_cterm_mods_append_to_existing(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        annotation.append_cterm_mod("Hydroxyl", inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide][Hydroxyl]")

    def test_add_cterm_mods_append_to_empty(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.append_cterm_mod("Amide", inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide]")

    def test_add_cterm_mods_list_input(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.set_cterm_mods(["Amide", "Hydroxyl"], inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide][Hydroxyl]")

    def test_add_cterm_mods_none_input_no_append(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        annotation.set_cterm_mods(None, inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: add_mod_by_type
    """

    def test_add_mod_by_type_string_parameter(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = annotation.set_nterm_mods("Acetyl", inplace=True)

        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_mod_by_type_mod_type_parameter(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = annotation.set_nterm_mods("Acetyl", inplace=True)

        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_mod_by_type_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = annotation.set_nterm_mods("Acetyl", inplace=False)

        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Acetyl]-PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_mod_by_type_all_types(self):
        # Test each mod type
        test_cases = [
            ("Glycan", "labile", "{Glycan}PEPTIDE"),
            ("Unknown", "unknown", "[Unknown]?PEPTIDE"),
            ("Acetyl", "nterm", "[Acetyl]-PEPTIDE"),
            ("Amide", "cterm", "PEPTIDE-[Amide]"),
            ("15N", "isotope", "<15N>PEPTIDE"),
            ("57@C", "static", "<57@C>PEPTIDE"),
        ]

        for mod_val, mod_type, expected_serialization in test_cases:
            with self.subTest(mod_type=mod_type):
                annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
                # dispatch to the corresponding setter to match the 'labile' pattern
                if mod_type == "labile":
                    annotation.set_labile_mods(mod_val, inplace=True)
                elif mod_type == "unknown":
                    annotation.set_unknown_mods(mod_val, inplace=True)
                elif mod_type == "nterm":
                    annotation.set_nterm_mods(mod_val, inplace=True)
                elif mod_type == "cterm":
                    annotation.set_cterm_mods(mod_val, inplace=True)
                elif mod_type == "isotope":
                    annotation.set_isotope_mods(mod_val, inplace=True)
                elif mod_type == "static":
                    annotation.set_static_mods(mod_val, inplace=True)
                self.assertEqual(annotation.serialize(), expected_serialization)

    def test_add_mod_by_type_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotation.set_charge(2, inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE/2")

    def test_add_mod_by_type_charge_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        annotation.append_charge_adduct("Na+", inplace=True)

        self.assertEqual(annotation.serialize(), "PEPTIDE/2[Na+]")

    def test_add_mod_by_type_internal(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        # append an internal mod (new API appends a single internal mod by default at position 0)
        annotation.append_internal_mod_at_index(0, "Phospho", inplace=True)

        # Internal mods are added at position 0 by default in add_mod_by_type
        self.assertEqual(annotation.serialize(), "P[Phospho]EPTIDE")

    def test_add_mod_by_type_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        interval = pt.Interval(1, 3, ambiguous=False, mods=[pt.Mod("Phospho", 1)])
        # set interval-style mods using a specific setter
        annotation.set_intervals(interval, inplace=True)

        self.assertEqual(annotation.serialize(), "P(EP)[Phospho]TIDE")

    def test_add_mod_by_type_append_behavior(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        annotation.append_nterm_mod("Formyl", inplace=True)

        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE")

    """
    TESTS FOR: add_mod_dict
    """

    def test_add_mod_dict_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mod_dict = {"nterm": "Acetyl", "cterm": "Amide", "labile": "Glycan"}
        result = annotation.set_mod_dict(mod_dict, append=False, inplace=True)

        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "{Glycan}[Acetyl]-PEPTIDE-[Amide]")

    def test_add_mod_dict_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mod_dict = {"nterm": "Acetyl", "cterm": "Amide"}
        result = annotation.set_mod_dict(mod_dict, append=False, inplace=False)

        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Acetyl]-PEPTIDE-[Amide]")
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_mod_dict_with_internal_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mod_dict = {
            "nterm": "Acetyl",
            1: "Phospho",  # Internal mod at position 1
            3: "Methyl",  # Internal mod at position 3
        }
        annotation.set_mod_dict(mod_dict, append=False, inplace=True)

        self.assertEqual(annotation.serialize(), "[Acetyl]-PE[Phospho]PT[Methyl]IDE")

    def test_add_mod_dict_append_behavior(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        mod_dict = {"nterm": "Formyl", "cterm": "Amide"}
        annotation.set_mod_dict(mod_dict, append=True, inplace=True)

        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE-[Amide]")

    def test_add_mod_dict_replace_behavior(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        mod_dict = {"nterm": "Formyl", "cterm": "Hydroxyl"}
        annotation.set_mod_dict(mod_dict, append=False, inplace=True)

        self.assertEqual(annotation.serialize(), "[Formyl]-PEPTIDE-[Hydroxyl]")

    """
    TESTS FOR: edge cases and validation
    """

    def test_add_methods_preserve_sequence(self):
        original_sequence = "PEPTIDE"

        # Test each add method preserves the sequence
        test_cases = [
            ("set_labile_mods", "Glycan"),
            ("set_unknown_mods", "Unknown"),
            ("set_nterm_mods", "Acetyl"),
            ("set_cterm_mods", "Amide"),
        ]

        for method_name, mod_val in test_cases:
            with self.subTest(method=method_name):
                annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
                method = getattr(annotation, method_name)
                method(mod_val, append=False, inplace=True)
                self.assertEqual(annotation.sequence, original_sequence)

    def test_add_methods_return_self_when_inplace(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")

        methods = [
            ("set_labile_mods", "Glycan"),
            ("set_unknown_mods", "Unknown"),
            ("set_nterm_mods", "Acetyl"),
            ("set_cterm_mods", "Amide"),
        ]

        for method_name, mod_val in methods:
            with self.subTest(method=method_name):
                method = getattr(annotation, method_name)
                result = method(mod_val, append=False, inplace=True)
                self.assertIs(result, annotation)

    def test_add_methods_return_copy_when_not_inplace(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")

        methods = [
            ("set_labile_mods", "Glycan"),
            ("set_unknown_mods", "Unknown"),
            ("set_nterm_mods", "Acetyl"),
            ("set_cterm_mods", "Amide"),
        ]

        for method_name, mod_val in methods:
            with self.subTest(method=method_name):
                method = getattr(annotation, method_name)
                result = method(mod_val, append=False, inplace=False)
                self.assertIsNot(result, annotation)

    def test_add_methods_with_multipliers(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mod_with_multiplier = pt.Mod("Phospho", 2)

        annotation.set_nterm_mods(mod_with_multiplier, inplace=True)

        self.assertEqual(annotation.serialize(), "[Phospho]^2-PEPTIDE")

    def test_add_empty_list_to_existing_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        original_serialization = annotation.serialize()

        # use extend for appending lists (append single mods, extend lists)
        annotation.extend_nterm_mods([], inplace=True)

        # Adding empty list should not change the serialization
        self.assertEqual(annotation.serialize(), original_serialization)

    def test_complex_mod_dict_scenario(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        complex_mod_dict = {
            "labile": ["Glycan", "HexNAc"],
            "nterm": "Acetyl",
            "cterm": ["Amide", "Hydroxyl"],
            0: "Phospho",  # Internal mod at position 0
            2: ["Methyl", "Oxidation"],  # Multiple internal mods at position 2
            "charge": 2,
            "isotope": "15N",
        }

        annotation.set_mod_dict(complex_mod_dict, append=False, inplace=True)

        # Verify the complex modification produces expected serialization
        self.assertEqual(
            annotation.serialize(),
            "{Glycan}{HexNAc}<15N>[Acetyl]-P[Phospho]EP[Methyl][Oxidation]TIDE-[Amide][Hydroxyl]/2",
        )


if __name__ == "__main__":
    unittest.main()

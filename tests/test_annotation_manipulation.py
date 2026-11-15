"""Coverage for peptacular.annotation.manipulation: coverage arrays, residue
counting, peptidoform condensing, and subsequence search helpers.
"""

import pytest

import peptacular as pt


class TestCoverage:
    def test_ignores_ambiguous_regions_by_default(self):
        a = pt.parse("(?PEP)TIDE")
        subs = [pt.parse("PEP"), pt.parse("TIDE")]
        assert a.coverage(subs) == [0, 0, 0, 1, 1, 1, 1]

    def test_ignore_ambiguity_true_still_counts_matches(self):
        a = pt.parse("(?PEP)TIDE")
        subs = [pt.parse("PEP"), pt.parse("TIDE")]
        result = a.coverage(subs, ignore_ambiguity=True)
        assert result[0:3] == [1, 1, 1]


class TestCondenseToPeptidoform:
    def test_interval_with_mods_moves_to_unknown(self):
        a = pt.parse("(PEP)[Phospho]TIDE")
        result = a.condense_to_peptidoform(inplace=False)
        assert result.serialize() == "[Phospho]?PEPTIDE"

    def test_interval_without_mods_is_dropped(self):
        a = pt.parse("(PEP)TIDE")
        result = a.condense_to_peptidoform(inplace=False)
        assert result.serialize() == "PEPTIDE"


class TestCountResidues:
    def test_charge_is_not_part_of_a_key(self):
        assert pt.parse("PEPTIDE/2").count_residues() == pt.parse("PEPTIDE").count_residues()

    def test_charge_keeps_the_modified_residue_key(self):
        counts = pt.parse("PEM[Oxidation]TIDE/2").count_residues()
        assert counts["M[Oxidation]"] == 1
        assert not any("/" in residue for residue in counts)

    def test_charge_adducts_are_not_part_of_a_key(self):
        assert pt.parse("PEPTIDE/[M+2H]").count_residues() == pt.parse("PEPTIDE").count_residues()

    def test_labile_modification_is_not_part_of_a_key(self):
        counts = pt.parse("{+79.966}PEPTIDE").count_residues()
        assert counts == {"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}

    def test_unlocalised_modification_is_not_part_of_a_key(self):
        counts = pt.parse("[Phospho]?PEPTIDE").count_residues()
        assert counts == {"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}

    def test_names_are_not_part_of_a_key(self):
        plain = {"P": 2, "E": 2, "T": 1, "I": 1, "D": 1}
        assert pt.parse("(>Peptide)PEPTIDE/2").count_residues() == plain
        assert pt.parse("(>>Ion)PEPTIDE/2").count_residues() == plain
        assert pt.parse("(>>>Compound)(>Peptide)PEPTIDE/2").count_residues() == plain

    def test_terminal_modifications_stay_with_their_residue(self):
        counts = pt.parse("[Acetyl]-PEPTIDE-[Amidated]").count_residues()
        assert counts == {"[Acetyl]-P": 1, "E": 1, "P": 1, "T": 1, "I": 1, "D": 1, "E-[Amidated]": 1}

    def test_include_mods_false_ignores_every_annotation(self):
        assert pt.parse("PEM[Oxidation]TIDE/2").count_residues(include_mods=False) == {
            "P": 1,
            "E": 2,
            "M": 1,
            "T": 1,
            "I": 1,
            "D": 1,
        }

    def test_percent_residues_follows_the_same_keys(self):
        assert pt.parse("PEPTIDE/2").percent_residues() == pt.parse("PEPTIDE").percent_residues()

    def test_global_isotope_label_stays_in_the_key(self):
        # A global isotope label applies to every residue, so unlike the charge and the
        # labile, unlocalised and name annotations it belongs in the key. The docstring
        # of count_residues promises this; without a test the promise cannot fail.
        assert pt.parse("<13C>PEPTIDE").count_residues() == {"<13C>P": 2, "<13C>E": 2, "<13C>T": 1, "<13C>I": 1, "<13C>D": 1}
        assert pt.parse("<13C>PEPTIDE/2").count_residues() == pt.parse("<13C>PEPTIDE").count_residues()


class TestPercentResidues:
    def test_empty_sequence_returns_empty_dict(self):
        a = pt.parse("")
        assert a.percent_residues() == {}


class TestIsSubsequence:
    def test_sequence_not_present_returns_false(self):
        a = pt.parse("XYZ")
        other = pt.parse("ABCDEF")
        assert a.is_subsequence(other) is False

    def test_exact_full_annotation_match_without_ignoring_intervals(self):
        sub = pt.parse("PEP")
        whole = pt.parse("PEPTIDE")
        assert sub.is_subsequence(whole, ignore_intervals=False) is True

    def test_ignore_mods_true_matches_on_sequence_only(self):
        sub = pt.parse("PEP")
        whole = pt.parse("P[Oxidation]EPTIDE")
        assert sub.is_subsequence(whole, ignore_mods=True) is True


class TestFindIndices:
    def test_non_annotation_other_raises_type_error(self):
        sub = pt.parse("PEP")
        with pytest.raises(TypeError, match="must be a ProFormaAnnotation"):
            sub.find_indices("PEPTIDE")  # type: ignore[arg-type]

    def test_finds_overlapping_and_non_overlapping_matches(self):
        sub = pt.parse("PEP")
        whole = pt.parse("PEPPEPTIDE")
        assert sub.find_indices(whole) == [0, 3]


_PARITY_CASES = [
    ("plain", "PEPTIDE", ["PEP", "TID"]),
    ("mod_prot_unmod_pep", "P[+80]EPTIDE", ["PEP", "TID"]),
    ("mod_pep_mod_prot", "P[+80]EPTIDE", ["P[+80]EP", "PEP", "TID"]),
    ("multiple_occurrence", "PEPTIDEPEP", ["PEP"]),
    ("overlapping", "IIIII", ["II"]),
    ("ambig_sub", "PEPTIDE", ["P(?EP)"]),
    ("ambig_prot", "P(?EP)TIDE", ["PEP"]),
]

# The coverage vector each case returned before the string-search fast path, measured on
# f9907bc (the commit before #101) for every case and option combination below. These are
# literals on purpose: comparing the annotation method against the sequence wrapper only
# shows that the two entry points agree, and both of them go through the fast path, so such
# a comparison cannot see the fast path drift away from the answer it replaced.
_PRE_FAST_PATH_COVERAGE: dict[tuple[str, bool, bool, bool], list[int]] = {
    ("plain", False, False, False): [1, 1, 1, 1, 1, 1, 0],
    ("plain", False, False, True): [1, 1, 1, 1, 1, 1, 0],
    ("plain", False, True, False): [1, 1, 1, 1, 1, 1, 0],
    ("plain", False, True, True): [1, 1, 1, 1, 1, 1, 0],
    ("plain", True, False, False): [1, 1, 1, 1, 1, 1, 0],
    ("plain", True, False, True): [1, 1, 1, 1, 1, 1, 0],
    ("plain", True, True, False): [1, 1, 1, 1, 1, 1, 0],
    ("plain", True, True, True): [1, 1, 1, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", False, False, False): [0, 0, 0, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", False, False, True): [0, 0, 0, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", False, True, False): [1, 1, 1, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", False, True, True): [1, 1, 1, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", True, False, False): [0, 0, 0, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", True, False, True): [0, 0, 0, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", True, True, False): [1, 1, 1, 1, 1, 1, 0],
    ("mod_prot_unmod_pep", True, True, True): [1, 1, 1, 1, 1, 1, 0],
    ("mod_pep_mod_prot", False, False, False): [1, 1, 1, 1, 1, 1, 0],
    ("mod_pep_mod_prot", False, False, True): [1, 1, 1, 1, 1, 1, 0],
    ("mod_pep_mod_prot", False, True, False): [1, 1, 1, 1, 1, 1, 0],
    ("mod_pep_mod_prot", False, True, True): [1, 1, 1, 1, 1, 1, 0],
    ("mod_pep_mod_prot", True, False, False): [1, 1, 1, 1, 1, 1, 0],
    ("mod_pep_mod_prot", True, False, True): [1, 1, 1, 1, 1, 1, 0],
    ("mod_pep_mod_prot", True, True, False): [2, 2, 2, 1, 1, 1, 0],
    ("mod_pep_mod_prot", True, True, True): [2, 2, 2, 1, 1, 1, 0],
    ("multiple_occurrence", False, False, False): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("multiple_occurrence", False, False, True): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("multiple_occurrence", False, True, False): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("multiple_occurrence", False, True, True): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("multiple_occurrence", True, False, False): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("multiple_occurrence", True, False, True): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("multiple_occurrence", True, True, False): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("multiple_occurrence", True, True, True): [1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    ("overlapping", False, False, False): [1, 1, 1, 1, 1],
    ("overlapping", False, False, True): [1, 1, 1, 1, 1],
    ("overlapping", False, True, False): [1, 1, 1, 1, 1],
    ("overlapping", False, True, True): [1, 1, 1, 1, 1],
    ("overlapping", True, False, False): [1, 2, 2, 2, 1],
    ("overlapping", True, False, True): [1, 2, 2, 2, 1],
    ("overlapping", True, True, False): [1, 2, 2, 2, 1],
    ("overlapping", True, True, True): [1, 2, 2, 2, 1],
    ("ambig_sub", False, False, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_sub", False, False, True): [1, 1, 1, 0, 0, 0, 0],
    ("ambig_sub", False, True, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_sub", False, True, True): [1, 1, 1, 0, 0, 0, 0],
    ("ambig_sub", True, False, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_sub", True, False, True): [1, 1, 1, 0, 0, 0, 0],
    ("ambig_sub", True, True, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_sub", True, True, True): [1, 1, 1, 0, 0, 0, 0],
    ("ambig_prot", False, False, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_prot", False, False, True): [1, 1, 1, 0, 0, 0, 0],
    ("ambig_prot", False, True, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_prot", False, True, True): [1, 1, 1, 0, 0, 0, 0],
    ("ambig_prot", True, False, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_prot", True, False, True): [1, 1, 1, 0, 0, 0, 0],
    ("ambig_prot", True, True, False): [1, 0, 0, 0, 0, 0, 0],
    ("ambig_prot", True, True, True): [1, 1, 1, 0, 0, 0, 0],
}


class TestCoverageFastPathParity:
    """The string search fast path returns what annotation slicing returned, option for option."""

    def test_expected_table_covers_every_case_and_option(self):
        assert len(_PRE_FAST_PATH_COVERAGE) == len(_PARITY_CASES) * 8
        assert {key[0] for key in _PRE_FAST_PATH_COVERAGE} == {case[0] for case in _PARITY_CASES}

    @pytest.mark.parametrize(("name", "prot_str", "peps_str"), _PARITY_CASES)
    @pytest.mark.parametrize("accumulate", [False, True])
    @pytest.mark.parametrize("ignore_mods", [False, True])
    @pytest.mark.parametrize("ignore_ambiguity", [False, True])
    def test_parity_across_options(
        self,
        name: str,
        prot_str: str,
        peps_str: list[str],
        accumulate: bool,
        ignore_mods: bool,
        ignore_ambiguity: bool,
    ):
        expected = _PRE_FAST_PATH_COVERAGE[name, accumulate, ignore_mods, ignore_ambiguity]

        prot = pt.parse(prot_str)
        peps = [pt.parse(p) for p in peps_str]

        # The annotation method, against the value measured before the fast path existed.
        cov = prot.coverage(
            peps,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )
        assert len(cov) == len(prot)
        assert cov == expected

        # The public sequence wrapper reaches the same answer.
        cov_public = pt.coverage(
            prot_str,
            peps_str,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )
        assert cov_public == expected

    def test_overlapping_repeats_find_all_positions(self):
        # Zero-width lookahead ensures overlapping occurrences are all identified
        assert pt.parse("II").find_indices(pt.parse("IIIII")) == [0, 1, 2, 3]
        assert pt.coverage("IIIII", ["II"]) == [1, 1, 1, 1, 1]
        assert pt.coverage("IIIII", ["II"], accumulate=True) == [1, 2, 2, 2, 1]

    def test_ambiguous_intervals_zeroed_out(self):
        # Ambiguous interval in subpeptide stays at 0
        assert pt.coverage("PEPTIDE", ["P(?EP)"]) == [1, 0, 0, 0, 0, 0, 0]
        # Ambiguous interval in protein stays at 0
        assert pt.coverage("P(?EP)TIDE", ["PEP"]) == [1, 0, 0, 0, 0, 0, 0]


class TestCoverageIgnoreModsSkipsSlicing:
    """`ignore_mods=True` no longer slices the parent, so a cut it refused is no longer reached.

    Before the fast path, `coverage` sliced the parent annotation once per match even when the
    caller had asked for modifications to be ignored, and the slice raised whenever it would cut
    an interval or a localization group. Measured on f9907bc, each call below raised
    `PeptacularError`; the fast path reaches the answer without a slice. With `ignore_mods=False`
    the slow path still runs and still refuses, which is what keeps the refusal meaningful.
    """

    @pytest.mark.parametrize(
        ("prot_str", "peps_str", "expected"),
        [
            ("PE(?PT)IDE", ["PEP"], [1, 1, 0, 0, 0, 0, 0]),
            ("PE(PT)IDE", ["PEP"], [1, 1, 1, 0, 0, 0, 0]),
            ("PEPT[Phospho#g1]IDEK[#g1]", ["PEPT"], [1, 1, 1, 1, 0, 0, 0, 0]),
        ],
    )
    def test_ignore_mods_returns_a_vector_where_slicing_refused(self, prot_str: str, peps_str: list[str], expected: list[int]):
        assert pt.coverage(prot_str, peps_str, ignore_mods=True) == expected

    @pytest.mark.parametrize(
        ("prot_str", "peps_str", "message"),
        [
            ("PE(?PT)IDE", ["PEP"], "Interval"),
            ("PE(PT)IDE", ["PEP"], "Interval"),
            ("PEPT[Phospho#g1]IDEK[#g1]", ["PEPT"], "Localization group"),
        ],
    )
    def test_ignore_mods_false_still_refuses_the_cut(self, prot_str: str, peps_str: list[str], message: str):
        with pytest.raises(pt.PeptacularError, match=message):
            pt.coverage(prot_str, peps_str, ignore_mods=False)

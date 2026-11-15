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


class TestCoverageFastPathParity:
    """Verify behavioral parity of the string search fast path on coverage across all parameter combinations."""

    @pytest.mark.parametrize(
        ("name", "prot_str", "peps_str"),
        [
            ("plain", "PEPTIDE", ["PEP", "TID"]),
            ("mod_prot_unmod_pep", "P[+80]EPTIDE", ["PEP", "TID"]),
            ("mod_pep_mod_prot", "P[+80]EPTIDE", ["P[+80]EP", "PEP", "TID"]),
            ("multiple_occurrence", "PEPTIDEPEP", ["PEP"]),
            ("overlapping", "IIIII", ["II"]),
            ("ambig_sub", "PEPTIDE", ["P(?EP)"]),
            ("ambig_prot", "P(?EP)TIDE", ["PEP"]),
        ],
    )
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
        prot = pt.parse(prot_str)
        peps = [pt.parse(p) for p in peps_str]

        # Calculate via annotation method
        cov = prot.coverage(
            peps,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )
        assert len(cov) == len(prot)

        # Calculate via public sequence API
        cov_public = pt.coverage(
            prot_str,
            peps_str,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )
        assert cov == cov_public

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

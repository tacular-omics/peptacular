"""Tests for peptacular.annotation.ambiguity.

Covers the implementation that backs the public ``annotate_ambiguity`` /
``condense_ambiguity_to_xnotation`` wrappers exposed on ``ProFormaAnnotation``
and in ``peptacular.sequence.basic``, plus module-level helpers that have no
public wrapper (``group_by_ambiguity``, ``unique_fragments``) and a couple of
private helpers that are only reachable directly (see notes on the specific
tests below).
"""

import pytest

import peptacular as pt
from peptacular.annotation.ambiguity import (
    _apply_mass_shift,
    _validate_coverage_lengths,
    annotate_ambiguity,
    condense_ambiguity_to_xnotation,
    group_by_ambiguity,
    unique_fragments,
)

SEQ = "PEPTIDE"


class TestAnnotateAmbiguityNotInplace:
    def test_default_call_returns_new_annotation_and_leaves_original_untouched(self) -> None:
        original = pt.parse(SEQ)
        result = annotate_ambiguity(
            original,
            [0, 1, 1, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0],
            None,
            add_mods_to_intervals=False,
            sort_mods=True,
            inplace=False,
        )
        assert result.serialize() == "(?PE)PTI(?DE)"
        assert original.serialize() == SEQ
        assert result is not original


class TestAnnotateAmbiguityValidation:
    def test_raises_when_annotation_already_has_intervals(self) -> None:
        already_annotated = annotate_ambiguity(
            pt.parse(SEQ),
            [0, 1, 1, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0],
            None,
            add_mods_to_intervals=False,
            sort_mods=True,
            inplace=False,
        )
        assert already_annotated.has_intervals

        with pytest.raises(ValueError, match="should not contain intervals"):
            annotate_ambiguity(
                already_annotated,
                [0] * 7,
                [0] * 7,
                None,
                add_mods_to_intervals=False,
                sort_mods=True,
                inplace=True,
            )

    def test_raises_when_forward_and_reverse_coverage_disagree_with_sequence_length(self) -> None:
        with pytest.raises(ValueError, match="Coverage length does not match sequence length"):
            annotate_ambiguity(
                pt.parse(SEQ),
                [1, 2, 3],
                [1, 2],
                None,
                add_mods_to_intervals=False,
                sort_mods=True,
                inplace=True,
            )


class TestValidateCoverageLengths:
    """Exercises the private length-check helper directly.

    Note: ``len(a) != len(b) != seq_len`` is a chained comparison, so it only
    raises when forward/reverse *also* disagree with each other. A case where
    forward and reverse coverage have equal length to each other but both
    differ from the sequence length (e.g. forward=[1,1,1], reverse=[1,1,1],
    seq_len=7) silently passes validation instead of raising -- this looks
    like a latent bug in ``_validate_coverage_lengths``, flagged separately
    rather than encoded here as expected behavior.
    """

    def test_raises_on_mismatched_lengths(self) -> None:
        with pytest.raises(ValueError, match=r"3 != 2 != 3"):
            _validate_coverage_lengths([1, 2, 3], [1, 2], 3)

    def test_does_not_raise_when_all_lengths_match(self) -> None:
        assert _validate_coverage_lengths([1, 2, 3], [1, 2, 3], 3) is None


class TestAnnotateAmbiguityAddModsToIntervals:
    def test_folds_internal_mod_into_enclosing_ambiguous_interval(self) -> None:
        # No forward coverage at all, and only position 1 unreached by reverse
        # coverage: mass shift lands (as a singleton) on position 0, which
        # itself falls inside the (?PE) ambiguous interval.
        without_condense = annotate_ambiguity(
            pt.parse(SEQ),
            [0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 1, 1, 1],
            120,
            add_mods_to_intervals=False,
            sort_mods=True,
            inplace=False,
        )
        assert without_condense.serialize() == "(?P[+120]E)PTIDE"

        with_condense = annotate_ambiguity(
            pt.parse(SEQ),
            [0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 1, 1, 1],
            120,
            add_mods_to_intervals=True,
            sort_mods=True,
            inplace=False,
        )
        assert with_condense.serialize() == "(?PE)[+120]PTIDE"


class TestAnnotateAmbiguityFullCoverage:
    def test_full_coverage_on_both_strands_produces_no_ambiguity_intervals(self) -> None:
        # Every position is covered from both directions, so there is no
        # overlap of "uncovered" zero-runs at all (the combine step's
        # common-indices set is empty from the start).
        result = annotate_ambiguity(
            pt.parse(SEQ),
            [1] * len(SEQ),
            [1] * len(SEQ),
            None,
            add_mods_to_intervals=False,
            sort_mods=True,
            inplace=False,
        )
        assert result.serialize() == SEQ
        assert not result.has_intervals


class TestAnnotateAmbiguityMassShiftPlacement:
    def test_unlocalizable_mass_shift_becomes_labile_mod(self) -> None:
        result = annotate_ambiguity(
            pt.parse(SEQ),
            [0, 1, 1, 1, 1, 0, 0],
            [0, 0, 1, 1, 1, 1, 0],
            120,
            add_mods_to_intervals=False,
            sort_mods=True,
            inplace=False,
        )
        assert result.serialize() == "{+120}(?PE)PTI(?DE)"

    def test_mass_shift_gap_that_already_matches_an_ambiguous_interval_reuses_it(self) -> None:
        result = annotate_ambiguity(
            pt.parse(SEQ),
            [0, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0],
            120,
            add_mods_to_intervals=False,
            sort_mods=True,
            inplace=False,
        )
        assert result.serialize() == "(?PE)P(?TI)[+120](?DE)"
        # Exactly the ambiguous interval got the mod; no separate interval was created.
        assert len(result.intervals) == 3


class TestApplyMassShiftDirect:
    """Exercises ``_apply_mass_shift`` directly.

    Given how ``annotate_ambiguity`` always builds its ambiguity intervals
    from the very same forward/reverse coverage vectors before calling this
    helper, the "no existing interval matches -> create a new one" branch is
    unreachable through the public entry point (the matching ambiguous
    interval has always already been inserted). Calling the helper directly
    on an annotation with no pre-existing intervals is the only way to
    exercise that branch.
    """

    def test_creates_new_interval_when_no_matching_interval_exists(self) -> None:
        # Regression test: the newly created interval used to hardcode `mods=None`
        # instead of the "interval already exists" branch's `found_int.append_mod
        # (mass_shift)`, silently dropping the shift. It's now applied via
        # `append_mod` after construction, same as the existing-interval branch.
        annotation = pt.parse(SEQ)
        assert not annotation.has_intervals

        _apply_mass_shift(annotation, [1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1], 79.966)

        assert annotation.has_intervals
        assert len(annotation.intervals) == 1
        new_interval = annotation.intervals[0]
        assert (new_interval.start, new_interval.end) == (3, 6)
        assert new_interval.mods.serialize() == "[+79.966]"

    def test_adds_labile_mod_when_no_localization_interval_exists(self) -> None:
        annotation = pt.parse(SEQ)
        # Forward coverage extends past where reverse coverage begins, so the
        # highest-forward/lowest-reverse boundaries cross and no interval can
        # be localized -- this is the ``_get_mass_shift_interval(...) is None`` case.
        _apply_mass_shift(annotation, [0, 1, 1, 1, 1, 0, 0], [0, 0, 1, 1, 1, 1, 0], 42.0)
        assert annotation.has_labile_mods

    def test_adds_internal_mod_at_index_for_singleton_gap(self) -> None:
        annotation = pt.parse(SEQ)
        _apply_mass_shift(annotation, [1, 1, 1, 1, 1, 1, 0], [0, 0, 0, 0, 0, 0, 0], 15.994)
        assert annotation.serialize() == "PEPTIDE[+15.994]"


class TestCondenseAmbiguityToXnotationNotInplace:
    def test_default_call_returns_new_annotation_and_leaves_original_untouched(self) -> None:
        ambiguous = annotate_ambiguity(
            pt.parse(SEQ),
            [0, 1, 1, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0],
            None,
            add_mods_to_intervals=False,
            sort_mods=True,
            inplace=False,
        )
        condensed = condense_ambiguity_to_xnotation(ambiguous, inplace=False)

        assert condensed is not ambiguous
        assert ambiguous.serialize() == "(?PE)PTI(?DE)"
        assert condensed.sequence == "XPTIX"
        assert "?" not in condensed.serialize()


class TestGroupByAmbiguity:
    def test_precision_below_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="Precision must be an integer between 0 and 10"):
            group_by_ambiguity([pt.parse(SEQ)], precision=-1)

    def test_precision_above_ten_raises(self) -> None:
        with pytest.raises(ValueError, match="Precision must be an integer between 0 and 10"):
            group_by_ambiguity([pt.parse(SEQ)], precision=11)

    def test_empty_input_returns_no_groups(self) -> None:
        assert group_by_ambiguity([]) == []

    def test_single_annotation_forms_its_own_group(self) -> None:
        annotation = pt.parse(SEQ)
        groups = group_by_ambiguity([annotation])
        assert groups == [(annotation,)]

    def test_identical_annotations_are_grouped_together(self) -> None:
        a1 = pt.parse(SEQ)
        a2 = pt.parse(SEQ)
        groups = group_by_ambiguity([a1, a2])
        assert len(groups) == 1
        assert groups[0] == (a1, a2)

    def test_annotations_with_disjoint_fragment_masses_form_separate_groups(self) -> None:
        a1 = pt.parse(SEQ)
        a2 = pt.parse(SEQ)
        modified = pt.parse("PEPT[+80]IDE")
        groups = group_by_ambiguity([a1, a2, modified])

        assert len(groups) == 2
        serialized_groups = [{a.serialize() for a in group} for group in groups]
        assert {"PEPTIDE"} in serialized_groups
        assert {"PEPT[+80]IDE"} in serialized_groups


class TestUniqueFragments:
    def test_precision_below_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="Precision must be an integer between 0 and 10"):
            unique_fragments([pt.parse(SEQ)], precision=-1)

    def test_precision_above_ten_raises(self) -> None:
        with pytest.raises(ValueError, match="Precision must be an integer between 0 and 10"):
            unique_fragments([pt.parse(SEQ)], precision=11)

    def test_empty_input_returns_empty_list(self) -> None:
        assert unique_fragments([]) == []

    def test_single_annotation_counts_all_its_masses_as_unique(self) -> None:
        assert unique_fragments([pt.parse(SEQ)]) == [len(SEQ) * 2 - 1]

    def test_identical_annotations_have_no_unique_masses(self) -> None:
        a1 = pt.parse(SEQ)
        a2 = pt.parse(SEQ)
        modified = pt.parse("PEPT[+80]IDE")

        counts = unique_fragments([a1, a2, modified])

        assert counts[0] == 0
        assert counts[1] == 0
        assert counts[2] > 0

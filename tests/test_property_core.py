"""Tests for peptacular.property.core (property/window calculation helpers)."""

import pytest

from peptacular.property.core import (
    aa_property_percentage,
    calc_property,
    calc_window_property,
    charge_at_ph,
    generate_partitions,
    secondary_structure,
)
from peptacular.property.types import AggregationMethod, MissingAAHandling


class TestCalcPropertyMissingAAHandling:
    """A scale missing 'G' forces calc_property through the missing-AA handling paths."""

    scale = {"A": 1.0, "C": 3.0}

    def test_zero_defaults_missing_aa_to_zero(self) -> None:
        assert calc_property("AG", scale=self.scale, missing_aa_handling="zero") == pytest.approx(0.5)

    def test_avg_defaults_missing_aa_to_mean_of_scale(self) -> None:
        assert calc_property("AG", scale=self.scale, missing_aa_handling="avg") == pytest.approx(1.5)

    def test_min_defaults_missing_aa_to_min_of_scale(self) -> None:
        assert calc_property("AG", scale=self.scale, missing_aa_handling="min") == pytest.approx(1.0)

    def test_max_defaults_missing_aa_to_max_of_scale(self) -> None:
        assert calc_property("AG", scale=self.scale, missing_aa_handling="max") == pytest.approx(2.0)

    def test_median_defaults_missing_aa_to_median_of_scale(self) -> None:
        assert calc_property("AG", scale=self.scale, missing_aa_handling="median") == pytest.approx(1.5)

    def test_skip_ignores_missing_aa_entirely(self) -> None:
        # Only 'A' contributes; 'G' is dropped rather than defaulted.
        assert calc_property("AG", scale=self.scale, missing_aa_handling="skip") == pytest.approx(1.0)

    def test_error_raises_on_missing_aa(self) -> None:
        with pytest.raises(ValueError, match="Invalid amino acid: G"):
            calc_property("AG", scale=self.scale, missing_aa_handling="error")

    def test_error_is_the_default_handling(self) -> None:
        with pytest.raises(ValueError, match="Invalid amino acid: G"):
            calc_property("AG", scale=self.scale)

    def test_enum_member_and_string_agree(self) -> None:
        assert calc_property("AG", scale=self.scale, missing_aa_handling=MissingAAHandling.ZERO) == calc_property(
            "AG", scale=self.scale, missing_aa_handling="zero"
        )

    def test_unknown_handling_string_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown MissingAAHandling"):
            calc_property("A", scale=self.scale, missing_aa_handling="bogus")


class TestCalcPropertyAmbiguousAndWildcardAminoAcids:
    """B/J/Z stand in for two constituent amino acids; X stands in for any."""

    scale = {"A": 1.0, "C": 3.0}

    @pytest.mark.parametrize("aa", ["B", "J", "Z"])
    def test_ambiguous_aa_uses_average_of_constituents(self, aa: str) -> None:
        # Neither constituent (D/N, L/I, E/Q) is in the scale, so each constituent
        # defaults to the scale's mean (2.0) under "avg" handling, and B/J/Z average those.
        assert calc_property(aa, scale=self.scale, missing_aa_handling="avg") == pytest.approx(2.0)

    def test_ambiguous_aa_raises_when_no_constituent_has_a_value(self) -> None:
        # "skip" handling makes both constituents resolve to None, leaving nothing to average.
        with pytest.raises(ValueError, match="No valid values found for ambiguous amino acid B"):
            calc_property("B", scale=self.scale, missing_aa_handling="skip")

    def test_x_uses_average_of_entire_scale(self) -> None:
        assert calc_property("X", scale=self.scale, missing_aa_handling="avg") == pytest.approx(2.0)


class TestCalcPropertyNormalization:
    scale = {"A": 1.0, "C": 3.0}

    def test_normalize_scales_value_into_zero_one_range(self) -> None:
        assert calc_property("AC", scale=self.scale, normalize=True) == pytest.approx(0.5)

    def test_normalize_returns_zero_when_scale_has_no_spread(self) -> None:
        flat_scale = {"A": 5.0, "C": 5.0}
        assert calc_property("AC", scale=flat_scale, normalize=True) == pytest.approx(0.0)


class TestCalcPropertyScaleAndAggregation:
    def test_named_scale_string_looks_up_property_scales(self) -> None:
        from_name = calc_property("ACDE", scale="hphob_kyte_doolittle", missing_aa_handling="error")
        from_dict = calc_property(
            "ACDE",
            scale={"A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5},
            missing_aa_handling="error",
        )
        assert from_name == pytest.approx(from_dict)

    def test_unknown_named_scale_raises(self) -> None:
        with pytest.raises(ValueError, match="Scale 'NOT_A_SCALE' not found"):
            calc_property("A", scale="NOT_A_SCALE")

    def test_sum_aggregation(self) -> None:
        scale = {"A": 1.0, "C": 3.0}
        assert calc_property("AC", scale=scale, aggregation_method="sum") == pytest.approx(4.0)

    def test_avg_aggregation_is_default(self) -> None:
        scale = {"A": 1.0, "C": 3.0}
        assert calc_property("AC", scale=scale, aggregation_method=AggregationMethod.AVG) == pytest.approx(2.0)

    def test_invalid_aggregation_method_raises(self) -> None:
        scale = {"A": 1.0, "C": 3.0}
        with pytest.raises(ValueError, match="Invalid aggregation method"):
            calc_property("AC", scale=scale, aggregation_method="bogus")


class TestCalcWindowProperty:
    scale = {"A": 1.0, "C": 2.0, "D": 3.0}

    def test_generates_one_value_per_sliding_window(self) -> None:
        windows = calc_window_property("ACDACD", scale=self.scale, window_size=3)
        assert windows == pytest.approx([2.0, 2.0, 2.0, 2.0])

    def test_empty_sequence_raises(self) -> None:
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            calc_window_property("", scale=self.scale, window_size=3)

    def test_window_size_larger_than_sequence_raises(self) -> None:
        with pytest.raises(ValueError, match="cannot be greater than sequence length"):
            calc_window_property("AC", scale=self.scale, window_size=5)

    def test_non_positive_window_size_raises(self) -> None:
        with pytest.raises(ValueError, match="Window size must be positive"):
            calc_window_property("AC", scale=self.scale, window_size=0)


class TestAaPropertyPercentage:
    def test_computes_fraction_of_specified_residues(self) -> None:
        assert aa_property_percentage("ACDEF", residues=["A", "C"]) == pytest.approx(0.4)

    def test_empty_sequence_returns_zero(self) -> None:
        assert aa_property_percentage("", residues=["A"]) == pytest.approx(0.0)

    def test_residues_absent_from_sequence_contribute_nothing(self) -> None:
        assert aa_property_percentage("AAAA", residues=["W"]) == pytest.approx(0.0)


class TestChargeAtPh:
    def test_empty_sequence_has_zero_charge(self) -> None:
        assert charge_at_ph("", pH=7.0) == pytest.approx(0.0)

    def test_basic_residue_is_positive_at_low_ph(self) -> None:
        assert charge_at_ph("K", pH=2.0) > 0.0

    def test_acidic_residue_is_negative_at_neutral_ph(self) -> None:
        assert charge_at_ph("D", pH=7.0) < 0.0

    def test_unknown_nterminal_residue_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid amino acid: 1"):
            charge_at_ph("1ACDE", pH=7.0)

    def test_unknown_cterminal_residue_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid amino acid: 1"):
            charge_at_ph("ACDE1", pH=7.0)


class TestSecondaryStructure:
    def test_default_scale_returns_four_structure_types(self) -> None:
        result = secondary_structure("AAAA")
        assert set(result) == {"alpha_helix", "beta_sheet", "beta_turn", "coil"}

    def test_values_are_normalized_to_sum_to_one(self) -> None:
        result = secondary_structure("PEPTIDE")
        assert sum(result.values()) == pytest.approx(1.0)

    def test_levitt_scale_has_no_coil_component(self) -> None:
        result = secondary_structure("AAAA", scale="Levitt")
        assert set(result) == {"alpha_helix", "beta_sheet", "beta_turn"}
        assert sum(result.values()) == pytest.approx(1.0)


class TestGeneratePartitionsValidation:
    scale = {"A": 1.0, "C": 2.0, "D": 3.0, "E": 4.0}

    def test_non_positive_num_windows_raises(self) -> None:
        with pytest.raises(ValueError, match="num_windows must be positive"):
            generate_partitions("ACDE", scale=self.scale, num_windows=0)

    def test_negative_overlap_raises(self) -> None:
        with pytest.raises(ValueError, match="aa_overlap cannot be negative"):
            generate_partitions("ACDE", scale=self.scale, num_windows=2, aa_overlap=-1)

    def test_empty_sequence_raises(self) -> None:
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            generate_partitions("", scale=self.scale, num_windows=2)

    def test_too_many_non_overlapping_windows_raises(self) -> None:
        with pytest.raises(ValueError, match="Cannot create 5 non-overlapping windows"):
            generate_partitions("AC", scale=self.scale, num_windows=5, aa_overlap=0)


class TestGeneratePartitionsSingleWindow:
    def test_single_window_covers_the_whole_sequence(self) -> None:
        scale = {"A": 1.0, "C": 3.0, "D": 3.86, "E": 4.0}
        assert generate_partitions("ACDE", scale=scale, num_windows=1) == [calc_property("ACDE", scale=scale)]


class TestGeneratePartitionsShortSequence:
    """seq_len < num_windows forces heavily-overlapping (or repeated) windows."""

    scale = {"A": 1.0, "C": 2.0, "D": 3.0, "E": 4.0}

    def test_single_residue_sequence_repeats_same_value_for_every_window(self) -> None:
        result = generate_partitions("A", scale=self.scale, num_windows=3, aa_overlap=1)
        assert result == pytest.approx([1.0, 1.0, 1.0])

    def test_short_sequence_produces_num_windows_values(self) -> None:
        result = generate_partitions("AC", scale=self.scale, num_windows=4, aa_overlap=1)
        assert len(result) == 4
        assert result == pytest.approx([1.0, 1.0, 1.5, 1.5])


class TestGeneratePartitionsNormalCases:
    scale = {"A": 1.0, "C": 3.0, "D": 3.86, "E": 4.0}

    def test_two_windows_are_placed_at_start_and_end(self) -> None:
        result = generate_partitions("ACDEACDEAC", scale=self.scale, num_windows=2, aa_overlap=2)
        assert len(result) == 2
        assert result == pytest.approx([2.643333333333333, 2.643333333333333])

    def test_windows_are_centered_when_they_fit_with_room_to_spare(self) -> None:
        result = generate_partitions("ACDEACDEACDEACDE", scale=self.scale, num_windows=3, aa_overlap=0)
        assert len(result) == 3
        assert result == pytest.approx([2.643333333333333, 3.2866666666666666, 3.2866666666666666])

    def test_oversized_window_size_is_clamped_to_sequence_length(self) -> None:
        # A large aa_overlap relative to a short sequence drives the ideal window_size
        # up to seq_len, which must be clamped back down. With window_size == seq_len
        # there is only one valid window position, so every window covers the full
        # sequence and all three values are identical.
        result = generate_partitions("ACDEA", scale=self.scale, num_windows=3, aa_overlap=10)
        assert len(result) == 3
        assert result == pytest.approx([2.572, 2.572, 2.572])

    def test_num_windows_matches_output_length_across_overlap_values(self) -> None:
        for overlap in range(0, 4):
            result = generate_partitions("ACDEACDEACDE", scale=self.scale, num_windows=4, aa_overlap=overlap)
            assert len(result) == 4


class TestGeneratePartitionsUnevenSpacingEdgeCase:
    """Regression test for a fixed bug in uneven window spacing.

    Previously, when the requested step size marched an interior window's start all
    the way to (or past) the end of the sequence before the final window, the "ensure
    valid window" fallback (`end = min(seq_len, start + 1)`) was a no-op because
    `start` already equalled `seq_len`, producing an empty window; `calc_property`
    then silently returned 0.0 for it instead of a real property value. Every
    window's start is now clamped to `max(0, seq_len - window_size)` (not just the
    last window's), so interior windows can no longer run past the sequence end.
    """

    def test_middle_window_no_longer_collapses_to_zero(self) -> None:
        scale = {"A": 1.0, "C": 3.0, "D": 3.86, "E": 4.0}
        result = generate_partitions("ACDEAC", scale=scale, num_windows=5, aa_overlap=0)
        assert len(result) == 5
        assert result == pytest.approx([2.0, 3.93, 2.0, 2.0, 2.0])

"""Tests for window weighting schemes (peptacular.property.weights)."""

import pytest

from peptacular.property.weights import WeightingMethods, get_weights


class TestWeightShapes:
    def test_uniform_all_equal(self):
        w = get_weights(5, "uniform", min_weight=0.1, max_weight=1.0)
        assert list(w) == [1.0, 1.0, 1.0, 1.0, 1.0]

    def test_linear_is_center_peaked_and_symmetric(self):
        w = list(get_weights(5, "linear", min_weight=0.1, max_weight=1.0))
        assert w[2] == pytest.approx(1.0)  # center is the peak
        assert w[0] == pytest.approx(0.1)  # edges are the minimum
        assert w[-1] == pytest.approx(0.1)
        assert w == list(reversed(w))  # symmetric

    def test_gaussian_center_peaked_and_symmetric(self):
        w = list(get_weights(7, "gaussian", min_weight=0.1, max_weight=1.0))
        assert w[3] == pytest.approx(1.0)
        assert w == pytest.approx(list(reversed(w)))
        assert max(w) == pytest.approx(1.0)

    def test_all_schemes_in_range(self):
        for method in WeightingMethods:
            w = list(get_weights(9, method.value, min_weight=0.2, max_weight=0.9))
            assert all(0.0 <= x <= 1.0001 for x in w), f"{method.value}: {w}"
            assert len(w) == 9


class TestWeightEdgeCases:
    def test_length_one(self):
        for method in WeightingMethods:
            w = list(get_weights(1, method.value))
            assert len(w) == 1

    def test_length_zero_is_empty(self):
        assert list(get_weights(0, "linear")) == []

    def test_length_zero_is_empty_for_every_scheme(self):
        for method in WeightingMethods:
            assert list(get_weights(0, method.value)) == []

    def test_non_sequence_weights_raises_type_error(self):
        with pytest.raises(TypeError, match="weights must be a sequence"):
            get_weights(5, weights=123)  # type: ignore[arg-type]

    def test_mismatched_weight_list_length_raises_value_error(self):
        with pytest.raises(ValueError, match="does not match sequence length"):
            get_weights(5, weights=[1, 2, 3])

    def test_length_two_linear(self):
        # Documented current behavior: a 2-wide window is uniform (no interior peak).
        assert list(get_weights(2, "linear")) == [1.0, 1.0]

    def test_string_and_enum_agree(self):
        assert list(get_weights(6, "gaussian")) == list(get_weights(6, WeightingMethods.GAUSSIAN))

    def test_explicit_weight_list_passed_through(self):
        custom = [0.5, 0.6, 0.7]
        assert list(get_weights(3, custom)) == custom  # type: ignore[arg-type]

    def test_unknown_scheme_raises(self):
        with pytest.raises((ValueError, KeyError)):
            get_weights(5, "not_a_real_scheme")


class TestWindowPropertyIntegration:
    def test_weighted_window_property_runs(self):
        import peptacular as pt

        # Exercise the weighting path end-to-end through the public property API.
        vals = pt.hydrophobicity("PEPTIDEKMWLV")
        assert isinstance(vals, float)

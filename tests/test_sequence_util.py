"""Tests for peptacular.sequence.util (shared helpers for the functional API)."""

import pytest

import peptacular as pt
from peptacular.sequence.util import get_annotation_input, is_sequence_valid, round_to_precision, sequence_to_annotation


class TestRoundToPrecision:
    def test_rounds_when_precision_given(self):
        assert round_to_precision(3.14159, 2) == 3.14

    def test_passthrough_when_precision_none(self):
        assert round_to_precision(3.14159, None) == 3.14159

    def test_default_precision_is_none(self):
        assert round_to_precision(3.14159) == 3.14159

    def test_zero_precision(self):
        assert round_to_precision(3.6, 0) == 4.0


class TestGetAnnotationInput:
    def test_string_input_parses(self):
        result = get_annotation_input("PEPTIDE")
        assert isinstance(result, pt.ProFormaAnnotation)
        assert result.serialize() == "PEPTIDE"

    def test_annotation_input_copy_true_returns_new_object(self):
        a = pt.parse("PEPTIDE")
        result = get_annotation_input(a, copy=True)
        assert result is not a
        assert result == a

    def test_annotation_input_copy_false_returns_same_object(self):
        a = pt.parse("PEPTIDE")
        result = get_annotation_input(a, copy=False)
        assert result is a

    def test_default_copy_is_true(self):
        a = pt.parse("PEPTIDE")
        assert get_annotation_input(a) is not a

    def test_invalid_type_raises_type_error(self):
        with pytest.raises(TypeError, match="got int: 12345"):
            get_annotation_input(12345)  # type: ignore[arg-type]


class TestIsSequenceValid:
    def test_valid_string(self):
        assert is_sequence_valid("PEPTIDE") is True

    def test_invalid_string(self):
        assert is_sequence_valid("PEP[Oxidation") is False

    def test_annotation_input_is_always_valid(self):
        a = pt.parse("PEPTIDE")
        assert is_sequence_valid(a) is True


class TestSequenceToAnnotation:
    def test_parses_string(self):
        result = sequence_to_annotation("PEPTIDE")
        assert isinstance(result, pt.ProFormaAnnotation)
        assert result.serialize() == "PEPTIDE"

"""Coverage for peptacular.spans: Span helpers and digestion span builders."""

import pytest

from peptacular.spans import (
    Span,
    _grouped_left_semi_span_builder,
    _grouped_right_semi_span_builder,
    build_enzymatic_spans,
    build_left_semi_spans,
    build_non_enzymatic_spans,
    build_right_semi_spans,
    build_semi_spans,
    build_spans,
    calculate_span_coverage,
)


class TestSpan:
    def test_span_len(self):
        assert Span(0, 3, 0).span_len() == 3


class TestBuildNonEnzymaticSpans:
    def test_tuple_input_is_converted(self):
        result = list(build_non_enzymatic_spans((0, 3, 0)))
        assert result == [Span(0, 1, 0), Span(0, 2, 0), Span(0, 3, 0), Span(1, 2, 0), Span(1, 3, 0), Span(2, 3, 0)]

    def test_span_input(self):
        assert list(build_non_enzymatic_spans(Span(0, 3, 0))) == list(build_non_enzymatic_spans((0, 3, 0)))


class TestBuildLeftSemiSpans:
    def test_tuple_input_is_converted(self):
        result = list(build_left_semi_spans((0, 3, 0)))
        assert result == [Span(0, 2, 0), Span(0, 1, 0)]

    def test_span_input(self):
        assert list(build_left_semi_spans(Span(0, 3, 0))) == list(build_left_semi_spans((0, 3, 0)))


class TestBuildRightSemiSpans:
    def test_tuple_input_is_converted(self):
        result = list(build_right_semi_spans((0, 3, 0)))
        assert result == [Span(1, 3, 0), Span(2, 3, 0)]

    def test_span_input(self):
        assert list(build_right_semi_spans(Span(0, 3, 0))) == list(build_right_semi_spans((0, 3, 0)))


class TestBuildEnzymaticSpans:
    def test_default_min_max_len(self):
        result = list(build_enzymatic_spans(5, [3], 1))
        assert result == [Span(0, 3, 0), Span(0, 5, 1), Span(3, 5, 0)]

    def test_explicit_min_len(self):
        assert list(build_enzymatic_spans(5, [3], 1, min_len=5)) == [Span(0, 5, 1)]

    def test_explicit_max_len(self):
        assert list(build_enzymatic_spans(5, [3], 1, max_len=3)) == [Span(0, 3, 0), Span(3, 5, 0)]


class TestGroupedLeftSemiSpanBuilder:
    def test_default_min_len(self):
        result = list(_grouped_left_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)]))
        assert result == [Span(0, 4, 1), Span(0, 2, 0), Span(0, 1, 0), Span(3, 4, 0)]

    def test_span_len_at_or_below_min_len_breaks_group(self):
        # span_len(1) <= min_len(1) -> breaks before yielding anything for this group.
        assert list(_grouped_left_semi_span_builder([(0, 1, 0)], min_len=1)) == []


class TestGroupedRightSemiSpanBuilder:
    def test_default_min_len(self):
        result = list(_grouped_right_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)]))
        assert result == [Span(1, 3, 0), Span(2, 3, 0), Span(1, 5, 1), Span(2, 5, 1), Span(4, 5, 0)]

    def test_span_len_below_min_len_breaks_group(self):
        # span_len(1) < min_len(2) -> breaks before yielding anything for this group.
        assert list(_grouped_right_semi_span_builder([(1, 2, 0)], min_len=2)) == []


class TestBuildSemiSpans:
    def test_combines_left_and_right(self):
        result = list(build_semi_spans([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=1, max_len=5))
        assert result == [
            Span(0, 4, 1),
            Span(0, 2, 0),
            Span(0, 1, 0),
            Span(3, 4, 0),
            Span(1, 3, 0),
            Span(2, 3, 0),
            Span(1, 5, 1),
            Span(2, 5, 1),
            Span(4, 5, 0),
        ]


class TestBuildSpans:
    def test_non_enzymatic_case(self):
        # every position is an enzyme site -> falls back to pure non-enzymatic spans
        result = list(build_spans(5, [0, 1, 2, 3, 4, 5], 0))
        assert Span(0, 1, 0) in result
        assert Span(0, 5, 0) in result

    def test_enzymatic_non_semi(self):
        result = list(build_spans(5, [3], 1, semi=False))
        assert result == [Span(0, 3, 0), Span(0, 5, 1), Span(3, 5, 0)]

    def test_enzymatic_semi(self):
        result = list(build_spans(5, [3], 1, semi=True))
        assert Span(0, 3, 0) in result
        assert Span(1, 3, 0) in result  # a semi span not present in the non-semi output


class TestCalculateSpanCoverage:
    def test_non_overlapping(self):
        assert calculate_span_coverage([(0, 3, 0), (3, 6, 0), (6, 9, 0)], 9) == [1] * 9

    def test_overlapping_without_accumulate(self):
        assert calculate_span_coverage([(0, 3, 0), (0, 3, 0), (6, 9, 0)], 9) == [1, 1, 1, 0, 0, 0, 1, 1, 1]

    def test_overlapping_with_accumulate(self):
        assert calculate_span_coverage([(0, 3, 0), (0, 3, 0), (6, 9, 0)], 9, accumulate=True) == [2, 2, 2, 0, 0, 0, 1, 1, 1]

    def test_max_index_larger_than_spans(self):
        assert calculate_span_coverage([(0, 3, 0), (3, 6, 0), (6, 9, 0)], 12) == [1] * 9 + [0, 0, 0]

    def test_max_index_smaller_than_spans_raises(self):
        with pytest.raises(IndexError):
            calculate_span_coverage([(0, 3, 0), (3, 6, 0), (6, 9, 0)], 6)

    def test_span_object_input(self):
        assert calculate_span_coverage([Span(0, 3, 0)], 3) == [1, 1, 1]

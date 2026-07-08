"""Coverage for peptacular.property.types: string-to-enum helpers."""

import pytest

from peptacular.property.types import AggregationMethod, MissingAAHandling, WeightingMethods


class TestMissingAAHandling:
    def test_from_str_case_insensitive(self):
        assert MissingAAHandling.from_str("ZERO") is MissingAAHandling.ZERO

    def test_from_str_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown MissingAAHandling"):
            MissingAAHandling.from_str("bogus")


class TestAggregationMethod:
    def test_from_str_case_insensitive(self):
        assert AggregationMethod.from_str("SUM") is AggregationMethod.SUM

    def test_from_str_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown AggregationMethod"):
            AggregationMethod.from_str("bogus")


class TestWeightingMethods:
    def test_from_str_case_insensitive(self):
        assert WeightingMethods.from_str("UNIFORM") is WeightingMethods.UNIFORM

    def test_from_str_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown WeightingMethods"):
            WeightingMethods.from_str("bogus")

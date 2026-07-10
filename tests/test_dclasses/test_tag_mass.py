"""
Tests for parsing mass delta modification tags.
"""

import pytest

import peptacular as pt


class TestTagMass:
    """Tests for parsing mass delta modifications"""

    def test_positive_mass(self):
        """Test parsing positive mass delta"""
        result = pt.ModificationTags.from_string("+15.995").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(15.995)  # type: ignore
        assert result.cv is None

    def test_negative_mass(self):
        """Test parsing negative mass delta"""
        result = pt.ModificationTags.from_string("-18.010").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(-18.010)  # type: ignore

    def test_mass_without_sign_invalid(self):
        """Test that mass without sign is treated as a name (not valid ProForma)"""
        result = pt.ModificationTags.from_string("42.010").tags[0]
        # Without a sign, it's not recognized as a mass, so it becomes a TagName
        assert isinstance(result, pt.TagName)
        assert result.name == "42.010"

    def test_integer_mass(self):
        """Test parsing integer mass with required sign"""
        result = pt.ModificationTags.from_string("+16").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == 16.0

    def test_mass_with_cv_prefix(self):
        """Test parsing mass delta with single-letter CV prefix (Rule 2)"""
        result = pt.ModificationTags.from_string("U:+15.995").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(15.995)  # type: ignore
        assert result.cv == pt.CV.UNIMOD

        # Test other single-letter prefixes
        result2 = pt.ModificationTags.from_string("M:-18.010").tags[0]
        assert isinstance(result2, pt.TagMass)
        assert result2.mass == pytest.approx(-18.010)  # type: ignore
        assert result2.cv == pt.CV.PSI_MOD

    def test_observed_mass_basic(self):
        """Test parsing observed mass delta with Obs prefix"""
        result = pt.ModificationTags.from_string("Obs:+17.05685").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(17.05685)  # type: ignore
        assert result.cv == pt.CV.OBSERVED

    def test_observed_mass_negative(self):
        """Test parsing negative observed mass delta"""
        result = pt.ModificationTags.from_string("Obs:-18.010").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(-18.010)  # type: ignore
        assert result.cv == pt.CV.OBSERVED

    def test_observed_mass_case_insensitive(self):
        """Test that Obs prefix is case insensitive"""
        result = pt.ModificationTags.from_string("obs:+17.05685").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(17.05685)  # type: ignore
        assert result.cv == pt.CV.OBSERVED

        result2 = pt.ModificationTags.from_string("OBS:+17.05685").tags[0]
        assert isinstance(result2, pt.TagMass)
        assert result2.mass == pytest.approx(17.05685)  # type: ignore
        assert result2.cv == pt.CV.OBSERVED

    def test_custom_mass_prefix_round_trip(self):
        """The ``C:`` custom-mass prefix must survive a parse -> serialize round trip.

        Regression: ``C:`` was parsed to ``cv=None`` (indistinguishable from a bare
        mass), so it serialized back to ``+100`` and the custom marker was lost.
        """
        result = pt.ModificationTags.from_string("C:+100.0").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.mass == pytest.approx(100.0)  # type: ignore
        assert result.cv == pt.CV.CUSTOM
        assert str(result) == "C:+100"

    def test_bare_mass_still_has_no_cv(self):
        """A prefix-less mass must remain ``cv=None`` (not treated as custom)."""
        result = pt.ModificationTags.from_string("+100.0").tags[0]
        assert isinstance(result, pt.TagMass)
        assert result.cv is None
        assert str(result) == "+100"

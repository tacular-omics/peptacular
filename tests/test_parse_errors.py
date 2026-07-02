"""Tests that malformed ProForma input is rejected with clear, actionable errors."""

import pytest

import peptacular as pt


class TestUnclosedBrackets:
    """Regression: an unterminated modification bracket must not be silently completed."""

    @pytest.mark.parametrize("seq", ["PEP[Oxidation", "[Acetyl-PEPTIDE", "PEPTIDE-[Amide"])
    def test_unclosed_bracket_rejected(self, seq):
        with pytest.raises(ValueError, match="Unclosed"):
            pt.parse(seq)


class TestEmptyModifications:
    """Regression: empty ``[]`` modifications must be rejected at parse time, not lazily."""

    @pytest.mark.parametrize("seq", ["PEP[]TIDE", "[]-PEPTIDE", "PEPTIDE-[]"])
    def test_empty_modification_rejected(self, seq):
        with pytest.raises(ValueError, match="Empty modification"):
            pt.parse(seq)


class TestChargeSeparator:
    """Regression: a dangling '/' charge separator must be rejected."""

    def test_trailing_slash_rejected(self):
        with pytest.raises(ValueError, match="Expected a charge value after"):
            pt.parse("PEPTIDE/")

    @pytest.mark.parametrize("seq", ["PEPTIDE/abc", "PEPTIDE/x"])
    def test_non_numeric_charge_rejected(self, seq):
        with pytest.raises(ValueError, match="Invalid charge after"):
            pt.parse(seq)


class TestActionableMessages:
    """Errors should name the offending value and what is expected (useful for agents)."""

    def test_lowercase_residue_message(self):
        with pytest.raises(ValueError, match="must be uppercase"):
            pt.parse("peptide")

    def test_unknown_mod_name_message(self):
        # names the value and gives a hint on how to specify a valid modification
        with pytest.raises(ValueError, match="Unknown modification name 'Notarealmod'.*Specify one of"):
            pt.mass("PEP[Notarealmod]TIDE")

    def test_unknown_accession_message(self):
        with pytest.raises(ValueError, match="Unknown modification accession 'UNIMOD:9999999'"):
            pt.mass("PEP[UNIMOD:9999999]TIDE")

    def test_wrong_type_reports_type_and_value(self):
        with pytest.raises(TypeError, match=r"got int: 12345"):
            pt.mass(12345)  # type: ignore[arg-type]


class TestStillValid:
    """Guard against over-eager rejection of valid input."""

    @pytest.mark.parametrize(
        "seq",
        [
            "PEP[Oxidation]TIDE",
            "PEPTIDE/2",
            "PEPTIDE/-1",
            "[Acetyl]-PEPTIDE-[Amidated]",
            "PEM[Oxidation]TIDE",
        ],
    )
    def test_valid_round_trips(self, seq):
        assert pt.parse(seq).serialize() == seq

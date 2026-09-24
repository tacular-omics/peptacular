"""Typed errors on the main public paths (4.2). Every class is still a ``ValueError``."""

import warnings

import pytest

import peptacular as pt

ERROR_CLASSES = [
    pt.ProFormaFormatError,
    pt.UnknownModificationError,
    pt.CompositionError,
    pt.InvalidAdjustmentError,
    pt.UnsupportedOperationError,
    pt.InvalidPositionError,
    pt.FastaFormatError,
]


@pytest.mark.parametrize("cls", ERROR_CLASSES)
def test_errors_share_peptacular_base_and_stay_value_errors(cls):
    assert issubclass(cls, pt.PeptacularError)
    assert issubclass(cls, ValueError)


def test_peptacular_error_is_value_error():
    assert issubclass(pt.PeptacularError, ValueError)


def test_parse_error_caught_by_base():
    with pytest.raises(pt.PeptacularError):
        pt.parse("PEP[")


@pytest.mark.parametrize("n", ["a", 1.5, None])
def test_shift_rejects_non_integer(n):
    with pytest.raises(TypeError, match="n must be an int"):
        pt.shift("PEPTIDE", n)


def test_shift_still_accepts_int():
    assert pt.shift("PEPTIDE", 2) == "PTIDEPE"


def test_mass_of_empty_sequence_is_composition_error():
    with pytest.raises(pt.CompositionError, match="empty sequence"):
        pt.mass("")


def test_unknown_ion_type_lists_valid_types():
    with pytest.raises(pt.UnsupportedOperationError, match=r"Unknown ion type 'q'.*'b'.*'y'"):
        pt.fragment("PEPTIDE", ion_types=["q"])


def test_unknown_ion_type_on_frag_and_mass():
    annot = pt.parse("PEPTIDE")
    with pytest.raises(pt.UnsupportedOperationError, match="Unknown ion type"):
        annot.frag(ion_type="q")
    with pytest.raises(pt.UnsupportedOperationError, match="Unknown ion type"):
        pt.mass("PEPTIDE", ion_type="q")


def test_out_of_range_slice_is_position_error():
    annot = pt.parse("PEPTIDE")
    with pytest.raises(pt.InvalidPositionError, match="exceeds the sequence length"):
        annot[2:40]
    with pytest.raises(pt.InvalidPositionError):
        annot[-40:]


def test_out_of_range_frag_position_is_position_error():
    with pytest.raises(pt.InvalidPositionError, match="position"):
        pt.parse("PEPTIDE").frag(ion_type="b", position=40)


@pytest.mark.parametrize(
    ("text", "match"),
    [
        ("garbage", "Sequence data before header"),
        (">x\n", "No valid FASTA sequences"),
        ("", "Empty input"),
        (">\nPEP", "Empty header"),
    ],
)
def test_parse_fasta_text_errors_are_typed(text, match):
    with pytest.raises(pt.FastaFormatError, match=match):
        pt.parse_fasta_text(text)


def test_digest_warns_on_unknown_plain_enzyme_name():
    with pytest.warns(UserWarning, match="notanenzyme"):
        result = pt.digest("PEPTIDEK", "notanenzyme")
    assert [seq for seq, *_ in result] == ["PEPTIDEK"]


@pytest.mark.parametrize("enzyme", ["trypsin", "Trypsin", "([KR])", "(?<=K)"])
def test_digest_known_protease_or_regex_does_not_warn(enzyme):
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        warnings.filterwarnings("ignore", message="The regex pattern has a non-zero-length match")
        pt.digest("PEPTIDEKAAR", enzyme)

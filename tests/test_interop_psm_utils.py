import pytest

from peptacular import ProFormaAnnotation
from peptacular.interop.psm_utils import from_psm_utils, to_psm_utils

psm_utils = pytest.importorskip("psm_utils")


@pytest.mark.parametrize("sequence", ["PEPTIDE", "[Acetyl]-PEM[Oxidation]TIDE/2"])
def test_psm_utils_round_trip(sequence):
    annotation = ProFormaAnnotation.parse(sequence)

    converted = to_psm_utils(annotation)

    assert isinstance(converted, psm_utils.Peptidoform)
    assert from_psm_utils(converted) == annotation


def test_psm_utils_rejects_wrong_input_type():
    with pytest.raises(TypeError, match="Expected psm_utils.Peptidoform"):
        from_psm_utils("PEPTIDE")


def test_psm_utils_negative_charge_never_silently_changes():
    from peptacular.interop import InteropConversionError

    original = ProFormaAnnotation.parse("PEPTIDE/-2")
    try:
        converted = to_psm_utils(original)
    except InteropConversionError as exc:
        assert "round trip" in str(exc)
    else:
        assert from_psm_utils(converted).to_dict() == original.to_dict()

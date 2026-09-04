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

import pytest

pytest.importorskip("pyteomics")

from pyteomics import mass, proforma

from peptacular import ChargedFormula, ProFormaAnnotation
from peptacular.interop import InteropConversionError
from peptacular.interop.pyteomics import (
    from_pyteomics,
    from_pyteomics_composition,
    to_pyteomics,
    to_pyteomics_composition,
)


@pytest.mark.parametrize(
    "sequence",
    [
        "PEPTIDE",
        "[Acetyl]-PEM[Oxidation]TIDE/2",
        "PEPTIDE/[H:z+1,Na:z+1]",
    ],
)
def test_pyteomics_annotation_round_trip(sequence):
    annotation = ProFormaAnnotation.parse(sequence)

    converted = to_pyteomics(annotation)

    assert isinstance(converted, proforma.ProForma)
    assert from_pyteomics(converted) == annotation


def test_pyteomics_rejects_wrong_input_type():
    with pytest.raises(TypeError, match="Expected pyteomics.proforma.ProForma"):
        from_pyteomics("PEPTIDE")


def test_pyteomics_composition_round_trip_with_isotope():
    formula = ChargedFormula.from_composition({"C": 4, "13C": 2, "H": 7})

    converted = to_pyteomics_composition(formula)

    assert isinstance(converted, mass.Composition)
    assert dict(converted) == {"C": 4, "C[13]": 2, "H": 7}
    assert from_pyteomics_composition(converted).get_dict_composition() == formula.get_dict_composition()


def test_pyteomics_composition_rejects_charge_metadata():
    formula = ChargedFormula.from_composition({"H": 1}, charge=1)

    with pytest.raises(InteropConversionError, match="charge metadata"):
        to_pyteomics_composition(formula)


def test_pyteomics_composition_rejects_special_keys():
    with pytest.raises(InteropConversionError, match="not representable"):
        from_pyteomics_composition({"H+": 1})

"""Tests for ProtParam property calculations: instability index, aliphatic index, and molar extinction coefficient."""

import pytest

import peptacular as pt
from peptacular.diagnostics import PeptacularError
from peptacular.property.core import (
    ALIPHATIC_INDEX_ILE_LEU_FACTOR,
    ALIPHATIC_INDEX_SCALE,
    ALIPHATIC_INDEX_VAL_FACTOR,
    EXTINCTION_COEFFICIENT_CYSTINE,
    EXTINCTION_COEFFICIENT_TRP,
    EXTINCTION_COEFFICIENT_TYR,
    INSTABILITY_INDEX_SCALE,
    aliphatic_index,
    instability_index,
    molar_extinction_coefficient,
)


def test_protparam_constants():
    assert INSTABILITY_INDEX_SCALE == 10.0
    assert ALIPHATIC_INDEX_VAL_FACTOR == 2.9
    assert ALIPHATIC_INDEX_ILE_LEU_FACTOR == 3.9
    assert ALIPHATIC_INDEX_SCALE == 100.0
    assert EXTINCTION_COEFFICIENT_TRP == 5500
    assert EXTINCTION_COEFFICIENT_TYR == 1490
    assert EXTINCTION_COEFFICIENT_CYSTINE == 125


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("PEPTIDE", 60.9143),
        ("SAMPLER", 71.3429),
        ("MKWVTFISLLLLFSSAYS", 17.5667),
        ("AVILAVIL", 56.9),
        ("ACDEFGHIKLMNPQRSTVWY", 84.74),
    ],
)
def test_instability_index_biopython_benchmark_values(sequence: str, expected: float):
    assert instability_index(sequence) == pytest.approx(expected, abs=1e-4)
    assert pt.instability_index(sequence) == pytest.approx(expected, abs=1e-4)


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("PEPTIDE", 55.7143),
        ("AVILAVIL", 292.5),
        ("MKWVTFISLLLLFSSAYS", 130.0),
        ("ACDEFGHIKLMNPQRSTVWY", 58.5),
    ],
)
def test_aliphatic_index_biopython_benchmark_values(sequence: str, expected: float):
    assert aliphatic_index(sequence) == pytest.approx(expected, abs=1e-4)
    assert pt.aliphatic_index(sequence) == pytest.approx(expected, abs=1e-4)


def test_molar_extinction_coefficient_biopython_benchmark_values():
    assert molar_extinction_coefficient("WYCCWY", cystines=False) == 13980
    assert molar_extinction_coefficient("WYCCWY", cystines=True) == 14105
    assert molar_extinction_coefficient("PEPTIDE", cystines=False) == 0
    assert molar_extinction_coefficient("PEPTIDE", cystines=True) == 0

    assert pt.molar_extinction_coefficient("WYCCWY") == 13980
    assert pt.molar_extinction_coefficient("WYCCWY", cystines=True) == 14105
    assert pt.molar_extinction_coefficient("PEPTIDE") == 0


def test_edge_cases_empty_and_short_sequences():
    assert instability_index("") == 0.0
    assert instability_index("A") == 0.0
    assert pt.instability_index("") == 0.0
    assert pt.instability_index("M") == 0.0

    assert aliphatic_index("") == 0.0
    assert pt.aliphatic_index("") == 0.0

    assert molar_extinction_coefficient("") == 0
    assert pt.molar_extinction_coefficient("") == 0


def test_instability_index_unknown_dipeptide_raises():
    with pytest.raises(PeptacularError, match="PX"):
        instability_index("PEPXIDE")

    with pytest.raises(PeptacularError, match="PX"):
        pt.instability_index("PEPXIDE")

    with pytest.raises(PeptacularError, match="XP"):
        instability_index("XPEPTIDE")


def test_annotation_properties_stripped_sequence():
    annot = pt.parse("PEPTIDE")
    assert annot.prop.instability_index == pytest.approx(60.9143, abs=1e-4)
    assert annot.prop.aliphatic_index == pytest.approx(55.7143, abs=1e-4)
    assert annot.prop.molar_extinction_coefficient() == 0

    mod_annot = pt.parse("[Acetyl]-PEP[Phospho]TIDE-[Amidated]/2")
    assert mod_annot.prop.instability_index == pytest.approx(60.9143, abs=1e-4)
    assert mod_annot.prop.aliphatic_index == pytest.approx(55.7143, abs=1e-4)
    assert mod_annot.prop.molar_extinction_coefficient() == 0

    wy_annot = pt.parse("W[Phospho]YCCWY/3")
    assert wy_annot.prop.molar_extinction_coefficient(cystines=False) == 13980
    assert wy_annot.prop.molar_extinction_coefficient(cystines=True) == 14105


def test_sequence_wrapper_accepts_annotation_and_list():
    annot = pt.parse("MKWVTFISLLLLFSSAYS")
    assert pt.instability_index(annot) == pytest.approx(17.5667, abs=1e-4)
    assert pt.aliphatic_index(annot) == pytest.approx(130.0, abs=1e-4)
    assert pt.molar_extinction_coefficient(annot) == 5500 + 1490

    seqs = ["PEPTIDE", "SAMPLER", "AVILAVIL"]
    ii_list = pt.instability_index(seqs, n_workers=1)
    assert ii_list == pytest.approx([60.9143, 71.3429, 56.9], abs=1e-4)

    ai_list = pt.aliphatic_index(seqs, n_workers=1)
    assert ai_list == pytest.approx([55.7143, 70.0, 292.5], abs=1e-4)

    ec_list = pt.molar_extinction_coefficient(["WYCCWY", "PEPTIDE"], cystines=True, n_workers=1)
    assert ec_list == [14105, 0]

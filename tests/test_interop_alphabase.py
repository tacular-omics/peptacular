import pandas as pd
import pytest

from peptacular import ProFormaAnnotation
from peptacular.interop import (
    InteropConversionError,
    LossPolicy,
    LossyConversionWarning,
    from_alphabase_dataframe,
    from_alphabase_row,
    to_alphabase_dataframe,
    to_alphabase_row,
)

alphabase_modification = pytest.importorskip("alphabase.constants.modification")
alphabase_spectral_library = pytest.importorskip("alphabase.spectral_library.base")


def test_alphabase_localized_modification_row_round_trip():
    annotation = ProFormaAnnotation.parse("[Acetyl]-PEM[Oxidation]TIDE-[Amidated]/2")

    row = to_alphabase_row(annotation)

    assert row == {
        "sequence": "PEMTIDE",
        "mods": "Acetyl@Any_N-term;Oxidation@M;Amidated@Any_C-term",
        "mod_sites": "0;3;-1",
        "charge": 2,
    }
    assert from_alphabase_row(row) == annotation

    masses = alphabase_modification.calc_modification_mass(
        len(row["sequence"]),
        row["mods"].split(";"),
        [int(site) for site in row["mod_sites"].split(";")],
    )
    assert masses[0] != 0
    assert masses[2] != 0
    assert masses[-1] != 0


def test_alphabase_dataframe_is_native_and_refined():
    annotations = [
        ProFormaAnnotation.parse("PEPTIDE/2"),
        ProFormaAnnotation.parse("AC[Carbamidomethyl]DE/3"),
    ]

    dataframe = to_alphabase_dataframe(annotations)

    assert isinstance(dataframe, pd.DataFrame)
    assert list(dataframe.columns) == ["sequence", "mods", "mod_sites", "charge", "nAA"]
    assert dataframe["nAA"].tolist() == [4, 7]
    assert from_alphabase_dataframe(dataframe) == [annotations[1], annotations[0]]

    library = alphabase_spectral_library.SpecLibBase()
    library.precursor_df = dataframe.copy()
    pd.testing.assert_frame_equal(library.precursor_df, dataframe)


def test_alphabase_uncharged_dataframe_omits_charge_column():
    dataframe = to_alphabase_dataframe([ProFormaAnnotation.parse("PEPTIDE")])

    assert "charge" not in dataframe.columns
    assert from_alphabase_dataframe(dataframe) == [ProFormaAnnotation.parse("PEPTIDE")]


def test_alphabase_dataframe_rejects_mixed_charge_presence():
    with pytest.raises(InteropConversionError, match="cannot mix charged and uncharged"):
        to_alphabase_dataframe([ProFormaAnnotation.parse("PEPTIDE"), ProFormaAnnotation.parse("PEPTIDE/2")])


def test_alphabase_expands_fixed_modifications():
    annotation = ProFormaAnnotation.parse("<[Carbamidomethyl]@C>ACDC")

    row = to_alphabase_row(annotation)

    assert row["mods"] == "Carbamidomethyl@C;Carbamidomethyl@C"
    assert row["mod_sites"] == "2;4"
    assert from_alphabase_row(row).serialize() == "AC[Carbamidomethyl]DC[Carbamidomethyl]"


def test_alphabase_rejects_unsupported_features_by_default():
    annotation = ProFormaAnnotation.parse("{Oxidation}PEPTIDE")

    with pytest.raises(InteropConversionError, match="labile modifications"):
        to_alphabase_row(annotation)


def test_alphabase_rejects_cross_link_modifications():
    annotation = ProFormaAnnotation.parse("PEPTK[XLMOD:02001#XL1]IDE")

    with pytest.raises(InteropConversionError, match="XLMOD:02001#XL1"):
        to_alphabase_row(annotation)


def test_alphabase_warn_policy_discards_unsupported_features():
    annotation = ProFormaAnnotation.parse("{Oxidation}PEPTIDE")

    with pytest.warns(LossyConversionWarning, match="labile modifications"):
        row = to_alphabase_row(annotation, loss_policy=LossPolicy.WARN)

    assert row == {"sequence": "PEPTIDE", "mods": "", "mod_sites": "", "charge": None}


def test_alphabase_drop_policy_is_explicit_and_quiet():
    annotation = ProFormaAnnotation.parse("{Oxidation}PEPTIDE")

    row = to_alphabase_row(annotation, loss_policy="drop")

    assert row == {"sequence": "PEPTIDE", "mods": "", "mod_sites": "", "charge": None}


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ({}, "missing required 'sequence'"),
        ({"sequence": "PEPTIDE", "mods": "Oxidation@M", "mod_sites": ""}, "equal lengths"),
        ({"sequence": "PEPTIDE", "mods": "Oxidation@M", "mod_sites": "x"}, "Invalid AlphaBase modification site"),
        ({"sequence": "PEPTIDE", "mods": "Oxidation", "mod_sites": "1"}, "expected 'Name@Target'"),
        ({"sequence": "PEPTIDE", "mods": "Oxidation@M", "mod_sites": "1"}, "residue 1 is 'P'"),
        ({"sequence": "PEPTIDE", "mods": "Oxidation@M", "mod_sites": "20"}, "outside sequence length"),
        ({"sequence": "PEPTIDE", "mods": "", "mod_sites": "", "charge": "bad"}, "Invalid AlphaBase charge"),
    ],
)
def test_alphabase_rejects_invalid_rows(row, message):
    with pytest.raises(InteropConversionError, match=message):
        from_alphabase_row(row)


def test_alphabase_dataframe_rejects_missing_columns():
    with pytest.raises(InteropConversionError, match="mod_sites, mods"):
        from_alphabase_dataframe(pd.DataFrame({"sequence": ["PEPTIDE"]}))


def test_alphabase_dataframe_rejects_wrong_type():
    with pytest.raises(TypeError, match="Expected pandas.DataFrame"):
        from_alphabase_dataframe([])

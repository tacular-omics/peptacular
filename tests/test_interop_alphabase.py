import pytest

from peptacular import ProFormaAnnotation
from peptacular.interop import (
    AlphaBasePeptide,
    InteropConversionError,
    LossPolicy,
    LossyConversionWarning,
    from_alphabase,
    to_alphabase,
)

alphabase_modification = pytest.importorskip("alphabase.constants.modification")


def test_alphabase_localized_modification_round_trip():
    annotation = ProFormaAnnotation.parse("[Acetyl]-PEM[Oxidation]TIDE-[Amidated]/2")

    converted = to_alphabase(annotation)

    assert converted == AlphaBasePeptide(
        sequence="PEMTIDE",
        mods="Acetyl@Any_N-term;Oxidation@M;Amidated@Any_C-term",
        mod_sites="0;3;-1",
        charge=2,
    )
    assert from_alphabase(**converted.as_dict()) == annotation

    masses = alphabase_modification.calc_modification_mass(
        len(converted.sequence),
        converted.mods.split(";"),
        [int(site) for site in converted.mod_sites.split(";")],
    )
    assert masses[0] != 0
    assert masses[2] != 0
    assert masses[-1] != 0


def test_alphabase_expands_fixed_modifications():
    annotation = ProFormaAnnotation.parse("<[Carbamidomethyl]@C>ACDC")

    converted = to_alphabase(annotation)

    assert converted.mods == "Carbamidomethyl@C;Carbamidomethyl@C"
    assert converted.mod_sites == "2;4"
    assert from_alphabase(**converted.as_dict()).serialize() == "AC[Carbamidomethyl]DC[Carbamidomethyl]"


def test_alphabase_rejects_unsupported_features_by_default():
    annotation = ProFormaAnnotation.parse("{Oxidation}PEPTIDE")

    with pytest.raises(InteropConversionError, match="labile modifications"):
        to_alphabase(annotation)


def test_alphabase_rejects_cross_link_modifications():
    annotation = ProFormaAnnotation.parse("PEPTK[XLMOD:02001#XL1]IDE")

    with pytest.raises(InteropConversionError, match="XLMOD:02001#XL1"):
        to_alphabase(annotation)


def test_alphabase_warn_policy_discards_unsupported_features():
    annotation = ProFormaAnnotation.parse("{Oxidation}PEPTIDE")

    with pytest.warns(LossyConversionWarning, match="labile modifications"):
        converted = to_alphabase(annotation, loss_policy=LossPolicy.WARN)

    assert converted == AlphaBasePeptide(sequence="PEPTIDE")


def test_alphabase_drop_policy_is_explicit_and_quiet():
    annotation = ProFormaAnnotation.parse("{Oxidation}PEPTIDE")

    converted = to_alphabase(annotation, loss_policy="drop")

    assert converted == AlphaBasePeptide(sequence="PEPTIDE")


@pytest.mark.parametrize(
    ("mods", "sites", "message"),
    [
        ("Oxidation@M", "", "equal lengths"),
        ("Oxidation@M", "x", "Invalid AlphaBase modification site"),
        ("Oxidation", "1", "expected 'Name@Target'"),
        ("Oxidation@M", "1", "residue 1 is 'P'"),
        ("Oxidation@M", "20", "outside sequence length"),
    ],
)
def test_alphabase_rejects_invalid_columns(mods, sites, message):
    with pytest.raises(InteropConversionError, match=message):
        from_alphabase("PEPTIDE", mods, sites)

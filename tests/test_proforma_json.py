import json

import pytest
from tacular import AminoAcid, Element

from peptacular import (
    PROFORMA_JSON_SCHEMA_ID,
    PROFORMA_JSON_SCHEMA_VERSION,
    ChargedFormula,
    CompoundPeptidoformIon,
    FormulaElement,
    GlobalChargeCarrier,
    IsotopeReplacement,
    ModificationAmbiguousPrimary,
    ModificationAmbiguousSecondary,
    ModificationCrossLinker,
    ModificationTags,
    Peptidoform,
    PeptidoformIon,
    ProFormaAnnotation,
    SequenceElement,
    TagName,
    from_proforma_dict,
    get_proforma_json_schema,
)
from peptacular.annotation import Interval


def test_annotation_json_round_trip_preserves_all_state():
    annotation = ProFormaAnnotation(
        sequence="PEPTIDE",
        compound_name="compound",
        ion_name="ion",
        peptide_name="peptide",
        isotope_mods={"13C": 1},
        static_mods={"[Carbamidomethyl]@C": 1},
        labile_mods={"Glycan:Hex": 2},
        unknown_mods={"Phospho#g1": 1},
        nterm_mods={"Acetyl": 1},
        cterm_mods={"Amidated": 1},
        internal_mods={1: {"Oxidation": 2}, 5: {"#g1(0.25)": 1}},
        intervals=[Interval(1, 4, ambiguous=True, mods={"Phospho#g1(0.75)": 1})],
        charge=["Na:z+1", "H:z+1"],
    )

    encoded = annotation.to_dict()
    restored = ProFormaAnnotation.from_dict(encoded)

    assert encoded["schema_version"] == PROFORMA_JSON_SCHEMA_VERSION
    assert encoded["$schema"] == PROFORMA_JSON_SCHEMA_ID
    assert restored.to_dict() == encoded
    assert ProFormaAnnotation.from_json(annotation.to_json()).to_dict() == encoded
    json.dumps(encoded, allow_nan=False)


def test_compound_component_json_round_trip_is_lossless():
    phospho = ModificationTags((TagName("Phospho"),))
    crosslink = ModificationTags((TagName("Disulfide"),))
    primary = ModificationAmbiguousPrimary(label="g1", tags=phospho, score=0.7)
    secondary = ModificationAmbiguousSecondary(label="g1", score=0.3)
    first = Peptidoform(
        sequence=(
            SequenceElement(AminoAcid("P"), (primary,)),
            SequenceElement(AminoAcid("C"), (ModificationCrossLinker("XL1", crosslink),)),
        ),
        name="first",
    )
    second = Peptidoform(
        sequence=(
            SequenceElement(AminoAcid("C"), (ModificationCrossLinker("XL1"),)),
            SequenceElement(AminoAcid("E"), (secondary,)),
        ),
        name="second",
    )
    proton = GlobalChargeCarrier(ChargedFormula((FormulaElement(Element("H"), 1),), charge=1), occurance=2)
    compound = CompoundPeptidoformIon(
        peptidoform_ions=(PeptidoformIon((first, second), name="crosslinked", charge=(proton,)),),
        name="compound",
        isotope_replacement=(IsotopeReplacement(Element("C"), 13),),
    )

    encoded = compound.to_dict()
    restored = CompoundPeptidoformIon.from_json(compound.to_json(indent=2))

    assert restored == compound
    assert restored.to_dict() == encoded


def test_decoder_rejects_unknown_versions_fields_and_types():
    encoded = ProFormaAnnotation(sequence="PEPTIDE").to_dict()

    wrong_version = dict(encoded, schema_version="2.0")
    with pytest.raises(ValueError, match="Unsupported ProForma JSON schema version"):
        ProFormaAnnotation.from_dict(wrong_version)

    unknown_field = dict(encoded, surprise=True)
    with pytest.raises(ValueError, match="Unknown JSON field"):
        ProFormaAnnotation.from_dict(unknown_field)

    wrong_type = dict(encoded, **{"$type": "Peptidoform"})
    with pytest.raises(ValueError, match="Missing required JSON field"):
        from_proforma_dict(wrong_type)


def test_classmethod_enforces_expected_root_type():
    encoded = ProFormaAnnotation(sequence="PEPTIDE").to_dict()
    with pytest.raises(TypeError, match="Expected Peptidoform"):
        Peptidoform.from_dict(encoded)


def test_bundled_schema_matches_encoder_metadata():
    schema = get_proforma_json_schema()
    assert schema["$id"] == PROFORMA_JSON_SCHEMA_ID
    assert schema["properties"]["schema_version"]["const"] == PROFORMA_JSON_SCHEMA_VERSION

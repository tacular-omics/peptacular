import json
import runpy
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

import peptacular as pt
from peptacular.proforma_json import _component_types


def components():
    formula = pt.ChargedFormula.from_composition({"C": 2, "H": 3})
    tag = pt.ModificationTags((pt.TagName("Oxidation"),))
    rule = pt.PositionRule(pt.Terminal.N_TERM)
    residue = pt.SequenceElement(pt.AminoAcid.C, (tag,))
    peptide = pt.Peptidoform((residue,))
    ion = pt.PeptidoformIon((peptide,), charge=2)
    return [
        pt.FormulaElement(pt.Element.C, 2, 13),
        formula,
        rule,
        pt.TagAccession("35", pt.CV.UNIMOD),
        pt.TagMass("+15.9"),
        pt.PositionScore("g1", 0.5),
        pt.TagName("Oxidation"),
        pt.TagInfo("note"),
        pt.TagCustom("custom"),
        pt.GlycanComponent(203, 2),
        pt.GlycanTag((pt.GlycanComponent(formula, 1),)),
        pt.PositionTag((rule,)),
        pt.LimitTag(1),
        pt.ComkpTag(),
        pt.ComupTag(),
        pt.IsotopeReplacement(pt.Element.C, 13),
        pt.GlobalChargeCarrier(formula, 1),
        tag,
        pt.ModificationAmbiguousPrimary("g1", tag, 0.7, (rule,), 1, True, False),
        pt.ModificationAmbiguousSecondary("g1", 0.3),
        pt.ModificationCrossLinker("XL1", tag),
        pt.FixedModification(tag, (rule,)),
        residue,
        pt.SequenceRegion((residue,), (tag,), True),
        peptide,
        ion,
        pt.CompoundPeptidoformIon((ion,)),
    ]


@pytest.mark.parametrize("value", components(), ids=lambda value: type(value).__name__)
def test_every_component_round_trips_and_validates(value):
    document = value.to_dict()
    Draft202012Validator(pt.get_proforma_json_schema()).validate(document)
    assert type(value).from_json(value.to_json()) == value


def test_schema_covers_registry_and_matches_generator():
    assert {type(value).__name__ for value in components()} == set(_component_types())
    generator = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts/generate_json_schema.py"))
    schema = pt.get_proforma_json_schema()
    Draft202012Validator.check_schema(schema)
    assert schema == generator["build_schema"]()
    schema.clear()
    assert pt.get_proforma_json_schema()


@pytest.mark.parametrize("sequence", ["PEPTIDE", "PEM[Oxidation]TIDE/-2", "<13C>PEPTIDE", "(>name)PEPTIDE/[Na:z+1,H:z+1]"])
def test_annotation_json_preserves_calculations(sequence):
    original = pt.parse(sequence)
    document = original.to_dict()
    Draft202012Validator(pt.get_proforma_json_schema()).validate(document)
    restored = pt.ProFormaAnnotation.from_json(original.to_json())
    assert restored.serialize() == original.serialize()
    assert restored.mass() == original.mass()
    assert restored.comp() == original.comp()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("sequence", 42),
        ("names", []),
        ("charge", True),
        ("charge", [3]),
        (
            "intervals",
            [
                {
                    "start": True,
                    "end": 2,
                    "ambiguous": False,
                    "modifications": None,
                }
            ],
        ),
    ],
)
def test_annotation_rejects_wrong_field_types(field, value):
    data = pt.parse("PEPTIDE").to_dict()
    data[field] = value
    assert not Draft202012Validator(pt.get_proforma_json_schema()).is_valid(data)
    with pytest.raises(ValueError):
        pt.from_proforma_dict(data)


@pytest.mark.parametrize(
    ("field", "value"), [("occurance", True), ("element", "C"), ("element", {"$enum": "AminoAcid", "value": "C"}), ("isotope", 1.2), ("surprise", 1)]
)
def test_component_rejects_wrong_field_types(field, value):
    data = pt.FormulaElement(pt.Element.C, 1).to_dict()
    data[field] = value
    assert not Draft202012Validator(pt.get_proforma_json_schema()).is_valid(data)
    with pytest.raises(ValueError):
        pt.from_proforma_dict(data)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_json_rejects_nonfinite_numbers_in_both_directions(value):
    tag = pt.TagName("Phospho", score=value)
    with pytest.raises(ValueError):
        tag.to_dict()
    data = pt.TagName("Phospho").to_dict()
    data["score"] = value
    with pytest.raises(ValueError):
        pt.from_proforma_dict(data)
    with pytest.raises(ValueError):
        pt.from_proforma_json(json.dumps(data))


def test_json_rejects_duplicate_keys_and_enum_roots():
    text = pt.parse("PEPTIDE").to_json()
    text = text.replace('"sequence": "PEPTIDE"', '"sequence": "PEPTIDE", "sequence": "OTHER"')
    with pytest.raises(ValueError, match="Duplicate JSON field"):
        pt.from_proforma_json(text)
    data = {"$schema": pt.PROFORMA_JSON_SCHEMA_ID, "schema_version": "1.0", "$enum": "Element", "value": "C"}
    with pytest.raises(ValueError, match="root"):
        pt.from_proforma_dict(data)


def test_json_rejects_duplicate_internal_positions_and_wrong_names():
    data = pt.parse("PEM[Oxidation]TIDE").to_dict()
    data["modifications"]["internal"] *= 2
    with pytest.raises(ValueError, match="Duplicate internal"):
        pt.from_proforma_dict(data)
    data = pt.parse("PEPTIDE").to_dict()
    data["names"]["ion"] = 42
    with pytest.raises(ValueError, match="names.ion"):
        pt.from_proforma_dict(data)


def test_json_orders_internal_positions_deterministically():
    first = pt.ProFormaAnnotation("PEPTIDE", internal_mods={1: "Oxidation", 3: "Oxidation"})
    second = pt.ProFormaAnnotation("PEPTIDE", internal_mods={3: "Oxidation", 1: "Oxidation"})
    assert first.to_json() == second.to_json()

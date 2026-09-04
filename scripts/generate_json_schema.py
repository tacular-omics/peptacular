"""Regenerate the bundled schema from the supported JSON field types."""

import json
from pathlib import Path
from types import UnionType
from typing import Any, get_args, get_origin

from peptacular.proforma_json import (
    PROFORMA_JSON_SCHEMA_ID,
    PROFORMA_JSON_SCHEMA_VERSION,
    _component_types,
    _enum_types,
    _field_types,
)


def object_schema(properties: dict[str, Any], *, closed: bool = True) -> dict[str, Any]:
    result = {"type": "object", "properties": properties, "required": list(properties)}
    if closed:
        result["additionalProperties"] = False
    return result


def nullable(schema: dict[str, Any]) -> dict[str, Any]:
    return {"anyOf": [schema, {"type": "null"}]}


def type_schema(expected: Any) -> dict[str, Any]:
    origin = get_origin(expected)
    args = get_args(expected)
    if origin is UnionType:
        return {"anyOf": [type_schema(option) for option in args]}
    if origin is tuple:
        return {"type": "array", "items": type_schema(args[0])}
    primitives = {str: "string", int: "integer", float: "number", bool: "boolean", type(None): "null"}
    if expected in primitives:
        return {"type": primitives[expected]}
    if expected in _enum_types().values():
        return object_schema({"$enum": {"const": expected.__name__}, "value": {"enum": [item.value for item in expected]}})
    if expected in _component_types().values():
        return {"allOf": [{"$ref": f"#/$defs/{expected.__name__}"}], "unevaluatedProperties": False}
    raise TypeError(f"Unsupported schema field type: {expected}")


def build_schema() -> dict[str, Any]:
    counts = nullable({"type": "object", "additionalProperties": {"type": "integer"}})
    internal_entry = object_schema({"position": {"type": "integer"}, "modifications": counts})
    interval = object_schema(
        {
            "start": {"type": "integer"},
            "end": {"type": "integer"},
            "ambiguous": {"type": "boolean"},
            "modifications": counts,
        }
    )
    modifications = dict.fromkeys(["isotope", "fixed", "labile", "unlocalized", "n_terminal", "c_terminal"], counts)
    modifications["internal"] = nullable({"type": "array", "items": internal_entry})
    definitions = {
        "ProFormaAnnotation": object_schema(
            {
                "$type": {"const": "ProFormaAnnotation"},
                "sequence": type_schema(str | None),
                "names": object_schema({name: type_schema(str | None) for name in ["compound", "ion", "peptide"]}),
                "modifications": object_schema(modifications),
                "intervals": nullable({"type": "array", "items": interval}),
                "charge": {"anyOf": [{"type": "null"}, {"type": "integer"}, {"type": "array", "items": {"type": "string"}}]},
            },
            closed=False,
        ),
    }
    for name, component_type in _component_types().items():
        properties = {"$type": {"const": name}}
        properties.update({field: type_schema(expected) for field, expected in _field_types(component_type).items()})
        definitions[name] = object_schema(properties, closed=False)
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": PROFORMA_JSON_SCHEMA_ID,
        "title": "Peptacular ProForma JSON 1.0",
        "description": "Typed representation of annotations and components. Scientific validity is checked separately by Peptacular.",
        "type": "object",
        "properties": {"$schema": {"const": PROFORMA_JSON_SCHEMA_ID}, "schema_version": {"const": PROFORMA_JSON_SCHEMA_VERSION}},
        "required": ["$schema", "schema_version"],
        "anyOf": [{"$ref": f"#/$defs/{name}"} for name in definitions],
        "unevaluatedProperties": False,
        "$defs": definitions,
    }


if __name__ == "__main__":
    path = Path(__file__).resolve().parent.parent / "src/peptacular/schemas/proforma-json-v1.schema.json"
    path.write_text(json.dumps(build_schema(), indent=2) + "\n", encoding="utf-8")

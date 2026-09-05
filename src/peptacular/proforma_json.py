"""Versioned, lossless JSON encoding for Peptacular ProForma objects.

The decoder deliberately uses a closed registry of known value types.  JSON input can
therefore never select an arbitrary Python module or class.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from enum import Enum
from functools import lru_cache
from importlib.resources import files as resource_files
from math import isfinite
from types import UnionType
from typing import Any, cast, get_args, get_origin, get_type_hints

PROFORMA_JSON_SCHEMA_VERSION = "1.0"
PROFORMA_JSON_SCHEMA_ID = "https://peptacular.readthedocs.io/en/latest/proforma-json-v1.schema.json"


@lru_cache(maxsize=32)
def _field_types(component_type: type[Any]) -> dict[str, Any]:
    hints = get_type_hints(component_type)
    return {field.name: hints[field.name] for field in fields(component_type)}


def _matches_type(value: Any, expected: Any) -> bool:
    origin = get_origin(expected)
    args = get_args(expected)
    if origin is UnionType:
        return any(_matches_type(value, option) for option in args)
    if origin is tuple:
        return isinstance(value, tuple) and all(_matches_type(item, args[0]) for item in value)
    if expected is float:
        return type(value) in (int, float) and (not isinstance(value, float) or isfinite(value))
    return type(value) is expected


def _check_type(value: Any, expected: Any, field_name: str) -> None:
    if not _matches_type(value, expected):
        raise ValueError(f"{field_name} must have type {expected}")


@lru_cache(maxsize=1)
def _component_types() -> dict[str, type[Any]]:
    from .proforma_components import comps

    names = (
        "FormulaElement",
        "ChargedFormula",
        "PositionRule",
        "TagAccession",
        "TagMass",
        "PositionScore",
        "TagName",
        "TagInfo",
        "TagCustom",
        "GlycanComponent",
        "GlycanTag",
        "PositionTag",
        "LimitTag",
        "ComkpTag",
        "ComupTag",
        "IsotopeReplacement",
        "GlobalChargeCarrier",
        "ModificationTags",
        "ModificationAmbiguousPrimary",
        "ModificationAmbiguousSecondary",
        "ModificationCrossLinker",
        "FixedModification",
        "SequenceElement",
        "SequenceRegion",
        "Peptidoform",
        "PeptidoformIon",
        "CompoundPeptidoformIon",
    )
    return {name: getattr(comps, name) for name in names}


@lru_cache(maxsize=1)
def _enum_types() -> dict[str, type[Enum]]:
    from tacular import AminoAcid, Element, Monosaccharide

    from .constants import CV, Terminal

    return {enum_type.__name__: enum_type for enum_type in (AminoAcid, Element, Monosaccharide, CV, Terminal)}


def _mod_counts(value: Mapping[str, int] | None) -> dict[str, int] | None:
    return dict(value) if value is not None else None


def _encode_annotation(annotation: Any) -> dict[str, Any]:
    return {
        "$type": "ProFormaAnnotation",
        "sequence": annotation._sequence,
        "names": {
            "compound": annotation._compound_name,
            "ion": annotation._ion_name,
            "peptide": annotation._peptide_name,
        },
        "modifications": {
            "isotope": _mod_counts(annotation._isotope_mods),
            "fixed": _mod_counts(annotation._static_mods),
            "labile": _mod_counts(annotation._labile_mods),
            "unlocalized": _mod_counts(annotation._unknown_mods),
            "n_terminal": _mod_counts(annotation._nterm_mods),
            "c_terminal": _mod_counts(annotation._cterm_mods),
            "internal": (
                [{"position": position, "modifications": dict(modifications)} for position, modifications in sorted(annotation._internal_mods.items())]
                if annotation._internal_mods is not None
                else None
            ),
        },
        "intervals": (
            [
                {
                    "start": interval.start,
                    "end": interval.end,
                    "ambiguous": interval.ambiguous,
                    "modifications": _mod_counts(interval._mods),
                }
                for interval in annotation._intervals
            ]
            if annotation._intervals is not None
            else None
        ),
        "charge": list(annotation._charge) if isinstance(annotation._charge, list) else annotation._charge,
    }


def _encode(value: Any) -> Any:
    from .annotation import ProFormaAnnotation

    if isinstance(value, ProFormaAnnotation):
        return _encode_annotation(value)
    if isinstance(value, Enum):
        if _enum_types().get(type(value).__name__) is not type(value):
            raise TypeError(f"Unsupported ProForma JSON enum: {type(value).__name__}")
        return {"$enum": type(value).__name__, "value": value.value}
    if is_dataclass(value) and _component_types().get(type(value).__name__) is type(value):
        for name, expected in _field_types(type(value)).items():
            _check_type(getattr(value, name), expected, f"{type(value).__name__}.{name}")
        encoded = {"$type": type(value).__name__}
        encoded.update((field.name, _encode(getattr(value, field.name))) for field in fields(value))
        return encoded
    if isinstance(value, tuple):
        return [_encode(item) for item in value]
    if isinstance(value, float) and not isfinite(value):
        raise ValueError("JSON numbers must be finite")
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"Unsupported ProForma JSON value: {type(value).__name__}")


def _require_keys(data: Mapping[str, Any], required: set[str], optional: set[str] | None = None) -> None:
    optional = optional or set()
    missing = required - data.keys()
    unknown = data.keys() - required - optional
    if missing:
        raise ValueError(f"Missing required JSON field(s): {', '.join(sorted(missing))}")
    if unknown:
        raise ValueError(f"Unknown JSON field(s): {', '.join(sorted(unknown))}")


def _decode_mod_counts(value: Any, field_name: str) -> dict[str, int] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ValueError(f"{field_name} must be an object or null")
    result: dict[str, int] = {}
    for modification, count in value.items():
        if not isinstance(modification, str) or not isinstance(count, int) or isinstance(count, bool):
            raise ValueError(f"{field_name} must map modification strings to integer counts")
        result[modification] = count
    return result


def _decode_annotation(data: Mapping[str, Any]) -> Any:
    from .annotation import Interval, ProFormaAnnotation

    _require_keys(data, {"$type", "sequence", "names", "modifications", "intervals", "charge"})
    names = data["names"]
    modifications = data["modifications"]
    if not isinstance(names, Mapping) or not isinstance(modifications, Mapping):
        raise ValueError("names and modifications must be JSON objects")
    _require_keys(names, {"compound", "ion", "peptide"})
    _require_keys(modifications, {"isotope", "fixed", "labile", "unlocalized", "n_terminal", "c_terminal", "internal"})
    _check_type(data["sequence"], str | None, "sequence")
    for name, value in names.items():
        _check_type(value, str | None, f"names.{name}")

    internal_value = modifications["internal"]
    internal: dict[int, dict[str, int]] | None = None
    if internal_value is not None:
        if not isinstance(internal_value, list):
            raise ValueError("modifications.internal must be an array or null")
        internal = {}
        for entry in internal_value:
            if not isinstance(entry, Mapping):
                raise ValueError("Each internal modification entry must be an object")
            _require_keys(entry, {"position", "modifications"})
            position = entry["position"]
            if not isinstance(position, int) or isinstance(position, bool):
                raise ValueError("Internal modification position must be an integer")
            if position in internal:
                raise ValueError(f"Duplicate internal modification position: {position}")
            internal[position] = _decode_mod_counts(entry["modifications"], "internal modifications") or {}

    intervals_value = data["intervals"]
    intervals: list[Interval] | None = None
    if intervals_value is not None:
        if not isinstance(intervals_value, list):
            raise ValueError("intervals must be an array or null")
        intervals = []
        for entry in intervals_value:
            if not isinstance(entry, Mapping):
                raise ValueError("Each interval must be an object")
            _require_keys(entry, {"start", "end", "ambiguous", "modifications"})
            _check_type(entry["start"], int, "interval.start")
            _check_type(entry["end"], int, "interval.end")
            _check_type(entry["ambiguous"], bool, "interval.ambiguous")
            intervals.append(
                Interval(
                    start=entry["start"],
                    end=entry["end"],
                    ambiguous=entry["ambiguous"],
                    mods=_decode_mod_counts(entry["modifications"], "interval modifications"),
                )
            )

    charge = data["charge"]
    if charge is not None and not isinstance(charge, (int, list)):
        raise ValueError("charge must be an integer, an array of adduct strings, or null")
    if isinstance(charge, bool) or (isinstance(charge, list) and not all(isinstance(item, str) for item in charge)):
        raise ValueError("charge must be an integer, an array of adduct strings, or null")

    return ProFormaAnnotation(
        sequence=data["sequence"],
        compound_name=names["compound"],
        ion_name=names["ion"],
        peptide_name=names["peptide"],
        isotope_mods=_decode_mod_counts(modifications["isotope"], "isotope modifications"),
        static_mods=_decode_mod_counts(modifications["fixed"], "fixed modifications"),
        labile_mods=_decode_mod_counts(modifications["labile"], "labile modifications"),
        unknown_mods=_decode_mod_counts(modifications["unlocalized"], "unlocalized modifications"),
        nterm_mods=_decode_mod_counts(modifications["n_terminal"], "N-terminal modifications"),
        cterm_mods=_decode_mod_counts(modifications["c_terminal"], "C-terminal modifications"),
        internal_mods=internal,
        intervals=intervals,
        charge=charge,
    )


def _decode(value: Any) -> Any:
    if isinstance(value, list):
        return tuple(_decode(item) for item in value)
    if not isinstance(value, Mapping):
        if isinstance(value, float) and not isfinite(value):
            raise ValueError("JSON numbers must be finite")
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        raise ValueError(f"Unsupported JSON value: {type(value).__name__}")
    if "$enum" in value:
        _require_keys(value, {"$enum", "value"})
        enum_name = value["$enum"]
        _check_type(enum_name, str, "$enum")
        enum_type = _enum_types().get(enum_name)
        if enum_type is None:
            raise ValueError(f"Unknown ProForma JSON enum type: {enum_name!r}")
        return enum_type(value["value"])
    object_name = value.get("$type")
    _check_type(object_name, str, "$type")
    if object_name == "ProFormaAnnotation":
        return _decode_annotation(value)
    component_type = _component_types().get(object_name)
    if component_type is None:
        raise ValueError(f"Unknown ProForma JSON object type: {object_name!r}")
    field_types = _field_types(component_type)
    field_names = set(field_types)
    _require_keys(value, field_names | {"$type"})
    decoded = {name: _decode(value[name]) for name in field_names}
    for name, expected in field_types.items():
        _check_type(decoded[name], expected, f"{object_name}.{name}")
    return component_type(**decoded)


def to_proforma_dict(value: Any) -> dict[str, Any]:
    """Encode a supported ProForma object as a versioned JSON-compatible mapping."""
    encoded = _encode(value)
    if not isinstance(encoded, dict) or "$type" not in encoded:
        raise TypeError("The root ProForma JSON value must be a supported object")
    return {
        "$schema": PROFORMA_JSON_SCHEMA_ID,
        "schema_version": PROFORMA_JSON_SCHEMA_VERSION,
        **encoded,
    }


def from_proforma_dict(data: Mapping[str, Any], expected_type: type[Any] | None = None) -> Any:
    """Decode a versioned ProForma mapping, optionally enforcing its root type."""
    if not isinstance(data, Mapping):
        raise TypeError("ProForma JSON data must be a mapping")
    if "$schema" not in data:
        raise ValueError("Missing required JSON field: $schema")
    version = data.get("schema_version")
    if version != PROFORMA_JSON_SCHEMA_VERSION:
        raise ValueError(f"Unsupported ProForma JSON schema version: {version!r}")
    payload = dict(data)
    payload.pop("schema_version")
    schema_id = payload.pop("$schema", PROFORMA_JSON_SCHEMA_ID)
    if schema_id != PROFORMA_JSON_SCHEMA_ID:
        raise ValueError(f"Unsupported ProForma JSON schema: {schema_id!r}")
    if "$type" not in payload:
        raise ValueError("The root ProForma JSON value must be a supported object with $type")
    decoded = _decode(payload)
    if expected_type is not None and not isinstance(decoded, expected_type):
        raise TypeError(f"Expected {expected_type.__name__}, got {type(decoded).__name__}")
    return decoded


def to_proforma_json(value: Any, *, indent: int | None = None) -> str:
    """Encode a supported ProForma object as deterministic JSON text."""
    return json.dumps(to_proforma_dict(value), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True)


def from_proforma_json(data: str | bytes | bytearray, expected_type: type[Any] | None = None) -> Any:
    """Decode JSON text produced by :func:`to_proforma_json`."""
    parsed = json.loads(data, object_pairs_hook=_unique_object, parse_constant=_invalid_constant)
    if not isinstance(parsed, Mapping):
        raise ValueError("The root ProForma JSON value must be an object")
    return from_proforma_dict(parsed, expected_type=expected_type)


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON field: {key}")
        result[key] = value
    return result


def _invalid_constant(value: str) -> Any:
    raise ValueError(f"JSON numbers must be finite: {value}")


def get_proforma_json_schema() -> dict[str, Any]:
    """Return a fresh copy of the bundled JSON Schema for the stable representation."""
    schema_path = resource_files("peptacular").joinpath("schemas/proforma-json-v1.schema.json")
    return cast(dict[str, Any], json.loads(schema_path.read_text(encoding="utf-8")))


__all__ = [
    "PROFORMA_JSON_SCHEMA_ID",
    "PROFORMA_JSON_SCHEMA_VERSION",
    "from_proforma_dict",
    "from_proforma_json",
    "get_proforma_json_schema",
    "to_proforma_dict",
    "to_proforma_json",
]

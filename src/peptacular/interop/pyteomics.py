"""Conversions between Peptacular and Pyteomics objects."""

import re
from collections.abc import Mapping
from typing import Any, cast

from peptacular.annotation import ProFormaAnnotation
from peptacular.proforma_components import ChargedFormula

from ._errors import InteropConversionError
from ._optional import require_dependency

_PEPTACULAR_ISOTOPE = re.compile(r"^(\d+)([A-Z][a-z]?)$")
_PYTEOMICS_ISOTOPE = re.compile(r"^([A-Z][a-z]?)\[(\d+)\]$")


def to_pyteomics(annotation: ProFormaAnnotation) -> Any:
    """Convert an annotation to :class:`pyteomics.proforma.ProForma`.

    :param annotation: Peptacular annotation to convert.
    :type annotation: ProFormaAnnotation
    :return: Parsed Pyteomics ProForma object.
    :rtype: Any
    :raises MissingOptionalDependencyError: If Pyteomics is not installed.
    :raises InteropConversionError: If Pyteomics cannot represent the annotation.
    """
    proforma = require_dependency("pyteomics.proforma", "pyteomics")
    try:
        return proforma.ProForma.parse(annotation.serialize())
    except Exception as exc:
        raise InteropConversionError(f"Pyteomics could not parse the annotation {annotation.serialize()!r}: {exc}") from exc


def from_pyteomics(value: Any) -> ProFormaAnnotation:
    """Convert a Pyteomics ProForma object to a Peptacular annotation.

    :param value: A :class:`pyteomics.proforma.ProForma` object.
    :type value: Any
    :return: Peptacular annotation.
    :rtype: ProFormaAnnotation
    :raises MissingOptionalDependencyError: If Pyteomics is not installed.
    :raises TypeError: If *value* is not a Pyteomics ProForma object.
    :raises InteropConversionError: If Peptacular cannot parse the serialized object.
    """
    proforma = require_dependency("pyteomics.proforma", "pyteomics")
    if not isinstance(value, proforma.ProForma):
        raise TypeError(f"Expected pyteomics.proforma.ProForma, got {type(value).__name__}")
    try:
        return ProFormaAnnotation.parse(str(value))
    except Exception as exc:
        raise InteropConversionError(f"Peptacular could not parse the Pyteomics value {value!s}: {exc}") from exc


def to_pyteomics_composition(composition: Mapping[str, int] | ChargedFormula) -> Any:
    """Convert elemental counts to :class:`pyteomics.mass.Composition`.

    Isotopes are translated from Peptacular keys such as ``13C`` to Pyteomics
    keys such as ``C[13]``. Formula charge metadata is not part of a Pyteomics
    ``Composition`` and is therefore rejected.

    :param composition: Element-count mapping or neutral charged-formula value.
    :type composition: Mapping[str, int] | ChargedFormula
    :return: Pyteomics composition.
    :rtype: Any
    :raises InteropConversionError: If a charged formula is supplied.
    """
    mass = require_dependency("pyteomics.mass", "pyteomics")
    if isinstance(composition, ChargedFormula):
        if composition.charge not in (None, 0):
            raise InteropConversionError("Pyteomics Composition does not preserve ChargedFormula charge metadata")
        counts = composition.get_dict_composition()
    else:
        counts = dict(composition)

    converted: dict[str, int] = {}
    for key, count in counts.items():
        isotope = _PEPTACULAR_ISOTOPE.fullmatch(key)
        target_key = f"{isotope.group(2)}[{isotope.group(1)}]" if isotope else key
        converted[target_key] = converted.get(target_key, 0) + count
    try:
        return mass.Composition(converted)
    except Exception as exc:
        raise InteropConversionError(f"Pyteomics could not construct a composition from {converted!r}: {exc}") from exc


def from_pyteomics_composition(composition: Mapping[str, int]) -> ChargedFormula:
    """Convert a Pyteomics composition-like mapping to ``ChargedFormula``.

    :param composition: Pyteomics composition or compatible mapping.
    :type composition: Mapping[str, int]
    :return: Neutral Peptacular formula.
    :rtype: ChargedFormula
    :raises InteropConversionError: If a key is not an element or isotope.
    """
    require_dependency("pyteomics.mass", "pyteomics")
    converted: dict[str, int] = {}
    for key, count in composition.items():
        isotope = _PYTEOMICS_ISOTOPE.fullmatch(key)
        target_key = f"{isotope.group(2)}{isotope.group(1)}" if isotope else key
        if not re.fullmatch(r"(?:\d+)?[A-Z][a-z]?", target_key):
            raise InteropConversionError(f"Pyteomics composition key {key!r} is not representable as a Peptacular element")
        converted[target_key] = converted.get(target_key, 0) + count
    try:
        return ChargedFormula.from_composition(cast(Any, converted))
    except Exception as exc:
        raise InteropConversionError(f"Peptacular could not construct a formula from {converted!r}: {exc}") from exc

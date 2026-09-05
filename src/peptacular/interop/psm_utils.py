"""Conversions between Peptacular and psm_utils peptidoforms."""

from typing import Any

from peptacular.annotation import ProFormaAnnotation

from ._errors import InteropConversionError
from ._optional import require_dependency
from ._roundtrip import check_round_trip


def to_psm_utils(annotation: ProFormaAnnotation) -> Any:
    """Convert an annotation to :class:`psm_utils.Peptidoform`.

    :param annotation: Peptacular annotation to convert.
    :type annotation: ProFormaAnnotation
    :return: psm_utils peptidoform.
    :rtype: Any
    :raises MissingOptionalDependencyError: If psm_utils is not installed.
    :raises InteropConversionError: If psm_utils cannot represent the annotation.
    """
    psm_utils = require_dependency("psm_utils", "psm-utils")
    try:
        converted = psm_utils.Peptidoform(annotation.serialize())
        check_round_trip(annotation, converted.proforma, "psm_utils")
        return converted
    except Exception as exc:
        raise InteropConversionError(f"psm_utils could not parse the annotation {annotation.serialize()!r}: {exc}") from exc


def from_psm_utils(value: Any) -> ProFormaAnnotation:
    """Convert a :class:`psm_utils.Peptidoform` to Peptacular.

    PSM metadata is intentionally outside the scope of this conversion.

    :param value: psm_utils peptidoform.
    :type value: Any
    :return: Peptacular annotation.
    :rtype: ProFormaAnnotation
    :raises MissingOptionalDependencyError: If psm_utils is not installed.
    :raises TypeError: If *value* is not a psm_utils Peptidoform.
    :raises InteropConversionError: If Peptacular cannot parse the ProForma text.
    """
    psm_utils = require_dependency("psm_utils", "psm-utils")
    if not isinstance(value, psm_utils.Peptidoform):
        raise TypeError(f"Expected psm_utils.Peptidoform, got {type(value).__name__}")
    try:
        return ProFormaAnnotation.parse(value.proforma)
    except Exception as exc:
        raise InteropConversionError(f"Peptacular could not parse the psm_utils value {value.proforma!r}: {exc}") from exc

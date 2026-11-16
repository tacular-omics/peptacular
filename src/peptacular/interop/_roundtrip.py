"""Guard text-based adapters against silent changes in target parsers."""

from peptacular.annotation import ProFormaAnnotation

from ..diagnostics import ProFormaFormatError
from ._errors import InteropConversionError


def check_round_trip(annotation: ProFormaAnnotation, serialized: str, package: str) -> None:
    original = ProFormaAnnotation.parse(annotation.serialize())
    try:
        restored = ProFormaAnnotation.parse(serialized)
    except ProFormaFormatError as exc:
        # The target library wrote something peptacular cannot read back, so the round trip failed.
        raise InteropConversionError(f"{package} does not preserve this annotation during a ProForma round trip: {exc}") from exc
    if original.to_dict() != restored.to_dict():
        raise InteropConversionError(f"{package} does not preserve this annotation during a ProForma round trip")

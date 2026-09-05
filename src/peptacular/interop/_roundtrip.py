"""Guard text-based adapters against silent changes in target parsers."""

from peptacular.annotation import ProFormaAnnotation

from ._errors import InteropConversionError


def check_round_trip(annotation: ProFormaAnnotation, serialized: str, package: str) -> None:
    original = ProFormaAnnotation.parse(annotation.serialize())
    restored = ProFormaAnnotation.parse(serialized)
    if original.to_dict() != restored.to_dict():
        raise InteropConversionError(f"{package} does not preserve this annotation during a ProForma round trip")

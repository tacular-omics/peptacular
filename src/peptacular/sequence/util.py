from typing import Protocol, runtime_checkable

from ..annotation import ProFormaAnnotation

__all__ = [
    "HasSequence",
    "as_sequence_input",
    "sequence_to_annotation",
    "round_to_precision",
    "get_annotation_input",
    "is_sequence_valid",
]


@runtime_checkable
class HasSequence(Protocol):
    """Any object with a ``sequence`` string attribute, such as a FASTA or PEFF entry.

    Sequence functions accept these directly and read ``obj.sequence`` as ProForma, so
    ``pt.digest(entry, "trypsin")`` works for a ``fastatacular.SequenceEntry`` without
    peptacular depending on fastatacular.
    """

    @property
    def sequence(self) -> str: ...


def as_sequence_input(sequence: object) -> str | ProFormaAnnotation | None:
    """Return a str or annotation for a supported single input, or ``None`` if it is not one."""
    if isinstance(sequence, (str, ProFormaAnnotation)):
        return sequence
    value = getattr(sequence, "sequence", None)
    return value if isinstance(value, str) else None


def sequence_to_annotation(sequence: str) -> ProFormaAnnotation:
    return ProFormaAnnotation.parse(sequence)


def round_to_precision(value: float, precision: int | None = None) -> float:
    if precision is not None:
        value = round(value, precision)
    return value


def get_annotation_input(
    sequence: str | ProFormaAnnotation | HasSequence,
    copy: bool = True,
) -> ProFormaAnnotation:
    if isinstance(sequence, ProFormaAnnotation):
        return sequence.copy() if copy else sequence
    value = as_sequence_input(sequence)
    if isinstance(value, str):
        return sequence_to_annotation(value)
    raise TypeError(
        "Input sequence must be a ProForma str, a ProFormaAnnotation, or an object with a str "
        f"'sequence' attribute (e.g. a FASTA entry), got {type(sequence).__name__}: {sequence!r}"
    )


def is_sequence_valid(sequence: str | ProFormaAnnotation | HasSequence) -> bool:
    """
    Checks if the input sequence is a valid ProForma sequence.

    :param sequence: The sequence or ProFormaAnnotation object to be validated.
    :type sequence: Union[str, ProFormaAnnotation]

    :return: True if the sequence is a valid ProForma sequence, False otherwise.
    :rtype: bool

    """

    value = as_sequence_input(sequence)
    if isinstance(value, str):
        try:
            _ = sequence_to_annotation(value)
        except ValueError:
            return False
    return True

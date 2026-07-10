from ..annotation import ProFormaAnnotation


def sequence_to_annotation(sequence: str) -> ProFormaAnnotation:
    return ProFormaAnnotation.parse_single(sequence)


def round_to_precision(value: float, precision: int | None = None) -> float:
    if precision is not None:
        value = round(value, precision)
    return value


def get_annotation_input(
    sequence: str | ProFormaAnnotation,
    copy: bool = True,
) -> ProFormaAnnotation:
    if isinstance(sequence, str):
        return sequence_to_annotation(sequence)
    elif isinstance(sequence, ProFormaAnnotation):
        return sequence.copy() if copy else sequence
    raise TypeError(f"Input sequence must be a ProForma str or a ProFormaAnnotation, got {type(sequence).__name__}: {sequence!r}")


def is_sequence_valid(sequence: str | ProFormaAnnotation) -> bool:
    """
    Checks if the input sequence is a valid ProForma sequence.

    :param sequence: The sequence or ProFormaAnnotation object to be validated.
    :type sequence: Union[str, ProFormaAnnotation]

    :return: True if the sequence is a valid ProForma sequence, False otherwise.
    :rtype: bool

    """

    if isinstance(sequence, str):
        try:
            _ = sequence_to_annotation(sequence)
        except ValueError:
            return False
    return True

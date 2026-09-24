from tacular import FRAGMENT_ION_LOOKUP, FragmentIonInfo, IonType, IonTypeLiteral

from ..diagnostics import InvalidPositionError, UnsupportedOperationError

__all__ = [
    "to_ion_type",
    "validate_position",
]


def to_ion_type(ion_type: IonType | IonTypeLiteral | str) -> IonType:
    """Return ``ion_type`` as an :class:`IonType`, raising a clear error for unknown values.

    :raises UnsupportedOperationError: ``ion_type`` is not a known ion type.
    """
    if isinstance(ion_type, IonType):
        return ion_type
    try:
        return IonType(ion_type)
    except ValueError:
        valid = ", ".join(repr(t.value) for t in IonType)
        raise UnsupportedOperationError(f"Unknown ion type {ion_type!r}. Valid ion types: {valid}") from None


def _validate_position(ion_type: IonType, position: int | tuple[int, int] | None, sequence_length: int) -> tuple[int, int] | None:
    """Validate fragment position based on ion type."""

    ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]

    if position is None:
        return None  # whole sequence, no position-specific adjustments

    if isinstance(position, int):
        if ion_info.is_forward:
            return (0, position)
        if ion_info.is_backward:
            return (sequence_length - position, sequence_length)
        if ion_info.ion_type == IonType.IMMONIUM:
            return (position - 1, position)
        raise InvalidPositionError(f"Integer position {position} not valid for ion type {ion_type.name}, which is not a forward or backward ion.")

    if isinstance(position, tuple):
        if ion_info.is_internal:
            return (position[0] - 1, position[1])
        raise InvalidPositionError(f"Tuple position {position} not valid for ion type {ion_type.name}, which is not an internal ion.")

    raise InvalidPositionError(f"Invalid position type: {type(position)} for ion type {ion_type.name}.")


def validate_position(ion_type: IonType | IonTypeLiteral, position: int | tuple[int, int] | None, sequence_length: int) -> tuple[int, int] | None:
    """
    Validates the passed position and return the true zero index position tuple for the fragment.
    """
    pos: tuple[int, int] | None = _validate_position(to_ion_type(ion_type), position, sequence_length)
    if pos is None:
        return None
    start, end = pos
    if start < 0 or end > sequence_length or start >= end:
        raise InvalidPositionError(
            f"Invalid fragment position {position!r} for {IonType(ion_type).name} on a sequence of length {sequence_length} (resolved to residues {pos})."
        )
    return pos

from __future__ import annotations

import math
import statistics
from collections.abc import Generator, Iterable, Mapping, Sequence

from .data import (
    NEGATIVE_AMINO_ACIDS,
    POSITIVE_AMINO_ACIDS,
    PROPERTY_SCALES,
    ChargeScale,
    SecondaryStructureMethod,
    secondary_structure_scales_by_name,
)
from .types import (
    AggregationMethod,
    AggregationMethodLiteral,
    MissingAAHandling,
    MissingAAHandlingLiteral,
    WeightingMethods,
    WeightingMethodsLiteral,
)
from .weights import get_weights

# Handle ambiguous amino acids
AMIGUOUS_AMINO_ACID_MAP: dict[str, tuple[str, ...]] = {
    "B": ("D", "N"),  # Aspartic acid or Asparagine
    "J": ("L", "I"),  # Leucine or Isoleucine
    "Z": ("E", "Q"),  # Glutamic acid or Glutamine
}


def _calculate_median(values: Sequence[float]) -> float:
    """Calculate median of a list of values"""
    return statistics.median(values)


def _get_default_value(handling: MissingAAHandling, aa_data: Mapping[str, float], aa: str) -> float | None:
    """Calculate default value based on missing AA handling strategy"""

    match handling:
        case MissingAAHandling.ZERO:
            return 0.0
        case MissingAAHandling.AVG:
            return sum(aa_data.values()) / len(aa_data)
        case MissingAAHandling.MIN:
            return min(aa_data.values())
        case MissingAAHandling.MAX:
            return max(aa_data.values())
        case MissingAAHandling.MEDIAN:
            return _calculate_median(list(aa_data.values()))
        case MissingAAHandling.ERROR:
            raise ValueError(f"Invalid amino acid: {aa}")
        case MissingAAHandling.SKIP:
            return None
        case _:
            valid_options = ", ".join([f"'{opt.value}'" for opt in MissingAAHandling])
            raise ValueError(f"Invalid missing_aa_handling: {handling}. Choose from {valid_options}")


def _normalize_value(value: float, aa_data: Mapping[str, float]) -> float:
    """Normalize value to 0-1 range based on min/max in aa_data"""
    min_value = min(aa_data.values())
    max_value = max(aa_data.values())
    if max_value == min_value:
        return 0.0  # Avoid division by zero
    return (value - min_value) / (max_value - min_value)


def _apply_weighting_and_normalization(value: float, aa_data: Mapping[str, float], weighting_scheme: float, normalize: bool) -> float:
    """Apply normalization and weighting to a value"""
    if normalize:
        value = _normalize_value(value, aa_data)
    return value * weighting_scheme


def _get_ambiguous_aa_value(
    aa: str,
    constituent_aas: Sequence[str],
    aa_data: Mapping[str, float],
    missing_aa_handling: MissingAAHandling,
    weighting_scheme: float,
    normalize: bool,
) -> float:
    """Get average value for ambiguous amino acids (B, J, Z)"""
    values: list[float] = []
    for constituent_aa in constituent_aas:
        val = _get_aa_value(constituent_aa, aa_data, missing_aa_handling, 1.0, False)
        if val is not None:
            values.append(val)

    if not values:
        raise ValueError(f"No valid values found for ambiguous amino acid {aa}")

    avg_value = sum(values) / len(values)
    return _apply_weighting_and_normalization(avg_value, aa_data, weighting_scheme, normalize)


def _get_aa_value(
    aa: str,
    aa_data: Mapping[str, float],
    missing_aa_handling: str,
    weighting_scheme: float = 1.0,
    normalize: bool = False,
) -> float | None:
    """
    Get amino acid value with configurable handling for missing values.

    Args:
        aa: Single letter amino acid code
        aa_data: Dictionary mapping amino acid codes to values
        missing_aa_handling: Strategy for handling missing amino acids
        weighting_scheme: Multiplier applied to final value
        normalize: Whether to normalize value to 0-1 range

    Returns:
        Amino acid value or None if skipped

    Raises:
        ValueError: If amino acid is invalid and handling is ERROR
    """

    missing_aa_handling = MissingAAHandling.from_str(missing_aa_handling)

    # Handle standard amino acids
    if aa in aa_data:
        value = aa_data[aa]
        return _apply_weighting_and_normalization(value, aa_data, weighting_scheme, normalize)

    if aa in AMIGUOUS_AMINO_ACID_MAP:
        return _get_ambiguous_aa_value(
            aa,
            AMIGUOUS_AMINO_ACID_MAP[aa],
            aa_data,
            missing_aa_handling,
            weighting_scheme,
            normalize,
        )

    # Handle X (any amino acid) - use average of all
    if aa == "X":
        avg_value = sum(aa_data.values()) / len(aa_data)
        return _apply_weighting_and_normalization(avg_value, aa_data, weighting_scheme, normalize)

    # Handle unknown amino acids
    default_value = _get_default_value(missing_aa_handling, aa_data, aa)
    if default_value is None:
        return None

    return _apply_weighting_and_normalization(default_value, aa_data, weighting_scheme, normalize)


def _generate_string_sliding_windows(
    sequence: str,
    window_size: int,
    reverse: bool = False,
) -> Generator[str]:
    """
    Generate sliding windows of a plain amino acid sequence.

    Args:
        sequence: The amino acid sequence string
        window_size: Size of each window
        reverse: If True, generate windows from right to left

    Yields:
        String subsequences representing each window

    Raises:
        ValueError: If window_size is invalid or sequence is empty
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")

    if window_size <= 0:
        raise ValueError("Window size must be positive")

    seq_len = len(sequence)
    if window_size > seq_len:
        raise ValueError(f"Window size {window_size} cannot be greater than sequence length {seq_len}.")

    if reverse:
        # Generate windows from right to left
        for start in range(seq_len - window_size, -1, -1):
            stop = start + window_size
            yield sequence[start:stop]
    else:
        # Generate windows from left to right
        for start in range(seq_len - window_size + 1):
            stop = start + window_size
            yield sequence[start:stop]


def calc_property(
    sequence: str,
    scale: str | Mapping[str, float],
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.ERROR,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods | Sequence[float]) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> float:
    """
    Generic function to calculate the average property value of a protein sequence.
    """
    weighting_scheme = get_weights(
        length=len(sequence),
        weights=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )

    if isinstance(scale, str):
        if scale not in PROPERTY_SCALES:
            raise ValueError(f"Scale '{scale}' not found in available property scales: {list(PROPERTY_SCALES.keys())}")
        aa_data = PROPERTY_SCALES[scale]
    else:
        aa_data = scale

    values: list[float] = []
    for i, aa in enumerate(sequence):
        val = _get_aa_value(
            aa=aa,
            aa_data=aa_data,
            missing_aa_handling=missing_aa_handling,
            weighting_scheme=weighting_scheme[i],
            normalize=normalize,
        )
        if isinstance(val, (int, float)):
            values.append(val)

    if aggregation_method == AggregationMethod.SUM:
        result = sum(values) if values else 0.0
    elif aggregation_method == AggregationMethod.AVG:
        result = sum(values) / len(values) if values else 0.0
    else:
        raise ValueError(f"Invalid aggregation method: {aggregation_method}. Choose '{AggregationMethod.SUM}' or '{AggregationMethod.AVG}'.")

    return result


def calc_window_property(
    sequence: str,
    scale: str | Mapping[str, float],
    window_size: int = 9,
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.ERROR,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods | Sequence[float]) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> list[float]:
    """Calculate property values over sliding windows"""
    weighting_scheme = get_weights(
        window_size,
        weights=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )

    results: list[float] = []
    for window_sequence in _generate_string_sliding_windows(sequence, window_size):
        if len(window_sequence) != window_size:
            raise ValueError(f"Window size {window_size} does not match sequence length {len(window_sequence)}.")

        # Calculate the property average for the current window
        window_value = calc_property(
            sequence=window_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )
        results.append(window_value)

    return results


def _count_residues(sequence: str) -> dict[str, int]:
    """Counts the occurrences of each amino acid residue in the sequence."""
    counts: dict[str, int] = {}
    for aa in sequence:
        counts[aa] = counts.get(aa, 0) + 1
    return counts


def _percent_residues(sequence: str) -> dict[str, float]:
    """Calculates the percentage of each amino acid residue in the sequence."""
    total = len(sequence)
    if total == 0:
        return {}

    counts = _count_residues(sequence)
    return {aa: count / total for aa, count in counts.items()}


def aa_property_percentage(
    sequence: str,
    residues: Iterable[str],
) -> float:
    """Calculates the percentage of specified amino acids in the sequence."""
    residue_perc = _percent_residues(sequence)
    val = sum(residue_perc.get(aa, 0) for aa in residues)
    return float(val)


def charge_at_ph(sequence: str, pH: float = 7.0) -> float:
    """Calculate net charge at given pH"""
    # Count amino acids
    aa_counts = _count_residues(sequence=sequence)

    # Get terminal residues
    if not sequence:
        return 0.0

    nterm, cterm = sequence[0], sequence[-1]

    # Calculate positive charge (basic groups)
    positive_charge = 0.0

    # N-terminal charge
    nterm_pK = _get_aa_value(
        aa=nterm,
        aa_data=PROPERTY_SCALES[ChargeScale.PK_NTERMINAL],
        missing_aa_handling=MissingAAHandling.ERROR,
    )

    if not isinstance(nterm_pK, (int, float)):
        raise ValueError(f"Invalid pK value for N-terminal amino acid '{nterm}': {nterm_pK}")

    partial_charge: float = 1.0 / (10 ** (pH - nterm_pK) + 1.0)
    positive_charge += partial_charge

    # Side chain positive charges
    for aa in POSITIVE_AMINO_ACIDS:
        count = float(aa_counts.get(aa, 0))
        if count > 0:
            pK = _get_aa_value(
                aa=aa,
                aa_data=PROPERTY_SCALES[ChargeScale.PK_SIDECHAIN],
                missing_aa_handling=MissingAAHandling.ERROR,
            )

            if not isinstance(pK, (int, float)) or pK <= 0:
                raise ValueError(f"Invalid pK value for side chain amino acid '{aa}': {pK}")

            partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
            positive_charge += count * partial_charge

    # Calculate negative charge (acidic groups)
    negative_charge = 0.0

    # C-terminal charge
    cterm_pK = _get_aa_value(
        aa=cterm,
        aa_data=PROPERTY_SCALES[ChargeScale.PK_CTERMINAL],
        missing_aa_handling=MissingAAHandling.ERROR,
    )

    if not isinstance(cterm_pK, (int, float)):
        raise ValueError(f"Invalid pK value for C-terminal amino acid '{cterm}': {cterm_pK}")

    partial_charge = 1.0 / (10 ** (cterm_pK - pH) + 1.0)
    negative_charge += partial_charge

    # Side chain negative charges
    for aa in NEGATIVE_AMINO_ACIDS:
        count = float(aa_counts.get(aa, 0))
        if count > 0:
            pK = _get_aa_value(
                aa=aa,
                aa_data=PROPERTY_SCALES[ChargeScale.PK_SIDECHAIN],
                missing_aa_handling=MissingAAHandling.ERROR,
            )

            if not isinstance(pK, (int, float)) or pK <= 0:
                raise ValueError(f"Invalid pK value for side chain amino acid '{aa}': {pK}")

            if pK > 0:  # Only calculate if pK exists (non-zero)
                partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
                negative_charge += count * partial_charge

    net_charge = positive_charge - negative_charge
    return net_charge


def secondary_structure(
    sequence: str,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
) -> dict[str, float]:
    """Calculate secondary structure propensities"""

    scale = SecondaryStructureMethod(scale)

    d: dict[str, float] = {}
    for structure_scale_name, structure_scale in secondary_structure_scales_by_name[scale].items():
        val = calc_property(
            sequence=sequence,
            scale=structure_scale,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        d[structure_scale_name] = val

    # Normalize the values to sum to 1
    total = sum(d.values())
    if total > 0:
        for key in d:
            d[key] /= total

    return d


def generate_partitions(
    sequence: str,
    scale: str | dict[str, float],
    num_windows: int = 5,
    aa_overlap: int = 0,
    missing_aa_handling: (MissingAAHandlingLiteral | MissingAAHandling) = MissingAAHandling.AVG,
    aggregation_method: (AggregationMethodLiteral | AggregationMethod) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (WeightingMethodsLiteral | WeightingMethods | Sequence[float]) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> list[float]:
    """
    Generate property values for sliding windows across the sequence.

    If the sequence cannot be evenly divided, windows are distributed as evenly as
    possible, with the first and last windows equally shortened.

    Args:
        sequence: Peptide sequence (string or ProFormaAnnotation)
        scale: Property scale name or dictionary mapping AA to values
        num_windows: Number of windows to divide the sequence into
        aa_overlap: Number of amino acids to overlap between adjacent windows
        missing_aa_handling: Strategy for handling missing amino acids
        aggregation_method: How to aggregate values within each window
        normalize: Whether to normalize values to 0-1 range
        weighting_scheme: How to weight positions within each window
        min_weight: Minimum weight for position-based weighting
        max_weight: Maximum weight for position-based weighting

    Returns:
        List of length num_windows containing property value per window

    Raises:
        ValueError: If parameters result in invalid window configuration
    """
    seq_len = len(sequence)

    if num_windows <= 0:
        raise ValueError("num_windows must be positive")

    if aa_overlap < 0:
        raise ValueError("aa_overlap cannot be negative")

    if seq_len == 0:
        raise ValueError("Sequence cannot be empty")

    if num_windows > seq_len and aa_overlap == 0:
        raise ValueError(f"Cannot create {num_windows} non-overlapping windows from sequence of length {seq_len}")

    if num_windows == 1:
        # Single window covers entire sequence
        value = calc_property(
            sequence=sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )
        return [value]

    # Handle case where sequence is very short
    if seq_len < num_windows:
        # For very short sequences, we need to allow windows to heavily overlap
        # or even repeat the same positions
        window_values: list[float] = []

        # Create windows that maximally overlap
        for i in range(num_windows):
            # Distribute windows as evenly as possible across the sequence
            # Each window will be at least 1 AA long
            if seq_len == 1:
                # Special case: single AA repeated for all windows
                start = 0
                end = 1
            else:
                # Interpolate position within sequence
                progress = i / max(1, (num_windows - 1))
                center = progress * (seq_len - 1)

                # Create a window around this center point
                half_window = max(0, aa_overlap) / 2.0
                start = max(0, int(center - half_window))
                end = min(seq_len, int(center + half_window + 1))

                # Ensure we have at least 1 character
                if end <= start:
                    end = start + 1
                if end > seq_len:
                    start = seq_len - 1
                    end = seq_len

            window_seq = sequence[start:end]

            value = calc_property(
                sequence=window_seq,
                scale=scale,
                missing_aa_handling=missing_aa_handling,
                aggregation_method=aggregation_method,
                normalize=normalize,
                weighting_scheme=weighting_scheme,
                min_weight=min_weight,
                max_weight=max_weight,
            )
            window_values.append(value)

        return window_values

    # Calculate ideal window size and step size for normal cases
    window_size = (seq_len + (num_windows - 1) * aa_overlap) / num_windows
    window_size = max(1, math.ceil(window_size))  # Round up, minimum 1
    step_size = max(1, window_size - aa_overlap)  # Ensure step_size is at least 1

    # If step_size is too small, we need to adjust
    if step_size < 1:
        step_size = 1
        aa_overlap = window_size - 1

    # For very small sequences, ensure window_size doesn't exceed sequence length
    if window_size > seq_len:
        window_size = seq_len
        step_size = max(1, (seq_len - aa_overlap) // max(1, (num_windows - 1)))
        if step_size < 1:
            step_size = 1

    # Calculate positions for each window
    window_positions: list[tuple[int, int]] = []

    if num_windows == 2:
        # Special case for 2 windows - place at start and end
        window_positions = [
            (0, min(window_size, seq_len)),
            (max(0, seq_len - window_size), seq_len),
        ]
    else:
        # Calculate total span if we used the ideal step size
        total_span = (num_windows - 1) * step_size + window_size

        if total_span <= seq_len:
            # Windows fit within sequence - center them
            offset = (seq_len - total_span) / 2.0
            for i in range(num_windows):
                start = int(i * step_size + offset)
                end = min(seq_len, start + window_size)
                window_positions.append((start, end))
        else:
            # Windows don't fit using the requested aa_overlap for every gap (the ideal
            # step_size would run past the end of the sequence). Use the exact requested
            # step_size for every window boundary, so aa_overlap is honored everywhere it
            # can be; any window whose ideal start would run past the end (not just the
            # last one — with enough windows requested, interior windows can overrun too)
            # is pulled left just enough to end exactly at the sequence's end, mirroring
            # the num_windows == 2 case above. window_size is already clamped to <= seq_len
            # at this point, so max(0, seq_len - window_size) is always a valid start.
            # A previous version only clamped the last window, which left interior windows
            # free to start past the sequence end and produce a phantom empty window.
            for i in range(num_windows):
                start = min(i * step_size, max(0, seq_len - window_size))
                end = min(seq_len, start + window_size)

                window_positions.append((start, end))

    # Process each window
    window_values: list[float] = []
    for _, (start, end) in enumerate(window_positions):
        window_seq = sequence[start:end]

        value = calc_property(
            sequence=window_seq,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )
        window_values.append(value)

    return window_values

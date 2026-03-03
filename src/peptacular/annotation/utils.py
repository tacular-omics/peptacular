from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any, TypeVar, overload

from tacular import (
    ELEMENT_LOOKUP,
    FRAGMENT_ION_LOOKUP,
    NEUTRAL_DELTA_LOOKUP,
    ElementInfo,
    FragmentIonInfo,
    IonType,
    IonTypeLiteral,
)

from ..constants import ELECTRON_MASS
from ..proforma_components.comps import ChargedFormula, GlobalChargeCarrier
from .cached_comps import DeltaInfo, IsotopeInfo
from .frag import Fragment
from .mod import Mods

H_ELEMENT_INFO = ELEMENT_LOOKUP["H"]


def adjust_mass_mz(
    base: float | Counter[ElementInfo],
    charge: Mods[GlobalChargeCarrier],
    ion_type: IonType | IonTypeLiteral | FragmentIonInfo,
    monoisotopic: bool,
    isotope: IsotopeInfo,
    delta: DeltaInfo,
    position: int | tuple[int, int] | None,
    parent_sequence: str,
    parent_sequence_length: int,
    internal_charge: int = 0,  # already includded in base mass (only effects mz calculation)
) -> Fragment:
    """Adjust base mass by charge carriers and ion type."""

    base_mass = 0.0
    if isinstance(base, Counter):
        for elem, count in base.items():
            base_mass += elem.get_mass(monoisotopic=monoisotopic) * count
    else:
        base_mass = base

    # Apply user specified isotopes
    base_mass += isotope.get_mass_delta(monoisotopic)

    # Apply losses
    base_mass += delta.get_mass_delta(monoisotopic)

    total_charge = charge.get_charge() + internal_charge
    # Correct for electron mass based on charge
    base_mass += charge.get_mass(monoisotopic)

    ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type] if not isinstance(ion_type, FragmentIonInfo) else ion_type
    base_mass += ion_info.get_mass(monoisotopic=monoisotopic)

    base_mass -= total_charge * ELECTRON_MASS

    # get adducts only if not protonated (charge_state == 0 means not protonated, None means protonated)
    if all(m.value.is_protonated for m in charge.mods):
        adducts = None
    else:
        adducts = tuple(m for m in charge._mods.keys()) if charge._mods else None

    return Fragment(
        ion_type=ion_info.ion_type,
        position=position,
        mass=base_mass,
        monoisotopic=monoisotopic,
        charge_state=total_charge,
        charge_adducts=adducts,
        isotopes=isotope.to_fragment_mapping,
        deltas=delta.to_fragment_mapping,
        composition=None,
        parent_sequence=parent_sequence,
        parent_sequence_length=parent_sequence_length,
    )


def adjust_comp(
    base_comp: Counter[ElementInfo],
    charge: Mods[GlobalChargeCarrier],
    ion_type: IonType | IonTypeLiteral | FragmentIonInfo,
    monoisotopic: bool,
    isotope: IsotopeInfo,
    delta: DeltaInfo,
    parent_sequence: str,
    parent_sequence_length: int,
    position: int | tuple[int, int] | None,
    inplace: bool = True,
    isotope_map: dict[ElementInfo, ElementInfo] | None = None,
    internal_charge: int = 0,
) -> Fragment:
    """Adjust base composition by charge carriers and ion type, returning a Fragment object."""

    if not inplace:
        base_comp = base_comp.copy()

    # corrects base_comp for user specified isotopes
    isotope.adjust_composition(base_comp)

    # Build fragment notation for losses
    delta.adjust_composition(base_comp)

    ion_info = FRAGMENT_ION_LOOKUP[ion_type] if not isinstance(ion_type, FragmentIonInfo) else ion_type

    base_comp += ion_info.composition

    # correct for global isotopes
    if isotope_map:
        for original_element, replaced_element in isotope_map.items():
            if original_element in base_comp:
                count = base_comp.pop(original_element)
                base_comp[replaced_element] += count

    for mod in charge.mods:
        base_comp += mod.get_composition()

    if any(count < 0 for count in base_comp.values()):
        raise ValueError(f"Negative element counts after charge adjustments: {base_comp}")

    total_charge = charge.get_charge() + internal_charge

    # Validate no negative counts
    if any(count < 0 for count in base_comp.values()):
        raise ValueError(f"Negative element counts after adjustments: {base_comp}")

    # Calculate mass from final composition
    base_mass = 0.0
    for elem, count in base_comp.items():
        base_mass += elem.get_mass(monoisotopic=monoisotopic) * count

    # Correct for electron mass based on charge
    if total_charge != 0:
        base_mass -= total_charge * ELECTRON_MASS

    # get adducts only if not protonated (charge_state == 0 means not protonated, None means protonated)
    if all(m.value.is_protonated for m in charge.mods):
        adducts = None
    else:
        adducts = tuple(m for m in charge._mods.keys()) if charge._mods else None

    return Fragment(
        ion_type=ion_info.ion_type,
        position=position,
        mass=base_mass,
        charge_state=total_charge,
        monoisotopic=monoisotopic,
        charge_adducts=adducts,
        isotopes=isotope.to_fragment_mapping,
        deltas=delta.to_fragment_mapping,
        composition=base_comp,
        parent_sequence=parent_sequence,
        parent_sequence_length=parent_sequence_length,
    )


def comp_frag(
    comp: Counter[ElementInfo],
    charge: Mods[GlobalChargeCarrier],
    ion_type: IonType | IonTypeLiteral | FragmentIonInfo,
    monoisotopic: bool,
    isotopes: IsotopeInfo,
    deltas: DeltaInfo,
    parent_sequence: str,
    parent_sequence_length: int,
    position: int | tuple[int, int] | None,
) -> Fragment:
    mass: float = 0
    for elem, count in comp.items():
        mass += elem.get_mass(monoisotopic=monoisotopic) * count

    return adjust_mass_mz(
        base=mass,
        charge=charge,
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotope=isotopes,
        delta=deltas,
        position=position,
        parent_sequence=parent_sequence,
        parent_sequence_length=parent_sequence_length,
    )


def process_losses(
    losses: str | ChargedFormula | float | Mapping[str | ChargedFormula | float, int],
) -> Mapping[ChargedFormula | float, int]:
    """Convert loss input to composition counter."""
    if isinstance(losses, str):
        return {ChargedFormula.from_composition(NEUTRAL_DELTA_LOOKUP[losses].composition): 1}
    if isinstance(losses, ChargedFormula):
        return {losses: 1}
    if isinstance(losses, float):
        return {losses: 1}

    if not isinstance(losses, dict):
        raise TypeError(f"Invalid losses type: {type(losses)}")

    # dict case
    total: Mapping[ChargedFormula | float, int] = Counter()
    for key, count in losses.items():
        if isinstance(key, str):
            try:
                loss = ChargedFormula.from_composition(NEUTRAL_DELTA_LOOKUP[key].composition)
                total[loss] += count
            except KeyError as e:
                loss = ChargedFormula.from_string(key, require_formula_prefix=False)
                if loss.charge:
                    raise ValueError(f"Loss formula cannot have charge: {key}") from e
                total[loss] += count
        elif isinstance(key, ChargedFormula):
            total[key] += count
        elif isinstance(key, float):
            total[key] += count
        else:
            raise TypeError(f"Invalid key type for loss: {type(key)}")

    return total


T = TypeVar("T", float, Counter[Any])


@overload
def cumsum(numbers: Sequence[float], reverse: bool = False) -> list[float]: ...


@overload
def cumsum(numbers: Sequence[Counter[Any]], reverse: bool = False) -> list[Counter[Any]]: ...


def cumsum(numbers: Sequence[float] | Sequence[Counter[Any]], reverse: bool = False) -> list[float] | list[Counter[Any]]:
    """Compute cumulative sum of a list of numbers or Counters."""
    match numbers:
        case [Counter(), *_]:
            # Counter case
            total: Counter[Any] = Counter()
            result: list[Counter[Any]] = []
            for counter in reversed(numbers) if reverse else numbers:
                total = total + counter  # type: ignore
                result.append(total.copy())
            return result

        case [float() | int(), *_] | []:
            # Numeric case
            total_num = 0.0
            result_num: list[float] = []
            for number in reversed(numbers) if reverse else numbers:
                total_num += number  # type: ignore
                result_num.append(total_num)
            return result_num

        case _:
            raise TypeError(f"cumsum expects sequence of float or Counter, got {type(numbers[0])}")


# Define rules: ion_type -> (position, required_aas, excluded_aas, specific_ion_map)
FRAGMENT_RULES: Any = {
    IonType.D: ("end", None, {"G", "A", "P", "I", "T"}, {"V": IonType.D_VALINE}),
    IonType.DA: (
        "end",
        {"I", "T"},
        None,
        {"I": IonType.DA_ISOLEUCINE, "T": IonType.DA_THREONINE},
    ),
    IonType.DB: (
        "start",
        {"I", "T"},
        None,
        {"I": IonType.DB_ISOLEUCINE, "T": IonType.DB_THREONINE},
    ),
    IonType.D_VALINE: ("end", {"V"}, None, None),
    IonType.DA_THREONINE: ("end", {"T"}, None, None),
    IonType.DA_ISOLEUCINE: ("end", {"I"}, None, None),
    IonType.DB_THREONINE: ("start", {"T"}, None, None),
    IonType.DB_ISOLEUCINE: ("start", {"I"}, None, None),
    IonType.W: ("start", None, {"G", "A", "P", "I", "T"}, {"V": IonType.W_VALINE}),
    IonType.WA: (
        "start",
        {"I", "T"},
        None,
        {"I": IonType.WA_ISOLEUCINE, "T": IonType.WA_THREONINE},
    ),
    IonType.WB: (
        "end",
        {"I", "T"},
        None,
        {"I": IonType.WB_ISOLEUCINE, "T": IonType.WB_THREONINE},
    ),
    IonType.W_VALINE: ("start", {"V"}, None, None),
    IonType.WA_THREONINE: ("start", {"T"}, None, None),
    IonType.WA_ISOLEUCINE: ("start", {"I"}, None, None),
    IonType.WB_THREONINE: ("end", {"T"}, None, None),
    IonType.WB_ISOLEUCINE: ("end", {"I"}, None, None),
}


def can_fragment_sequence(sequence: str, ion_type: IonType | IonTypeLiteral) -> IonType:
    """Check if a sequence can produce a fragment of the given ion type."""

    if not isinstance(ion_type, IonType):
        ion_type = IonType(ion_type)

    first_aa = sequence[0]
    last_aa = sequence[-1]

    if ion_type not in FRAGMENT_RULES:
        return ion_type

    position, required, excluded, specific_map = FRAGMENT_RULES[ion_type]
    aa = last_aa if position == "end" else first_aa

    if excluded and aa in excluded:
        raise ValueError(f"{ion_type.name} fragments cannot be produced from sequences {position}ing in {aa}.")

    if required and aa not in required:
        raise ValueError(f"{ion_type.name} fragments can only be produced from sequences {position}ing in {', or '.join(required)}.")

    if specific_map and aa in specific_map:
        return specific_map[aa]

    return ion_type

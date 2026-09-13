from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property, lru_cache
from math import isfinite
from typing import Self, cast

from tacular import ELEMENT_LOOKUP, NEUTRAL_DELTA_LOOKUP, ElementInfo

from ..diagnostics import CompositionError, InvalidAdjustmentError
from ..proforma_components import ChargedFormula, GlobalChargeCarrier

# ============================================================================
# IsotopeInfo
# ============================================================================


def _validate_isotope_count(count: int) -> None:
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise InvalidAdjustmentError("Isotope count must be an integer and cannot be negative")


@dataclass(frozen=True)
class IsotopeInfo:
    data: tuple[tuple[ElementInfo, int], ...]

    def __post_init__(self) -> None:
        for element, count in self.data:
            _validate_isotope_count(count)
            if element.mass_number is None:
                raise InvalidAdjustmentError(f"Specify an isotope such as 13C, not {element}")

    @staticmethod
    def from_input(
        isotopes: int | dict[ElementInfo | str, int] | None,
    ) -> Self:
        return cast(Self, get_isotopes(isotopes))

    @cached_property
    def _get_adjustments(self) -> tuple[tuple[ElementInfo, int], ...]:
        data = []
        for iso_elem_info, count in self.data:
            mono_elem_info = ELEMENT_LOOKUP[iso_elem_info.symbol]
            data.append((mono_elem_info, -count))
            data.append((iso_elem_info, count))
        return tuple(data)

    def adjust_composition(self, base_comp: Counter[ElementInfo]) -> None:
        for elem_info, count in self._get_adjustments:
            base_comp[elem_info] += count

        if any(v < 0 for v in base_comp.values()):
            negative_counts = {str(k): v for k, v in base_comp.items() if v < 0}
            raise InvalidAdjustmentError(f"Isotopic adjustment resulted in negative element counts: {negative_counts}")

    @property
    def composition(self) -> Counter[ElementInfo]:
        composition: Counter[ElementInfo] = Counter()
        for elem_info, count in self._get_adjustments:
            composition[elem_info] += count
        return composition

    @cached_property
    def monoisotopic_mass_delta(self) -> float:
        delta = 0.0
        for elem_info, count in self._get_adjustments:
            delta += elem_info.get_mass(monoisotopic=True) * count
        return delta

    @cached_property
    def average_mass_delta(self) -> float:
        delta = 0.0
        for elem_info, count in self._get_adjustments:
            delta += elem_info.get_mass(monoisotopic=False) * count
        return delta

    def get_mass_delta(self, monoisotopic: bool = False) -> float:
        if monoisotopic:
            return self.monoisotopic_mass_delta
        return self.average_mass_delta

    @property
    def to_fragment_mapping(self) -> Mapping[str, int] | int | None:
        # Plain property (not cached): returns a fresh dict so callers can't corrupt a
        # shared instance held by the @cache singleton factory.
        if len(self.data) == 1:
            elem_info, count = self.data[0]
            if str(elem_info) == "13C":
                return count
        return {str(elem_info): count for elem_info, count in self.data}


# ============================================================================
# ChargeCarrierInfo
# ============================================================================


@lru_cache(maxsize=1024)
def _get_charge_carrier_info(adducts: tuple[GlobalChargeCarrier, ...]) -> "ChargeCarrierInfo":
    """Cached factory for ChargeCarrierInfo."""
    return ChargeCarrierInfo(adducts)


@dataclass(frozen=True)
class ChargeCarrierInfo:
    adducts: tuple[GlobalChargeCarrier, ...]

    @cached_property
    def charge(self) -> int:
        return sum(adduct.get_charge() for adduct in self.adducts)

    def get_mass(self, monoisotopic: bool = False) -> float:
        if monoisotopic:
            return self.monoisotopic_mass
        return self.average_mass

    @cached_property
    def monoisotopic_mass(self) -> float:
        return sum(adduct.get_mass(monoisotopic=True) for adduct in self.adducts)

    @cached_property
    def average_mass(self) -> float:
        return sum(adduct.get_mass(monoisotopic=False) for adduct in self.adducts)

    @property
    def composition(self) -> Counter[ElementInfo]:
        # Plain property (not cached): this instance is a @cache singleton, so a cached
        # mutable Counter would be shared and could be corrupted by any caller.
        composition: Counter[ElementInfo] = Counter()
        for adduct in self.adducts:
            for element, count in adduct.get_composition().items():
                composition[element] += count
        return composition

    def adjust_composition(self, base_comp: Counter[ElementInfo]) -> None:
        for adduct in self.adducts:
            adduct_comp = adduct.get_composition()
            for elem_info, count in adduct_comp.items():
                base_comp[elem_info] += count

        if any(v < 0 for v in base_comp.values()):
            raise ValueError(f"Charge carrier adjustment resulted in negative element counts: {base_comp}")

    @staticmethod
    def from_input(
        charge: int | str | GlobalChargeCarrier | tuple[GlobalChargeCarrier, ...] | list[str] | None,
    ) -> Self:
        carriers = handle_charge_input(charge)
        # Sort for consistent caching - use serialized form as sort key
        carriers = tuple(sorted(carriers, key=lambda x: x.serialize()))
        return cast(Self, _get_charge_carrier_info(carriers))

    @property
    def to_fragment_mapping(self) -> Mapping[str, int] | None:
        composition: Counter[ElementInfo] = self.composition
        if not composition:
            return None
        if len(composition) == 1:
            elem_info, count = next(iter(composition.items()))
            if str(elem_info) == "H":
                return None

        return {str(elem_info): count for elem_info, count in composition.items()}

    @property
    def to_explicit_fragment_mapping(self) -> Mapping[str, int]:
        composition: Counter[ElementInfo] = self.composition
        if not composition:
            return {}
        return {str(elem_info): count for elem_info, count in composition.items()}

    @property
    def to_proforma_charge(self) -> Mapping[str, int] | int | None:
        composition: Counter[ElementInfo] = self.composition
        if not composition:
            return None
        if len(composition) == 1:
            elem_info, count = next(iter(composition.items()))
            if str(elem_info) == "H":
                return count

        return {str(elem_info): count for elem_info, count in composition.items()}

    @cached_property
    def to_tuple_or_int(self) -> tuple[str, ...] | int | None:
        composition: Counter[ElementInfo] = self.composition
        if not composition:
            return None

        elif len(composition) == 1:
            elem_info, count = next(iter(composition.items()))
            if str(elem_info) == "H":
                return count

        return tuple(str(elem_info) for elem_info in composition.keys())


# ============================================================================
# DeltaInfo
# ============================================================================


@lru_cache(maxsize=1024)
def _normalize_delta_key(key_str: str) -> ChargedFormula | float:
    """Normalize a delta key string to ChargedFormula or float."""
    # Try as named delta first
    try:
        return ChargedFormula.from_composition(NEUTRAL_DELTA_LOOKUP[key_str].composition)
    except KeyError:
        pass

    # Try as float
    try:
        return float(key_str)
    except ValueError:
        pass

    # Parse as formula
    return ChargedFormula.from_string(key_str, require_formula_prefix=False)


@lru_cache(maxsize=1024)
def _get_delta_info(items: tuple[tuple[str, int], ...] | None) -> "DeltaInfo":
    """Cached factory for DeltaInfo."""
    if items is None or len(items) == 0:
        return DeltaInfo({})

    result: dict[ChargedFormula | float, int] = {}
    for key_str, count in items:
        key = _normalize_delta_key(key_str)
        result[key] = result.get(key, 0) + count

    # Remove zero counts
    result = {k: v for k, v in result.items() if v != 0}

    return DeltaInfo(result)


@dataclass(frozen=True, init=False)
class DeltaInfo:
    _items: tuple[tuple[ChargedFormula | float, int], ...]

    def __init__(self, deltas: Mapping[ChargedFormula | float, int]):
        for key, count in deltas.items():
            if isinstance(count, bool) or not isinstance(count, int):
                raise InvalidAdjustmentError("Delta counts must be integers")
            if isinstance(key, ChargedFormula):
                if key.is_charged:
                    raise InvalidAdjustmentError("Delta formulas must be neutral (charge=0)")
            elif isinstance(key, bool) or not isinstance(key, (int, float)) or not isfinite(key):
                raise InvalidAdjustmentError("Delta masses must be finite real numbers")
        object.__setattr__(self, "_items", tuple((float(key) if isinstance(key, int) else key, count) for key, count in deltas.items()))

    @property
    def deltas(self) -> Mapping[ChargedFormula | float, int]:
        """Return a fresh mapping so cached instances cannot be changed by callers."""
        return dict(self._items)

    @property
    def has_floats(self) -> bool:
        return any(isinstance(k, float) for k in self.deltas.keys())

    @cached_property
    def monoisotopic_mass_delta(self) -> float:
        delta = 0.0
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                delta += key.get_mass(monoisotopic=True) * count
            else:
                delta += key * count
        return delta

    @cached_property
    def average_mass_delta(self) -> float:
        delta = 0.0
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                delta += key.get_mass(monoisotopic=False) * count
            else:
                delta += key * count
        return delta

    def get_mass_delta(self, monoisotopic: bool = False) -> float:
        if monoisotopic:
            return self.monoisotopic_mass_delta
        return self.average_mass_delta

    @property
    def composition(self) -> Counter[ElementInfo]:
        if self.has_floats:
            raise CompositionError("Cannot get composition when deltas include floats. Use mass() or mz() for mass-only adjustments.")
        composition: Counter[ElementInfo] = Counter()
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                key_comp = key.get_composition()
                for elem_info, elem_count in key_comp.items():
                    composition[elem_info] += elem_count * count
        return composition

    def adjust_composition(self, base_comp: Counter[ElementInfo]) -> None:
        if self.has_floats:
            raise CompositionError("Cannot adjust composition when deltas include floats. Use mass() or mz() for mass-only adjustments.")
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                key_comp = key.get_composition()
                for elem_info, elem_count in key_comp.items():
                    base_comp[elem_info] += elem_count * count

        if any(v < 0 for v in base_comp.values()):
            raise InvalidAdjustmentError(f"Delta adjustment resulted in negative element counts: {base_comp}")

    @staticmethod
    def from_input(
        deltas: str | ChargedFormula | float | dict[str | ChargedFormula | float, int] | None,
    ) -> Self:
        if deltas is None:
            return cast(Self, _get_delta_info(None))

        # Normalize to dict first
        normalized_dict = _handle_delta_input(deltas)

        # Check for charged formulas
        for key in normalized_dict.keys():
            if isinstance(key, ChargedFormula) and key.is_charged:
                raise ValueError("Delta formulas must be neutral (charge=0)")

        # Convert to sorted tuple of (str_key, count) for caching
        items: list[tuple[str, int]] = []
        for k, v in normalized_dict.items():
            if isinstance(k, float):
                key_str = str(k)
            elif isinstance(k, ChargedFormula):
                key_str = k.serialize(include_formula_prefix=False)
            else:
                key_str = str(k)
            items.append((key_str, v))

        # Sort for consistent caching
        items_tuple = tuple(sorted(items, key=lambda x: (x[0], x[1])))

        return cast(Self, _get_delta_info(items_tuple))

    def __add__(self, other: "DeltaInfo") -> "DeltaInfo":
        """Combine two DeltaInfo objects (add deltas together)."""
        combined: dict[ChargedFormula | float, int] = dict(self.deltas)
        for key, count in other.deltas.items():
            combined[key] = combined.get(key, 0) + count

        # Remove zero counts
        combined = {k: v for k, v in combined.items() if v != 0}

        # Use from_input to get cached version
        return DeltaInfo.from_input(combined)  # type: ignore

    def __sub__(self, other: "DeltaInfo") -> "DeltaInfo":
        """Subtract one DeltaInfo from another."""
        result: dict[ChargedFormula | float, int] = dict(self.deltas)
        for key, count in other.deltas.items():
            result[key] = result.get(key, 0) - count

        # Remove zero counts
        result = {k: v for k, v in result.items() if v != 0}

        # Use from_input to get cached version
        return DeltaInfo.from_input(result)  # type: ignore

    def __neg__(self) -> "DeltaInfo":
        """Negate all deltas (flip signs)."""
        negated = {k: -v for k, v in self.deltas.items()}
        return DeltaInfo.from_input(negated)

    def __mul__(self, scalar: int) -> "DeltaInfo":
        """Multiply all delta counts by a scalar."""
        multiplied = {k: v * scalar for k, v in self.deltas.items()}
        return DeltaInfo.from_input(multiplied)

    def __rmul__(self, scalar: int) -> "DeltaInfo":
        """Allow scalar * DeltaInfo."""
        return self.__mul__(scalar)

    def __str__(self) -> str:
        parts = []
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                part = f"Formula:{str(key)}"
            else:
                part = str(key)
            if count != 1:
                part += f" x{count}"
            parts.append(part)
        return " + ".join(parts) if parts else "No Delta"

    @property
    def to_fragment_mapping(self) -> Mapping[str | float, int] | None:
        if not self.deltas:
            return None
        result: dict[str | float, int] = {}
        for key, count in self.deltas.items():
            if isinstance(key, ChargedFormula):
                result[key.serialize(include_formula_prefix=False)] = count
            else:
                result[key] = count
        return result


# ============================================================================
# Helper functions
# ============================================================================


def _handle_delta_input(
    deltas: str | ChargedFormula | float | dict[str | ChargedFormula | float, int] | None,
) -> dict[ChargedFormula | float, int]:
    """Normalize delta input to standard dict format."""
    if deltas is None:
        return {}

    if isinstance(deltas, dict):
        result: dict[ChargedFormula | float, int] = {}
        for key, count in deltas.items():
            if isinstance(count, bool) or not isinstance(count, int):
                raise InvalidAdjustmentError("Delta counts must be integers")
            # Recurse on each key to normalize it
            normalized = _handle_delta_input(key)  # Returns single-item dict
            for k, v in normalized.items():
                result[k] = result.get(k, 0) + (v * count)
        return result

    # Single items - base cases
    if isinstance(deltas, str):
        try:
            return {ChargedFormula.from_composition(NEUTRAL_DELTA_LOOKUP[deltas].composition): 1}
        except KeyError:
            pass

        return {ChargedFormula.from_string(deltas, require_formula_prefix=False): 1}

    if isinstance(deltas, (float, int)):
        if isinstance(deltas, bool) or not isfinite(deltas):
            raise InvalidAdjustmentError("Delta masses must be finite real numbers")
        return {float(deltas): 1}

    if isinstance(deltas, ChargedFormula):
        return {deltas: 1}

    raise TypeError(f"Invalid delta type: {type(deltas)}")


def handle_charge_input(
    charge: int | str | GlobalChargeCarrier | tuple[GlobalChargeCarrier | str, ...] | list[str] | None,
) -> tuple[GlobalChargeCarrier, ...]:
    """Normalize charge input to tuple of GlobalChargeCarrier."""
    if isinstance(charge, int):
        adduct: GlobalChargeCarrier = get_charge_adducts(charge)
        return (adduct,)
    elif isinstance(charge, str):
        adduct = GlobalChargeCarrier.from_string(charge)
        return (adduct,)
    elif isinstance(charge, GlobalChargeCarrier):
        return (charge,)
    elif isinstance(charge, tuple):
        adducts = []
        for c in charge:
            adducts.extend(handle_charge_input(c))
        return tuple(adducts)
    elif charge is None:
        return ()
    elif isinstance(charge, list):
        adducts = []
        for adduct_str in charge:
            cformula = ChargedFormula.from_string(adduct_str, require_formula_prefix=False)
            adduct = GlobalChargeCarrier(cformula, occurance=1)
            adducts.extend([adduct] * 1)
        return tuple(adducts)
    else:
        raise TypeError(f"Invalid charge type: {type(charge)}")


@lru_cache(maxsize=1024)
def get_charge_adducts(charge_state: int) -> GlobalChargeCarrier:
    """Get a cached GlobalChargeCarrier for a given charge state."""
    return GlobalChargeCarrier.charged_proton(charge_state)


C13: ElementInfo = ELEMENT_LOOKUP["13C"]


@lru_cache(maxsize=1024)
def _get_isotopes(
    isotopes: int | tuple[tuple[ElementInfo, int], ...],
) -> IsotopeInfo:
    """Cached factory for IsotopeInfo."""
    if isinstance(isotopes, tuple):
        return IsotopeInfo(isotopes)

    if isotopes == 0:
        return IsotopeInfo(())

    return IsotopeInfo(((C13, isotopes),))


def get_isotopes(
    isotopes: int | dict[ElementInfo | str, int] | None,
) -> IsotopeInfo:
    """Normalize isotope input and return cached IsotopeInfo."""
    if isotopes is None:
        return _get_isotopes(0)

    if isinstance(isotopes, dict):
        elem_infos: list[tuple[ElementInfo, int]] = []
        for elem, count in isotopes.items():
            # Validate before cache lookup because True and 1.0 compare equal to 1.
            _validate_isotope_count(count)
            if isinstance(elem, str):
                elem_info: ElementInfo = ELEMENT_LOOKUP[elem]
            else:
                elem_info = elem
            elem_infos.append((elem_info, count))

        # Sort for consistent caching
        elem_infos_sorted = tuple(sorted(elem_infos, key=lambda x: str(x[0])))
        return _get_isotopes(elem_infos_sorted)

    if isinstance(isotopes, int):
        if isinstance(isotopes, bool):
            raise InvalidAdjustmentError("Isotope count must be an integer, not bool")
        if isotopes < 0:
            raise ValueError("Isotope count cannot be negative")
        return _get_isotopes(isotopes)

    raise TypeError(f"Invalid isotope type: {type(isotopes)}")


@lru_cache(maxsize=1024)
def _get_losses(
    losses: tuple[tuple[ElementInfo, int], ...],
) -> Mapping[ChargedFormula, int]:
    """Cached factory for neutral loss mappings."""
    composition: Counter[ElementInfo] = Counter()
    for element, count in losses:
        composition[element] += count
    charged_formula = ChargedFormula.from_composition(composition)
    return {charged_formula: 1}


def get_losses(
    losses: Counter[ElementInfo],
) -> Mapping[ChargedFormula, int] | None:
    """Get cached neutral loss mapping from composition."""
    # Return a fresh dict: _get_losses is @cache-backed and would otherwise hand out a
    # shared mutable mapping that a caller could corrupt for all future callers.
    return dict(_get_losses(tuple(losses.items()))) if losses else None

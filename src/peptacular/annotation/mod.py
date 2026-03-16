import sys
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from typing import Any, Protocol, Self, cast

from tacular import AA_LOOKUP, ElementInfo

from ..constants import ModType
from ..proforma_components import (
    FixedModification,
    GlobalChargeCarrier,
    IsotopeReplacement,
    MassPropertyMixin,
    ModificationTags,
    TagMass,
)

# Define your modification types
ModValue = IsotopeReplacement | FixedModification | GlobalChargeCarrier | ModificationTags


class ModificationProtocol(Protocol):
    """Protocol defining the interface all modifications must implement."""

    def get_mass(self, monoisotopic: bool = True) -> float: ...

    def get_composition(self) -> Counter[ElementInfo]: ...

    def get_charge(self) -> int | None: ...

    def validate(self) -> str | None: ...


# Use bound instead of constraints


@dataclass(frozen=True, slots=True)
class Mod[T: ModificationProtocol]:
    """A modification with its occurrence count."""

    value: T
    count: int

    def __post_init__(self):
        if self.count < 0:
            raise ValueError(f"Count must be non-negative, got {self.count}")

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        if self.count < 0:
            return f"Count must be non-negative, got {self.count}"

        if hasattr(self.value, "validate"):
            return self.value.validate()

        return None

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get total mass for this modification occurrence."""
        mass: int | float = self.value.get_mass(monoisotopic)
        return mass * self.count

    def get_composition(self) -> Counter[ElementInfo]:
        """Get total composition for this modification occurrence."""
        return self.value.get_composition().copy()

    def get_charge(self) -> int:
        """Get total charge for this modification occurrence."""
        charge: int | None = self.value.get_charge()
        if charge is None:
            return 0

        return charge * self.count

    def as_tuple(self) -> tuple[T, int]:
        """Return the modification as a (value, count) tuple."""
        return (self.value, self.count)

    def copy(self) -> Self:
        return self.__class__(
            value=self.value,
            count=self.count,
        )


# At module level, define the parser mapping
_MOD_PARSERS: dict[ModType, Callable[[str], Any]] = {
    ModType.ISOTOPE: IsotopeReplacement.from_string,
    ModType.STATIC: FixedModification.from_string,
    ModType.LABILE: ModificationTags.from_string,
    ModType.UNKNOWN: ModificationTags.from_string,
    ModType.NTERM: ModificationTags.from_string,
    ModType.CTERM: ModificationTags.from_string,
    ModType.CHARGE: GlobalChargeCarrier.from_string,
    ModType.INTERNAL: ModificationTags.from_string,
    ModType.INTERVAL: ModificationTags.from_string,
}


@dataclass(frozen=True)
class Mods[T: ModificationProtocol](MassPropertyMixin):
    """Collection of modifications of a specific type."""

    mod_type: ModType
    _mods: dict[str, int] | None

    def __post_init__(self):
        """Validate mod_type is supported."""
        if self.mod_type not in _MOD_PARSERS:
            raise ValueError(f"Unsupported mod_type: {self.mod_type}")

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        if self.mod_type not in _MOD_PARSERS:
            return f"Unsupported mod_type: {self.mod_type}"

        if self._mods is None:
            return None

        for mod in self.mods:
            error = mod.validate()
            if error is not None:
                return error

        return None

    def _parse_mod(self, mod_str: str) -> T:
        """Parse a modification string based on mod_type."""
        parser = _MOD_PARSERS[self.mod_type]
        # We've validated mod_type in __post_init__, so this cast is safe
        return cast(T, parser(mod_str))

    def parse_items(self) -> Iterable[tuple[T, int]]:
        """Get raw modification items as (mod, count) tuples."""
        if self._mods is None:
            return ()
        return ((mod.value, mod.count) for mod in self.mods)

    def parse_tuples(self) -> Iterable[tuple[T, int]]:
        """Get raw modification items as (mod, count) tuples."""
        if self._mods is None:
            return ()
        return ((mod.value, mod.count) for mod in self.mods)

    @property
    def mods(self) -> tuple[Mod[T], ...]:
        """Parse stored modifications into typed Mod objects."""

        if self._mods is None:
            return tuple()

        return tuple(Mod(value=self._parse_mod(mod_str), count=count) for mod_str, count in self._mods.items())

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get total mass for all modifications."""
        return sum(mod.get_mass(monoisotopic) for mod in self.mods)

    def get_composition(self) -> Counter[ElementInfo]:
        """Get total composition for all modifications."""
        return sum((mod.get_composition() for mod in self.mods), Counter())

    def get_composition_with_delta_mass_charge(self, monoisotopic: bool = True) -> tuple[Counter[ElementInfo], float, int]:
        """Get total composition and when not possible fall back to delta mass for MassTags."""

        total_composition = Counter[ElementInfo]()
        total_delta_mass = 0.0
        total_charge = 0
        for mod in self.mods:
            total_charge += mod.get_charge()
            try:
                comp = mod.get_composition()
                total_composition += comp
            except ValueError as e:
                if isinstance(mod.value, ModificationTags) and isinstance(mod.value.first_tag, TagMass):
                    total_delta_mass += mod.get_mass(monoisotopic=monoisotopic)
                else:
                    raise ValueError(f"Cannot get composition for mod {mod.value}, and no mass tag found.") from e
        return total_composition, total_delta_mass, total_charge

    def get_charge(self) -> int:
        """Get total charge for all modifications."""
        return sum(mod.get_charge() for mod in self.mods)

    def get_mass_charge(self, monoisotopic: bool = True) -> tuple[float, int]:
        """Get total mass and charge for all modifications."""
        mods = self.mods

        if len(mods) == 1:
            return mods[0].get_mass(monoisotopic), mods[0].get_charge()

        total_mass = sum(mod.get_mass(monoisotopic) for mod in mods)
        total_charge = sum(mod.get_charge() for mod in mods)
        return total_mass, total_charge

    def __len__(self) -> int:
        return len(self._mods) if self._mods is not None else 0

    def __iter__(self) -> Iterator[Mod[T]]:
        return iter(self.mods)

    def __contains__(self, item: Any) -> bool:
        if not isinstance(item, str):
            item = str(item)

        if self._mods is None:
            return False
        return item in self._mods

    def serialize(self) -> str:
        """Serialize modifications to string format."""
        mod_str_comps: list[str] = []

        match self.mod_type:
            case ModType.ISOTOPE:
                for mod_str, count in (self._mods or {}).items():
                    for _ in range(count):
                        mod_str_comps.append(f"<{mod_str}>")
            case ModType.STATIC:
                for mod_str, count in (self._mods or {}).items():
                    for _ in range(count):
                        mod_str_comps.append(f"<{mod_str}>")
            case ModType.LABILE:
                for mod_str, count in (self._mods or {}).items():
                    for _ in range(count):
                        mod_str_comps.append(f"{{{mod_str}}}")
            case ModType.UNKNOWN:
                for mod_str, count in (self._mods or {}).items():
                    if count == 1:
                        mod_str_comps.append(f"[{mod_str}]")
                    else:
                        mod_str_comps.append(f"[{mod_str}]^{count}")
                else:
                    mod_str_comps.append("?")
            case ModType.NTERM:
                for mod_str, count in (self._mods or {}).items():
                    for _ in range(count):
                        mod_str_comps.append(f"[{mod_str}]")
                else:
                    mod_str_comps.append("-")
            case ModType.CTERM:
                prefix_added = False
                for mod_str, count in (self._mods or {}).items():
                    if not prefix_added:
                        mod_str_comps.append("-")
                        prefix_added = True
                    for _ in range(count):
                        mod_str_comps.append(f"[{mod_str}]")
            case ModType.CHARGE:
                for mod_str, cnt in (self._mods or {}).items():
                    if cnt != 1:
                        raise RuntimeError(f"Charge modifications cannot have count > 1, got {cnt} for mod {mod_str}")
                    mod_str_comps.append(f"{mod_str}")
                return f"[{','.join(mod_str_comps)}]"
            case ModType.INTERNAL:
                for mod_str, count in (self._mods or {}).items():
                    for _ in range(count):
                        mod_str_comps.append(f"[{mod_str}]")
            case ModType.INTERVAL:
                for mod_str, count in (self._mods or {}).items():
                    for _ in range(count):
                        mod_str_comps.append(f"[{mod_str}]")
            case _:
                raise ValueError(f"Unsupported mod_type: {self.mod_type}")
        return "".join(mod_str_comps)

    def __str__(self) -> str:
        if not self._mods:
            return f"Mods@{self.mod_type.value.capitalize()}()"

        # Show all modifications in compact form
        mods_repr = "[" + ", ".join(f"{k}^{v}" if v > 1 else f"{k}" for k, v in self._mods.items()) + "]"
        return f"Mods@{self.mod_type.value.capitalize()}({mods_repr})"

    def __repr__(self) -> str:
        if not self._mods:
            return f"Mods(mod_type={self.mod_type!r}, _mods=None)"

        # Show all modifications in compact form
        mods_repr = ", ".join(f"{k!r}: {v}" for k, v in self._mods.items())
        return f"Mods(mod_type={self.mod_type!r}, _mods={{{mods_repr}}})"

    def copy(self) -> Self:
        return self.__class__(
            mod_type=self.mod_type,
            _mods=self._mods.copy() if self._mods else None,
        )

    def __hash__(self) -> int:
        if self._mods is None:
            return hash((self.mod_type, None))
        # Create a frozenset of the items for hashing
        mods_frozenset = frozenset(self._mods.items())
        return hash((self.mod_type, mods_frozenset))


# Valid Amino Acid codes including Ambiguous ones (B, Z, J, X) and potential future codes
VALID_AMINO_ACIDS = {aa.one_letter_code for aa in AA_LOOKUP.ordered_amino_acids}


def condense_mod_str(
    mod_dict: dict[str, int],
    start_bracket: str,
    end_bracket: str,
    allow_mult: bool = False,
    prefix: str = "",
    suffix: str = "",
) -> str:
    mod_strs: list[str] = []
    for mod_str, count in mod_dict.items():
        if count == 1:
            mod_strs.append(f"{start_bracket}{mod_str}{end_bracket}")
        else:
            if allow_mult:
                mod_strs.append(f"{start_bracket}{mod_str}{end_bracket}^{count}")
            else:
                for _ in range(count):
                    mod_strs.append(f"{start_bracket}{mod_str}{end_bracket}")
    return prefix + "".join(mod_strs) + suffix


def convert_moddict_input(mod: Any) -> dict[str, int]:
    # Convert mod input to string representation if needed
    d: dict[str, int] = {}
    if isinstance(mod, dict) or isinstance(mod, Counter):
        # if value is not string, convert to string
        d = {sys.intern(str(k)): v for k, v in mod.items()}
    elif isinstance(mod, Mods):
        return convert_moddict_input(mod._mods)
    elif isinstance(mod, str):
        d = {sys.intern(str(mod)): 1}
    elif isinstance(mod, (int, float)):
        # ensure it has +/- in front of number
        num_str = sys.intern(f"{mod:+}")
        d = {num_str: 1}
    elif hasattr(mod, "__iter__"):
        for m in mod:
            if hasattr(m, "__iter__") and not isinstance(m, str):
                mod_str, count = convert_single_mod_input(m)
                d[mod_str] = d.get(mod_str, 0) + count
            else:
                d[sys.intern(str(m))] = d.get(str(m), 0) + 1
    return d


def convert_single_mod_input(mod: Any) -> tuple[str, int]:
    # Convert single mod input to string representation if needed
    if isinstance(mod, str):
        return sys.intern(str(mod)), 1
    elif isinstance(mod, (int, float)):
        # ensure it has +/- in front of number
        return sys.intern(f"{mod:+}"), 1
    elif isinstance(mod, tuple) and len(mod) == 2:
        return sys.intern(str(mod[0])), int(mod[1])
    elif isinstance(mod, Mod):
        return sys.intern(str(mod.value)), mod.count
    else:
        return sys.intern(str(mod)), 1


EMPTYP_INTERVAL_MODS = Mods[ModificationTags](mod_type=ModType.INTERVAL, _mods=None)


class Interval:
    __slots__ = ("_start", "_end", "_ambiguous", "_mods", "_validate")

    def __init__(
        self,
        start: int,
        end: int,
        ambiguous: bool = False,
        mods: Any | None = None,
        validate: bool = False,
    ):
        self._start = start
        self._end = end
        self._ambiguous = ambiguous
        self._validate = validate
        self._mods = None
        self.set_mods(mods, validate=validate)

        if self._start < 0:
            raise ValueError(f"Start position must be non-negative, got {self.start}")
        if self._end <= self.start:
            raise ValueError(f"End position must be >= start position, got {self.end} <= {self.start}")

    @property
    def is_valid(self) -> bool:
        return self.validate() is None

    def validate(self) -> str | None:
        if self.end < self.start:
            return f"End ({self.end}) must be >= start ({self.start})"

        return self.mods.validate()

    @property
    def mods(self) -> Mods[ModificationTags]:
        if self._mods is None:
            return EMPTYP_INTERVAL_MODS
        return Mods[ModificationTags](mod_type=ModType.INTERVAL, _mods=self._mods)

    @mods.setter
    def mods(self, value: Any) -> None:
        self.set_mods(value)

    def set_mods(
        self,
        mods: dict[Any, int] | Mods[ModificationTags] | None,
        validate: bool | None = None,
    ) -> None:
        if validate is None:
            validate = self._validate

        if isinstance(mods, Mods):
            mods = mods._mods

        if mods is None:
            self._mods = None
            return

        converted_mods = convert_moddict_input(mods)

        if len(converted_mods) == 0:
            self._mods = None
            return

        if validate:
            for mod_str in converted_mods.keys():
                ModificationTags.from_string(mod_str)

        self._mods = converted_mods

    def copy(self) -> Self:
        return self.__class__(
            start=self._start,
            end=self._end,
            ambiguous=self._ambiguous,
            mods=self._mods.copy() if self._mods else None,
            validate=self._validate,
        )

    def append_mod(self, mod: Any, validate: bool | None = None, inplace: bool = True) -> None:
        if not inplace:
            return self.copy().append_mod(mod, validate=validate, inplace=True)

        if validate is None:
            validate = self._validate

        mod_str, count = convert_single_mod_input(mod)

        if self._mods is None:
            self._mods = {}

        if validate:
            ModificationTags.from_string(mod_str)

        self._mods[mod_str] = self._mods.get(mod_str, 0) + count

    def extend_mods(self, mods: Any, validate: bool | None = None) -> None:
        if validate is None:
            validate = self._validate

        converted_mods = convert_moddict_input(mods)

        if self._mods is None:
            self._mods = {}

        if validate:
            for mod_str in converted_mods.keys():
                ModificationTags.from_string(mod_str)

        for mod_str, count in converted_mods.items():
            self._mods[mod_str] = self._mods.get(mod_str, 0) + count

    def remove_mods(self, mods: Any) -> None:
        converted_mods = convert_moddict_input(mods)

        if self._mods is None:
            return

        for mod_str, count in converted_mods.items():
            if mod_str in self._mods:
                self._mods[mod_str] -= count
                if self._mods[mod_str] <= 0:
                    del self._mods[mod_str]

        if len(self._mods) == 0:
            self._mods = None

    def pop_mods(self) -> Mods[ModificationTags]:
        mods = self.mods
        self._mods = None
        return mods

    @property
    def has_mods(self) -> bool:
        return self._mods is not None and len(self._mods) > 0

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    @property
    def ambiguous(self) -> bool:
        return self._ambiguous

    def update(self, **kwargs: Any) -> Self:
        return self.__class__(
            start=kwargs.get("start", self.start),
            end=kwargs.get("end", self.end),
            ambiguous=kwargs.get("ambiguous", self.ambiguous),
            mods=kwargs.get("mods", self._mods.copy() if self._mods else None),
            validate=kwargs.get("validate", self._validate),
        )

    @property
    def mod_str(self) -> str:
        if self._mods is None:
            return ""
        return condense_mod_str(self._mods, "[", "]", allow_mult=False)

    def __str__(self) -> str:
        return f"{self.start}-{self.end}[{self.mod_str}]"

    def __repr__(self) -> str:
        return f"Interval(start={self.start}, end={self.end}, ambiguous={self.ambiguous}, mods={self._mods})"

    def __eq__(self, value: object) -> bool:
        if not isinstance(value, Interval):
            return False
        return self.start == value.start and self.end == value.end and self.ambiguous == value.ambiguous and self._mods == value._mods

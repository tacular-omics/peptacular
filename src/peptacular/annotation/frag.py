from collections import Counter
from collections.abc import Mapping
from typing import Any

from tacular import (
    ELEMENT_LOOKUP,
    ElementInfo,
    IonType,
)

from ..constants import ModType
from ..proforma_components import (
    ChargedFormula,
    GlobalChargeCarrier,
)
from .mod import Mods
from .positions import validate_position


class Fragment:
    def __init__(
        self,
        ion_type: IonType,
        position: int | tuple[int, int] | None,
        mass: float,
        monoisotopic: bool,
        charge_state: int,
        charge_adducts: tuple[str, ...] | None = None,
        isotopes: Mapping[str, int] | int | None = None,
        deltas: Mapping[str | float, int] | None = None,
        composition: Mapping[ElementInfo, int] | None = None,
        parent_sequence: str | None = None,
        parent_sequence_length: int | None = None,
    ) -> None:
        self.ion_type: IonType = ion_type
        self.position: int | tuple[int, int] | None = position
        self.mass: int | float = mass
        self.monoisotopic: bool = monoisotopic
        self.charge_state: int = charge_state
        # If None and charge_state != 0: means protonated
        self._charge_adducts: tuple[str, ...] | None = charge_adducts
        # int means 13C count
        self._isotopes: Mapping[str, int] | int | None = isotopes
        self._losses: Mapping[str | float, int] | None = deltas
        # Optional composition cache
        self._composition: Counter[ElementInfo] | None = composition
        self.parent_sequence: str | None = parent_sequence
        self.parent_sequence_length: int | None = parent_sequence_length

    @property
    def composition(self) -> Counter[ElementInfo] | None:
        if self._composition is not None:
            return self._composition

        if self.parent_sequence is None:
            raise ValueError("Cannot calculate composition without parent sequence or explicit composition")

        from .annotation import ProFormaAnnotation

        annot = ProFormaAnnotation.parse(self.parent_sequence)
        return annot.comp(isotopes=self.isotopes, deltas=self.losses)  # type: ignore

    @property
    def mz(self) -> float:
        return self.mass / self.charge_state if self.charge_state != 0 else self.mass

    @property
    def neutral_mass(self) -> float:
        # subract adduct masses
        total_adduct_mass = 0.0
        for adduct in self.charge_adducts:
            total_adduct_mass += adduct.get_mass(self.monoisotopic)
        return self.mass - total_adduct_mass

    @property
    def charge_adducts(self) -> Mods[GlobalChargeCarrier]:
        if self._charge_adducts is None:
            if self.charge_state != 0:
                return Mods[GlobalChargeCarrier](
                    mod_type=ModType.CHARGE,
                    _mods={GlobalChargeCarrier.charged_proton(self.charge_state).serialize(): 1},
                )
            # no adducts no charge
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={})

        # we have adducts, convert to Mods object
        else:
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={k: 1 for k in self._charge_adducts})

    @property
    def is_protonated(self) -> bool:
        if self._charge_adducts is None and self.charge_state != 0:
            return True
        return False

    @property
    def isotopes(self) -> Mapping[ElementInfo, int]:
        if self._isotopes is not None:
            if isinstance(self._isotopes, int):
                return {ELEMENT_LOOKUP["13C"]: self._isotopes}

            if isinstance(self._isotopes, dict):
                return {ELEMENT_LOOKUP[elem]: count for elem, count in self._isotopes.items()}

        return {}

    @property
    def is_c13(self) -> bool:
        if self._isotopes is not None and isinstance(self._isotopes, int):
            return True
        return False

    @property
    def losses(self) -> Mapping[ChargedFormula | float, int]:
        if self._losses is not None:
            losses = {}
            for loss_name, count in self._losses.items():
                if isinstance(loss_name, float | int):
                    losses[loss_name] = count
                    continue
                loss_formula = ChargedFormula.from_string(loss_name, require_formula_prefix=False)
                losses[loss_formula] = count
            return losses
        return {}

    def asdict(self) -> dict[str, Any]:
        return {
            "ion_type": self.ion_type,
            "position": self.position,
            "mass": self.mass,
            "charge_state": self.charge_state,
            "monoisotopic": self.monoisotopic,
            "charge_adducts": self.charge_adducts,
            "isotopes": self.isotopes,
            "losses": self.losses,
        }

    def __str__(self) -> str:
        parts = []
        parts.append(f"ion_type={self.ion_type}")

        if self.position is not None:
            parts.append(f"position={self.position}")

        parts.append(f"mass={self.mass:.4f}")
        parts.append(f"charge={self.charge_state}")

        # Only show charge_adducts if not simple protonation
        if not self.is_protonated and self.charge_state != 0:
            parts.append(f"charge_adducts={self.charge_adducts}")

        if self.isotopes:
            parts.append(f"isotopes={dict(self.isotopes)}")

        if self.losses:
            parts.append(f"losses={dict(self.losses)}")

        return f"Fragment({', '.join(parts)})"

    def __repr__(self) -> str:
        return (
            f"Fragment(ion_type={self.ion_type!r}, position={self.position!r}, "
            f"mass={self.mass}, monoisotopic={self.monoisotopic}, "
            f"charge_state={self.charge_state}, charge_adducts={self._charge_adducts!r}, "
            f"isotopes={self._isotopes!r}, deltas={self._losses!r}, "
            f"composition={self._composition!r}, parent_sequence={self.parent_sequence!r}, "
            f"parent_sequence_length={self.parent_sequence_length})"
        )

    @property
    def sequence(self) -> str | None:
        if self.parent_sequence is None:
            raise ValueError("Cannot determine fragment sequence without parent sequence")

        if self.parent_sequence_length is None:
            raise ValueError("Cannot determine fragment sequence without parent sequence length")

        pos = validate_position(self.ion_type, self.position, self.parent_sequence_length)
        if pos is None:
            return self.parent_sequence
        elif isinstance(pos, tuple):
            from .annotation import ProFormaAnnotation

            start, end = pos
            return ProFormaAnnotation.parse(self.parent_sequence)[slice(start, end)].serialize()

        raise ValueError("Invalid position format for fragment sequence extraction")

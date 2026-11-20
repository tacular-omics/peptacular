"""Mass, composition and charge engine behind :class:`~peptacular.annotation.ProFormaAnnotation`.

Every function takes the annotation as its first argument. The annotation's ``mass()``,
``mz()``, ``comp()``, ``charge_adducts`` and related members are short calls into this
module; use those methods, not these functions.
"""

from collections import Counter
from typing import TYPE_CHECKING

from tacular import AA_LOOKUP, Element, ElementInfo, IonType

from ..constants import HYDROGEN_BINDING_MASS, ModType
from ..diagnostics import CompositionError, PeptacularError
from ..proforma_components import (
    ChargedFormula,
    FormulaElement,
    GlobalChargeCarrier,
    ModificationTags,
    TagMass,
    add_composition,
)
from .mod import Mods
from .positions import to_ion_type
from .utils import H_ELEMENT_INFO, _adjust_mass_value, can_fragment_sequence

if TYPE_CHECKING:
    from .annotation import CHARGE_TYPE, CUSTOM_LOSS_TYPE, ION_TYPE, ISOTOPE_TYPE, ProFormaAnnotation


fe = FormulaElement(element=Element.H, occurance=1)
H_CHARGE_FORMULA = ChargedFormula(formula=(fe,), charge=1)
H_DECHARGE_FORMULA = ChargedFormula(formula=(FormulaElement(element=Element.H, occurance=-1),), charge=-1)

EMPTY_CHARGE_MODS = Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods=None)

# Residue masses are fixed reference data, independent of mutable annotations.
_MONOISOTOPIC_AA_MASSES = {aa: info.monoisotopic_mass for aa, info in AA_LOOKUP.items()}
# Per-residue compositions read once; get_sequence_composition only reads them.
_AA_COMPOSITIONS: dict[str, Counter[ElementInfo] | None] = {aa: info.composition for aa, info in AA_LOOKUP.items()}
_AVERAGE_AA_MASSES = {aa: info.average_mass for aa, info in AA_LOOKUP.items()}


def charge_state(annot: "ProFormaAnnotation") -> int:
    """Body of :attr:`ProFormaAnnotation.charge_state`."""
    charge = annot.charge
    if isinstance(charge, int):
        return charge
    elif isinstance(charge, Mods):
        return sum(mod.get_charge() for mod in charge.mods)
    elif charge is None:
        return 0
    else:
        raise PeptacularError(f"Invalid charge type: {type(charge)}")


def charge_adducts(annot: "ProFormaAnnotation") -> Mods[GlobalChargeCarrier]:
    """Body of :attr:`ProFormaAnnotation.charge_adducts`."""
    charge = annot.charge
    if isinstance(charge, int):
        if charge == 0:
            return EMPTY_CHARGE_MODS
        elif charge > 0:
            s = str(GlobalChargeCarrier(charged_formula=H_CHARGE_FORMULA, occurance=charge))
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={s: 1})
        else:  # charge < 0
            s = str(GlobalChargeCarrier(charged_formula=H_DECHARGE_FORMULA, occurance=-charge))
            return Mods[GlobalChargeCarrier](mod_type=ModType.CHARGE, _mods={s: 1})
    elif isinstance(charge, Mods):
        return charge
    return EMPTY_CHARGE_MODS


def sequence_composition(annot: "ProFormaAnnotation") -> Counter[ElementInfo]:
    """Body of :meth:`ProFormaAnnotation.get_sequence_composition`."""
    sequence_composition: Counter[ElementInfo] = Counter()
    # Count residues once, then scale each residue's composition by its count.
    for aa, n in Counter(annot.stripped_sequence).items():
        residue_comp = _AA_COMPOSITIONS[aa] if aa in _AA_COMPOSITIONS else AA_LOOKUP[aa].composition
        if residue_comp is None:
            raise CompositionError(f"Composition not available for amino acid: {aa}")
        for element, count in residue_comp.items():
            sequence_composition[element] += count * n
    return sequence_composition


def base_comp(annot: "ProFormaAnnotation", skip_labile: bool = False, monoisotopic: bool = True) -> tuple[Counter[ElementInfo], int, float]:
    """Body of ``ProFormaAnnotation._base_comp``: (composition, internal charge, mass-only delta)."""
    total_composition: Counter[ElementInfo] = sequence_composition(annot)
    total_charge = 0  # results from internal formula mods
    total_delta_mass = 0.0  # only from MassTags

    if annot.has_unknown_mods:
        unknown_mods = annot.unknown_mods
        composition, delta_mass, charge = unknown_mods.get_composition_with_delta_mass_charge(monoisotopic=monoisotopic)
        add_composition(total_composition, composition)
        total_delta_mass += delta_mass
        total_charge += charge

    if not skip_labile and annot.has_labile_mods:
        labile_mods = annot.labile_mods
        composition, delta_mass, charge = labile_mods.get_composition_with_delta_mass_charge(monoisotopic=monoisotopic)
        add_composition(total_composition, composition)
        total_delta_mass += delta_mass
        total_charge += charge

    if annot.has_nterm_mods:
        nterm_mods = annot.nterm_mods
        composition, delta_mass, charge = nterm_mods.get_composition_with_delta_mass_charge(monoisotopic=monoisotopic)
        add_composition(total_composition, composition)
        total_delta_mass += delta_mass
        total_charge += charge

    if annot.has_cterm_mods:
        cterm_mods = annot.cterm_mods
        composition, delta_mass, charge = cterm_mods.get_composition_with_delta_mass_charge(monoisotopic=monoisotopic)
        add_composition(total_composition, composition)
        total_delta_mass += delta_mass
        total_charge += charge

    if annot.has_static_mods:
        static_mod_map = annot.map_static_mods_to_indexes()
        for _, mods in static_mod_map.items():
            for mod in mods:
                try:
                    add_composition(total_composition, mod.get_composition())
                except ValueError as e:
                    if isinstance(mod.value, ModificationTags) and isinstance(mod.value.first_tag, TagMass):
                        # MassTag does not have composition, only delta mass
                        total_delta_mass += mod.get_mass(monoisotopic=monoisotopic)
                    else:
                        raise e
                total_charge += mod.get_charge()

    # Internal mods
    if annot.has_internal_mods:
        for mods in annot.internal_mods.values():
            composition, delta_mass, charge = mods.get_composition_with_delta_mass_charge(monoisotopic=monoisotopic)
            add_composition(total_composition, composition)
            total_delta_mass += delta_mass
            total_charge += charge

    # Intervals
    if annot.has_intervals:
        for interval in annot.intervals:
            composition, delta_mass, charge = interval.mods.get_composition_with_delta_mass_charge(monoisotopic=monoisotopic)
            add_composition(total_composition, composition)
            total_delta_mass += delta_mass
            total_charge += charge

    return total_composition, total_charge, total_delta_mass


def comp(
    annot: "ProFormaAnnotation",
    charge: "CHARGE_TYPE | None" = None,
    *,
    ion_type: "ION_TYPE" = IonType.PRECURSOR,
    isotopes: "ISOTOPE_TYPE | None" = None,
    deltas: "CUSTOM_LOSS_TYPE | None" = None,
) -> Counter[ElementInfo]:
    """Body of :meth:`ProFormaAnnotation.comp`."""
    frag = annot.frag(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=True,
        isotopes=isotopes,
        deltas=deltas,
        calculate_with_composition=True,
        _include_sequence=False,
    )

    if frag.composition is None:
        raise PeptacularError("Fragment composition could not be calculated.")

    return frag.composition


def base_mass(annot: "ProFormaAnnotation", monoisotopic: bool = True, skip_labile: bool = False) -> tuple[float, int]:
    """Body of ``ProFormaAnnotation._base_mass``: (residue + modification mass, internal charge)."""
    total_mass = 0.0
    total_charge = 0  # results from internal formula mods

    # Inline mass lookup to avoid function call overhead
    # Amino acids - hot path, optimize heavily
    aa_lookup = _MONOISOTOPIC_AA_MASSES if monoisotopic else _AVERAGE_AA_MASSES
    for aa in annot.stripped_sequence:
        mass = aa_lookup[aa]
        if mass is None:
            raise PeptacularError(f"Mass not available for amino acid: {aa}")
        total_mass += mass

    # Unknown mods
    if annot.has_unknown_mods:
        m, c = annot.unknown_mods.get_mass_charge(monoisotopic=monoisotopic)
        total_mass += m
        total_charge += c

    # Labile mods
    if not skip_labile and annot.has_labile_mods:
        m, c = annot.labile_mods.get_mass_charge(monoisotopic=monoisotopic)
        total_mass += m
        total_charge += c

    # N-terminal mods
    if annot.has_nterm_mods:
        m, c = annot.nterm_mods.get_mass_charge(monoisotopic=monoisotopic)
        total_mass += m
        total_charge += c

    # Internal mods
    if annot.has_internal_mods:
        for mods in annot.internal_mods.values():
            m, c = mods.get_mass_charge(monoisotopic=monoisotopic)
            total_mass += m
            total_charge += c

    # Interval mods
    if annot.has_intervals:
        for interval in annot.intervals:
            m, c = interval.mods.get_mass_charge(monoisotopic=monoisotopic)
            total_mass += m
            total_charge += c

    # C-terminal mods
    if annot.has_cterm_mods:
        m, c = annot.cterm_mods.get_mass_charge(monoisotopic=monoisotopic)
        total_mass += m
        total_charge += c

    # Static mods
    if annot.has_static_mods:
        static_mod_map = annot.map_static_mods_to_indexes()
        for mods in static_mod_map.values():
            for mod in mods:
                total_mass += mod.get_mass(monoisotopic=monoisotopic)
                total_charge += mod.get_charge()

    return total_mass, total_charge


def residue_mass_vector(annot: "ProFormaAnnotation", monoisotopic: bool = True) -> list[float]:
    """Body of ``ProFormaAnnotation._get_mass_vector``: the neutral mass of each one-residue slice."""
    # slice sequence into single aa slcices
    vec: list[float] = []
    for i in range(len(annot)):
        sub_annot = annot.slice(i, i + 1, inplace=False)
        vec.append(
            sub_annot.mass(
                ion_type=IonType.NEUTRAL,
                charge=None,
                monoisotopic=monoisotopic,
                isotopes=None,
            )
        )
    return vec


def residue_comp_vector(annot: "ProFormaAnnotation") -> list[Counter[ElementInfo]]:
    """Body of ``ProFormaAnnotation._get_comp_vector``: the neutral composition of each one-residue slice."""
    # slice sequence into single aa slcices
    vec: list[Counter[ElementInfo]] = []
    for i in range(len(annot)):
        sub_annot = annot.slice(i, i + 1, inplace=False)
        vec.append(
            sub_annot.comp(
                ion_type=IonType.NEUTRAL,
                charge=None,
                isotopes=None,
            )
        )
    return vec


def mass_and_charge(
    annot: "ProFormaAnnotation",
    ion_type: "ION_TYPE",
    charge: "CHARGE_TYPE | None",
    monoisotopic: bool,
    isotopes: "ISOTOPE_TYPE | None",
    deltas: "CUSTOM_LOSS_TYPE | None",
    calculate_with_composition: bool,
) -> tuple[float, int]:
    """Body of ``ProFormaAnnotation._mass_and_charge``: (mass, total charge) of one ion.

    Avoids fragment allocation and annotation copies for ordinary intact ions.
    """
    effective_charge = annot._charge if charge is None else charge
    if (
        ion_type in (IonType.PRECURSOR, IonType.NEUTRAL)
        and isotopes is None
        and deltas is None
        and not calculate_with_composition
        and not annot.has_isotope_mods
        and (effective_charge is None or type(effective_charge) is int and effective_charge >= 0)
    ):
        ion_type = to_ion_type(ion_type)
        can_fragment_sequence(annot.sequence, ion_type)
        total_mass, internal_charge = base_mass(annot, monoisotopic=monoisotopic)
        external_charge = effective_charge or 0
        total_charge = external_charge + internal_charge
        mass = _adjust_mass_value(
            total_mass,
            # H atoms here, the electrons come off below; the binding term lifts H - e to PROTON_MASS.
            (H_ELEMENT_INFO.get_mass(monoisotopic=monoisotopic) + (HYDROGEN_BINDING_MASS if monoisotopic else 0.0)) * external_charge,
            total_charge,
            ion_type,
            monoisotopic,
        )
        return mass, total_charge

    f = annot.frag(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        calculate_with_composition=calculate_with_composition,
        _include_sequence=False,
    )
    return f.mass, f.charge_state

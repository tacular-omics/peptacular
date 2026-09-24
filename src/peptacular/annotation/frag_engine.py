"""Fragmentation engine behind :class:`~peptacular.annotation.ProFormaAnnotation`.

Every function that takes an annotation takes it as its first argument. The annotation's
``frag()``, ``fragment()``, ``fast_fragment()`` and their private helpers are short calls
into this module; use those methods, not these functions.
"""

from __future__ import annotations

from collections.abc import Callable, Generator, Sequence
from functools import partial
from itertools import product
from typing import TYPE_CHECKING, Any

from tacular import (
    FRAGMENT_ION_LOOKUP,
    NEUTRAL_DELTA_LOOKUP,
    FragmentIonInfo,
    IonType,
    NeutralDeltaInfo,
)

from ..constants import ELECTRON_MASS, PROTON_MASS
from ..diagnostics import (
    CompositionError,
    InvalidAdjustmentError,
    PeptacularError,
    UnsupportedOperationError,
)
from ..proforma_components import ChargedFormula
from .cached_comps import DeltaInfo, IsotopeInfo
from .frag import _ION_TYPE_TO_MZPAF_SERIES, proton_binding_offset
from .mass import _AVERAGE_AA_MASSES, _MONOISOTOPIC_AA_MASSES
from .mass import base_comp as _base_comp
from .mass import base_mass as _base_mass
from .positions import to_ion_type, validate_position
from .utils import (
    H_ELEMENT_INFO,
    SATELLITE_TRIM_END,
    SATELLITE_TRIM_START,
    Fragment,
    _adjust_mass_value,
    _ion_mass,
    adjust_comp,
    adjust_mass_mz,
    can_fragment_sequence,
    validate_mass,
)

if TYPE_CHECKING:
    from .annotation import CHARGE_TYPE, CUSTOM_LOSS_TYPE, ION_TYPE, ISOTOPE_TYPE, LOSS_TYPE, ProFormaAnnotation

__all__ = [
    "get_loss_combinations",
    "build_mass_vector",
    "frag_one",
    "satellite_mod_error",
    "frag_impl",
    "frag",
    "series_mass_vector",
    "fragment_series",
    "fragment_ions",
    "default_fragment_charges",
    "fragment",
    "fast_fragment",
]


def _as_options(value: Any) -> Any:
    """Wrap a single str/int option (an ion type, charge, isotope or loss name) in a tuple."""
    if isinstance(value, str | int):
        return (value,)
    return value


def _carrier_mass(monoisotopic: bool) -> float:
    """Mass one default (protonated) charge adds: CODATA ``PROTON_MASS``, or average H minus an electron."""
    if monoisotopic:
        return PROTON_MASS
    return H_ELEMENT_INFO.get_mass(monoisotopic=False) - ELECTRON_MASS


def _unless_impossible_loss(ndelta: DeltaInfo, make: Callable[..., Fragment], **kwargs: Any) -> Fragment | None:
    """Build one ion for ``fragment()``, or None when that ion cannot exist.

    ``neutral_deltas`` offers a loss wherever one of its residues occurs (H3PO4 on any S/T), so
    a loss can ask for more of an element than the fragment has (no phosphorus on an
    unmodified S). An ion can also lack the atoms its own offset removes (the one-residue a1 of
    ``G-[Amidated]`` at charge -1), or the atoms the caller's ``isotopes`` swap (``{"15N": 3}``
    on b1). Such ions are skipped; ``fragment()`` raises when the caller's isotopes leave no
    ion at all. The error still propagates when the caller's own ``deltas`` are at fault: the
    ion exists without them. An explicit ``frag()`` does not come through here and always raises.
    """
    try:
        return make(**kwargs)
    except InvalidAdjustmentError as error:
        if ndelta._items:
            return None
        delta: DeltaInfo = kwargs["delta"]
        isotope: IsotopeInfo = kwargs["isotope"]
        if not delta.deltas and not isotope.data:
            return None
        no_isotope = IsotopeInfo.from_input(None)
        try:
            make(**{**kwargs, "delta": DeltaInfo.from_input(None), "isotope": no_isotope})
        except InvalidAdjustmentError:
            return None  # the ion itself cannot exist
        if delta.deltas:
            try:
                make(**{**kwargs, "isotope": no_isotope})
            except InvalidAdjustmentError:
                raise error from None  # the caller's deltas are at fault
        return None  # the caller's isotopes do not fit this ion


def get_loss_combinations(losses: dict[NeutralDeltaInfo, int], max_losses: int) -> list[DeltaInfo]:
    """Generate all combinations of losses up to max_losses."""
    if not losses:
        return [DeltaInfo.from_input(None)]

    # Generate all possible count combinations for each loss
    loss_items = list(losses.items())
    count_ranges = [range(count + 1) for _, count in loss_items]

    loss_combinations: list[dict[NeutralDeltaInfo, int] | None] = [None]

    for counts in product(*count_ranges):
        total_losses = sum(counts)
        # Skip if no losses or exceeds max
        if total_losses == 0 or total_losses > max_losses:
            continue

        # Build combination dict
        combo = {}
        for (loss, _), count in zip(loss_items, counts, strict=True):
            if count > 0:
                combo[loss] = count

        loss_combinations.append(combo)

    # Convert to DeltaInfo
    delta_combinations: list[DeltaInfo] = []
    for loss_combo in loss_combinations:
        if loss_combo is None:
            delta_combinations.append(DeltaInfo.from_input(None))
        else:
            formula_dict: dict[ChargedFormula, int] = {}
            for nd, count in loss_combo.items():
                nd_formula: ChargedFormula = ChargedFormula.from_composition(nd.composition)
                formula_dict[nd_formula] = count
            delta_combinations.append(DeltaInfo.from_input(formula_dict))  # type: ignore

    return delta_combinations


def build_mass_vector(annot: ProFormaAnnotation, monoisotopic: bool = True) -> list[float]:
    """Build a per-residue mass array without sequence slicing.

    Each element is the residue mass plus any modifications localised to that
    position (internal mods, N-/C-terminal mods at the respective ends, and
    static mods mapped to their target residues).

    :param monoisotopic: Use monoisotopic masses when ``True``, average masses when ``False``.
    :type monoisotopic: bool
    :return: List of per-residue masses, length == len(annot).
    :rtype: list[float]
    :raises PeptacularError: If the annotation contains unknown mods or interval mods.
    """
    if annot.has_unknown_mods or annot.has_intervals:
        raise UnsupportedOperationError(f"fast_fragment not supported for sequences with unknown modifications or intervals: {str(annot)}")

    aa_lookup = _MONOISOTOPIC_AA_MASSES if monoisotopic else _AVERAGE_AA_MASSES
    masses: list[float] = []
    for aa in annot.stripped_sequence:
        m = aa_lookup[aa]
        if m is None:
            raise PeptacularError(f"Mass not available for amino acid: {aa}")
        masses.append(m)

    if annot.has_nterm_mods:
        m, _ = annot.nterm_mods.get_mass_charge(monoisotopic=monoisotopic)
        masses[0] += m

    if annot.has_internal_mods:
        for pos, mods in annot.internal_mods.items():
            m, _ = mods.get_mass_charge(monoisotopic=monoisotopic)
            masses[pos] += m

    if annot.has_cterm_mods:
        m, _ = annot.cterm_mods.get_mass_charge(monoisotopic=monoisotopic)
        masses[-1] += m

    if annot.has_static_mods:
        static_mod_map = annot.map_static_mods_to_indexes()
        for pos, mods_list in static_mod_map.items():
            # map_static_mods_to_indexes uses -1 for N-term and -2 for C-term
            if pos == -1:
                pos = 0
            elif pos == -2:
                pos = len(masses) - 1
            for mod in mods_list:
                masses[pos] += mod.get_mass(monoisotopic=monoisotopic)

    return masses


def frag_one(
    annot: ProFormaAnnotation,
    ion_type: IonType,
    monoisotopic: bool,
    isotope: IsotopeInfo,
    delta: DeltaInfo,
    calculate_with_composition: bool,
    parent_sequence: str,
    parent_sequence_length: int,
    position: int | tuple[int, int] | None,
) -> Fragment:
    # Satellite ions: the residue whose side chain is cleaved is not in the residue sum;
    # the ion offset carries its remnant (mzPAF 1.0.1). Its modifications leave with it,
    # but a terminal modification sits on the backbone and stays: a full-length d ion
    # keeps the C-terminal mod, a full-length v/w ion keeps the N-terminal mod.
    blocked = satellite_mod_error(annot, ion_type)
    if blocked is not None:
        raise PeptacularError(blocked)
    ion_annot = annot
    if ion_type in SATELLITE_TRIM_END:
        ion_annot = annot.slice(0, len(annot) - 1, inplace=False)
        if annot.has_cterm_mods:
            ion_annot.set_cterm_mods(annot.cterm_mods, validate=False)
    elif ion_type in SATELLITE_TRIM_START:
        ion_annot = annot.slice(1, len(annot), inplace=False)
        if annot.has_nterm_mods:
            ion_annot.set_nterm_mods(annot.nterm_mods, validate=False)
    return frag_impl(
        ion_annot,
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotope=isotope,
        delta=delta,
        calculate_with_composition=calculate_with_composition,
        parent_sequence=parent_sequence,
        parent_sequence_length=parent_sequence_length,
        position=position,
    )


def satellite_mod_error(annot: ProFormaAnnotation, ion_type: IonType) -> str | None:
    """Why a d or w ion of this (sub)sequence is undefined, or None when it is defined.

    A d or w ion keeps part of the cleaved residue's side chain (its beta substituent), so
    it is not defined when that residue carries a modification, explicit or from a global
    fixed modification (as in paftacular). A v ion loses the whole side chain, and its
    modification with it, so v ions are always defined.
    """
    if ion_type == IonType.V or not annot:
        return None
    if ion_type in SATELLITE_TRIM_END:
        index = len(annot) - 1
    elif ion_type in SATELLITE_TRIM_START:
        index = 0
    else:
        return None
    if annot.has_internal_mods_at_index(index) or (annot.has_static_mods and index in annot.map_static_mods_to_indexes()):
        label = _ION_TYPE_TO_MZPAF_SERIES.get(ion_type, ion_type.value)
        return f"{label} ion is not defined when residue {annot.stripped_sequence[index]} carries a modification"
    return None


def frag_impl(
    annot: ProFormaAnnotation,
    ion_type: IonType,
    monoisotopic: bool,
    isotope: IsotopeInfo,
    delta: DeltaInfo,
    calculate_with_composition: bool,
    parent_sequence: str,
    parent_sequence_length: int,
    position: int | tuple[int, int] | None,
) -> Fragment:
    # Dont include labile mods for fragment ions
    skip_labile = True
    if ion_type == IonType.NEUTRAL or ion_type == IonType.PRECURSOR:
        skip_labile = False

    # Elemental adjustments share one order in both calculation modes.
    # Mass-only modifications remain additive and do not invent atom counts.
    formula_deltas: dict[ChargedFormula | float, int] = {key: count for key, count in delta.deltas.items() if isinstance(key, ChargedFormula)}
    charge_carriers = annot.charge_adducts
    removes_atoms = any(count < 0 for mod in charge_carriers for count in mod.get_composition().values())
    if annot.has_isotope_mods or calculate_with_composition or isotope.data or formula_deltas or removes_atoms:
        base_comp, base_charge, delta_mass = _base_comp(annot, skip_labile=skip_labile, monoisotopic=monoisotopic)
        if calculate_with_composition and (delta_mass != 0.0 or delta.has_floats):
            raise CompositionError("Cannot calculate composition with delta mass changes. Use mass() or mz() instead.")
        result = adjust_comp(
            base_comp=base_comp,
            charge=charge_carriers,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            delta=DeltaInfo(formula_deltas),
            inplace=True,
            isotope_map=annot.map_isotopes() if annot.has_isotope_mods else None,
            position=position,
            parent_sequence=parent_sequence,
            parent_sequence_length=parent_sequence_length,
            internal_charge=base_charge,
            isotope_as_mass=not calculate_with_composition,
        )
        if not calculate_with_composition and not annot.has_isotope_mods:
            # The composition above only validates the ion (atoms left for a formula loss,
            # an isotope swap or a deprotonation). The mass comes from the listed masses,
            # as for the plain ion, so a loss or isotope peak is exactly the plain ion plus
            # its delta. Summing the mods' compositions instead would move named mods off
            # their listed mass (Oxidation 15.994915 vs 15.9949146 from O).
            base_mass, _ = _base_mass(annot, monoisotopic=monoisotopic, skip_labile=skip_labile)
            mass = _adjust_mass_value(
                base_mass,
                charge_carriers.get_mass(monoisotopic=monoisotopic) + proton_binding_offset(charge_carriers, monoisotopic),
                result.charge_state,
                ion_type,
                monoisotopic,
                isotope.get_mass_delta(monoisotopic),
                delta.get_mass_delta(monoisotopic),
            )
            result = result._replace(mass=mass, _composition=None, _deltas=delta.to_fragment_mapping)
        elif not calculate_with_composition:
            # A global isotope label (<13C>) changes every atom's mass, so the labelled
            # composition is the mass; mass-only tags and float deltas are added on top.
            mass = result.mass + delta_mass + sum(key * count for key, count in delta.deltas.items() if isinstance(key, float))
            result = result._replace(mass=mass, _composition=None, _deltas=delta.to_fragment_mapping)
        else:
            result = result._replace(_deltas=delta.to_fragment_mapping)
        validate_mass(result.mass)
        return result

    base_mass, base_charge = _base_mass(annot, monoisotopic=monoisotopic, skip_labile=skip_labile)

    return adjust_mass_mz(
        base=base_mass,
        charge=charge_carriers,
        monoisotopic=monoisotopic,
        ion_type=ion_type,
        isotope=isotope,
        delta=delta,
        position=position,
        parent_sequence=parent_sequence,
        parent_sequence_length=parent_sequence_length,
        internal_charge=base_charge,
    )


def frag(
    annot: ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    deltas: CUSTOM_LOSS_TYPE | None = None,
    calculate_with_composition: bool = False,
    position: int | tuple[int, int] | None = None,
    _include_sequence: bool = True,
) -> Fragment:
    """Body of :meth:`ProFormaAnnotation.frag`: one ion, preferring ``charge`` over the annotation's charge."""
    ion_type = to_ion_type(ion_type)
    delta_info = DeltaInfo.from_input(deltas)

    inplace = False

    frag_annot = annot
    if charge is not None:  # update charge
        frag_annot = frag_annot.set_charge(charge, inplace=inplace)
        inplace = True

    parent_sequence = None
    if _include_sequence:
        parent_sequence = frag_annot.serialize(exclude_charge=False)
    else:
        parent_sequence = ""

    if position is None:
        ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]
        if ion_info.is_forward | ion_info.is_backward:
            position = len(annot)  # default to full length for terminal ions
        if ion_info.is_internal:
            if ion_info.ion_type == IonType.IMMONIUM and len(frag_annot) != 1:
                raise PeptacularError("Immonium ions must be single amino acids, or the position must be specified.")
            position = (
                1,
                len(annot),
            )  # use whole sequence for internal ions by default
        if ion_info.is_intact:
            position = None  # use whole sequence for intact ions

    _pos: tuple[int, int] | None = validate_position(ion_type, position, len(annot))

    iso_info = IsotopeInfo.from_input(isotopes)

    # get the appropriate fragment annotation based on the position parameter
    match _pos:
        case None:
            pass
        case tuple() as pios_tuple:
            pos_start, pos_end = pios_tuple
            frag_annot = frag_annot.slice(pos_start, pos_end, inplace=inplace)
        case _:
            raise PeptacularError(f"Invalid position type: {type(position)}")

    # Checked on the fragment itself: satellite ions depend on its terminal residue.
    ion_type = can_fragment_sequence(frag_annot.sequence, ion_type)

    return frag_one(
        frag_annot,
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotope=iso_info,
        delta=delta_info,
        calculate_with_composition=calculate_with_composition,
        parent_sequence=parent_sequence,
        parent_sequence_length=len(annot),
        position=position,
    )


def series_mass_vector(annot: ProFormaAnnotation, monoisotopic: bool, calculate_with_composition: bool) -> list[float] | None:
    """Per-residue masses for the terminal-series fast path, or None when it does not apply.

    Terminal mods sit on the first/last residue, so a prefix sum of length ``i`` equals the
    mass of ``annot.slice(0, i)`` (the C-terminal mods only join at ``i == len``) and a suffix
    sum equals ``annot[len - i:]``. Anything that needs the composition path, or that slicing
    treats specially, returns None so the caller slices instead.
    """
    if (
        calculate_with_composition
        or annot.has_isotope_mods
        or annot.has_static_mods
        or annot.has_unknown_mods
        or annot.has_intervals
        or any(count < 0 for mod in annot.charge_adducts for count in mod.get_composition().values())
    ):
        return None
    aa_lookup = _MONOISOTOPIC_AA_MASSES if monoisotopic else _AVERAGE_AA_MASSES
    masses: list[float] = []
    for aa in annot.stripped_sequence:
        m = aa_lookup[aa]
        if m is None:
            return None
        masses.append(m)
    if annot.has_nterm_mods:
        m, c = annot.nterm_mods.get_mass_charge(monoisotopic=monoisotopic)
        if c:
            return None
        masses[0] += m
    if annot.has_internal_mods:
        for pos, mods in annot.internal_mods.items():
            m, c = mods.get_mass_charge(monoisotopic=monoisotopic)
            if c:
                return None
            masses[pos] += m
    if annot.has_cterm_mods:
        m, c = annot.cterm_mods.get_mass_charge(monoisotopic=monoisotopic)
        if c:
            return None
        masses[-1] += m
    return masses


def fragment_series(
    annot: ProFormaAnnotation,
    ion_type: IonType,
    *,
    forward: bool,
    monoisotopic: bool,
    isotopes: list[IsotopeInfo],
    deltas: list[DeltaInfo],
    neutral_deltas: list[NeutralDeltaInfo],
    calculate_with_composition: bool,
    parent_sequence: str,
    parent_sequence_length: int,
    max_deltas: int,
    min_length: int | None,
    max_length: int | None,
    _fast: bool = True,
) -> Generator[Fragment, None, None]:
    """Yield one terminal ion series (forward: b1..bn, backward: y1..yn).

    Ions with plain mass adjustments are built from prefix sums without slicing the
    annotation; everything else (composition mode, formula deltas, isotope swaps,
    satellite ions, ...) slices and goes through :meth:`_frag`.
    """
    n = len(annot)
    stripped = annot.stripped_sequence
    masses = series_mass_vector(annot, monoisotopic, calculate_with_composition) if _fast else None

    cumulative: list[float] = []
    charge_mass = 0.0
    external_charge = 0
    adducts: tuple[str, ...] | None = None
    if masses is not None:
        total = 0.0
        cumulative.append(total)
        for m in masses if forward else reversed(masses):
            total += m
            cumulative.append(total)
        charge_carriers = annot.charge_adducts
        charge_mass = charge_carriers.get_mass(monoisotopic=monoisotopic) + proton_binding_offset(charge_carriers, monoisotopic)
        external_charge = charge_carriers.get_charge()
        if not all(m.value.is_protonated for m in charge_carriers.mods):
            adducts = tuple(key for key, count in charge_carriers._mods.items() for _ in range(count)) if charge_carriers._mods else None

    # Per-ion work that does not depend on the position is done once per series:
    # the (isotope, delta, loss) products are cached per loss-site count, and the
    # ion-type lookup per (possibly residue-specific) ion type.
    combo_cache: dict[tuple[tuple[NeutralDeltaInfo, int], ...], list[tuple[IsotopeInfo, DeltaInfo, DeltaInfo, bool, float, float]]] = {}
    ion_cache: dict[IonType, tuple[IonType, bool, float]] = {}
    loss_dict: dict[NeutralDeltaInfo, int] = {}
    for i in range(1, n + 1):
        if min_length is not None and i < min_length:
            continue
        if max_length is not None and i > max_length:
            break

        sub_sequence = stripped[:i] if forward else stripped[n - i :]
        try:
            frag_type = can_fragment_sequence(sub_sequence, ion_type)
        except ValueError:
            continue

        if neutral_deltas:
            loss_dict.clear()
            for nd in neutral_deltas:
                loss_dict[nd] = min(nd.calculate_loss_sites(sub_sequence), max_deltas)
        loss_key = tuple(loss_dict.items())
        products = combo_cache.get(loss_key)
        if products is None:
            products = []
            for isotope in isotopes:
                for delta in deltas:
                    for ndelta in get_loss_combinations(loss_dict, max_deltas):
                        combined_delta = delta + ndelta
                        plain = not isotope.data and not any(isinstance(k, ChargedFormula) for k in combined_delta.deltas)
                        iso_mass = isotope.get_mass_delta(monoisotopic) if plain else 0.0
                        delta_mass = combined_delta.get_mass_delta(monoisotopic) if plain else 0.0
                        products.append((isotope, combined_delta, ndelta, plain, iso_mass, delta_mass))
            combo_cache[loss_key] = products

        ion_entry = ion_cache.get(frag_type)
        if ion_entry is None:
            fast_type = masses is not None and frag_type not in SATELLITE_TRIM_END and frag_type not in SATELLITE_TRIM_START
            frag_ion_type = FRAGMENT_ION_LOOKUP[frag_type].ion_type
            ion_entry = (frag_ion_type, fast_type, _ion_mass(frag_ion_type, monoisotopic) if fast_type else 0.0)
            ion_cache[frag_type] = ion_entry
        frag_ion_type, fast_type, ion_mass = ion_entry
        sub_annot: ProFormaAnnotation | None = None

        for isotope, combined_delta, ndelta, plain, iso_mass, delta_mass in products:
            if fast_type and plain:
                # Same arithmetic order as adjust_mass_mz / _adjust_mass_value.
                mass = cumulative[i]
                mass += iso_mass
                mass += delta_mass
                mass += charge_mass
                mass += ion_mass
                mass -= external_charge * ELECTRON_MASS
                validate_mass(mass)
                yield Fragment(
                    ion_type=frag_ion_type,
                    position=i,
                    mass=mass,
                    monoisotopic=monoisotopic,
                    charge_state=external_charge,
                    charge_adducts=adducts,
                    external_charge=external_charge,
                    isotopes=isotope.to_fragment_mapping,
                    deltas=combined_delta.to_fragment_mapping,
                    composition=None,
                    parent_sequence=parent_sequence,
                    parent_sequence_length=parent_sequence_length,
                )
                continue
            if sub_annot is None:
                sub_annot = annot.slice(0, i, inplace=False) if forward else annot[n - i : n]
            if frag_type in SATELLITE_TRIM_END or frag_type in SATELLITE_TRIM_START:
                # A series skips d/w ions of a modified cleaved residue; an explicit
                # frag() of the same ion raises.
                if satellite_mod_error(sub_annot, frag_type) is not None:
                    break
            fragment = _unless_impossible_loss(
                ndelta,
                partial(frag_one, sub_annot),
                ion_type=frag_type,
                monoisotopic=monoisotopic,
                isotope=isotope,
                delta=combined_delta,
                calculate_with_composition=calculate_with_composition,
                parent_sequence=parent_sequence,
                parent_sequence_length=parent_sequence_length,
                position=i,
            )
            if fragment is not None:
                yield fragment


def fragment_ions(
    annot: ProFormaAnnotation,
    ion_type: IonType,
    monoisotopic: bool = True,
    *,
    isotopes: list[IsotopeInfo],
    deltas: list[DeltaInfo],
    neutral_deltas: list[NeutralDeltaInfo],
    calculate_with_composition: bool,
    parent_sequence: str,
    parent_sequence_length: int,
    max_deltas: int,
    min_length: int | None,
    max_length: int | None,
    _expand: bool = True,
) -> Generator[Fragment, None, None]:
    if annot.has_unknown_mods or annot.has_intervals:
        raise PeptacularError(f"Fragmentation not supported for sequences with unknown modifications or intervals: {str(annot)}")

    # "d" and "w" cover the generic ion and the residue-specific a/b variants
    # (d-valine, da-/db-threonine, ...); each variant only forms on its own residues.
    if _expand and ion_type in (IonType.D, IonType.W):
        family = (IonType.D, IonType.DA, IonType.DB) if ion_type == IonType.D else (IonType.W, IonType.WA, IonType.WB)
        for member in family:
            yield from fragment_ions(
                annot,
                member,
                monoisotopic,
                isotopes=isotopes,
                deltas=deltas,
                neutral_deltas=neutral_deltas,
                calculate_with_composition=calculate_with_composition,
                parent_sequence=parent_sequence,
                parent_sequence_length=parent_sequence_length,
                max_deltas=max_deltas,
                min_length=min_length,
                max_length=max_length,
                _expand=False,
            )
        return

    ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]

    loss_dict: dict[NeutralDeltaInfo, int] = {}
    # Terminal series: forward ions (b1, b2, ...) grow from the N-terminus,
    # backward ions (y1, y2, ...) from the C-terminus.
    if ion_info.is_forward or ion_info.is_backward:
        yield from fragment_series(
            annot,
            ion_type,
            forward=ion_info.is_forward,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            calculate_with_composition=calculate_with_composition,
            parent_sequence=parent_sequence,
            parent_sequence_length=parent_sequence_length,
            max_deltas=max_deltas,
            min_length=min_length,
            max_length=max_length,
        )

    elif ion_info.is_intact:
        if min_length is not None and len(annot) < min_length:
            return

        if max_length is not None and len(annot) > max_length:
            return

        if neutral_deltas:
            loss_dict.clear()
            for nd in neutral_deltas:
                loss_dict[nd] = min(nd.calculate_loss_sites(annot.sequence), max_deltas)

        neutral_delta_combinations = get_loss_combinations(loss_dict, max_deltas)

        for isotope in isotopes:
            for delta in deltas:
                for ndelta in neutral_delta_combinations:
                    combined_delta = delta + ndelta
                    fragment = _unless_impossible_loss(
                        ndelta,
                        partial(frag_one, annot),
                        ion_type=ion_type,
                        monoisotopic=monoisotopic,
                        isotope=isotope,
                        delta=combined_delta,
                        calculate_with_composition=calculate_with_composition,
                        parent_sequence=parent_sequence,
                        parent_sequence_length=parent_sequence_length,
                        # Intact precursor/neutral ions represent the whole sequence.
                        # Keep position unset so Fragment.composition/sequence do not
                        # try to validate an integer cleavage position for a non-series ion.
                        position=None,
                    )
                    if fragment is not None:
                        yield fragment
    elif ion_info.is_internal:
        if ion_info.ion_type == IonType.IMMONIUM:
            # Immonium ions are single residue fragments
            for i in range(1, len(annot) + 1):
                sub_annot = annot.slice(i - 1, i, inplace=False)

                if neutral_deltas:
                    loss_dict.clear()
                    for nd in neutral_deltas:
                        loss_dict[nd] = min(nd.calculate_loss_sites(sub_annot.sequence), max_deltas)

                neutral_delta_combinations = get_loss_combinations(loss_dict, max_deltas)

                for isotope in isotopes:
                    for delta in deltas:
                        for ndelta in neutral_delta_combinations:
                            combined_delta = delta + ndelta
                            fragment = _unless_impossible_loss(
                                ndelta,
                                partial(frag_one, sub_annot),
                                ion_type=ion_type,
                                monoisotopic=monoisotopic,
                                isotope=isotope,
                                delta=combined_delta,
                                calculate_with_composition=calculate_with_composition,
                                parent_sequence=parent_sequence,
                                parent_sequence_length=parent_sequence_length,
                                position=i,  # Position is the residue index
                            )
                            if fragment is not None:
                                yield fragment
        else:
            # gen all internal fragmetns from 1 to n-1
            for start in range(2, len(annot)):  # Start from position 1 to len-1
                for end in range(start, len(annot)):  # End before C-terminus
                    # Apply length filters
                    if min_length is not None and (end - start + 1) < min_length:
                        continue
                    if max_length is not None and (end - start + 1) > max_length:
                        continue
                    sub_annot = annot.slice(start - 1, end, inplace=False)

                    if neutral_deltas:
                        loss_dict.clear()
                        for nd in neutral_deltas:
                            loss_dict[nd] = min(nd.calculate_loss_sites(sub_annot.sequence), max_deltas)

                    neutral_delta_combinations = get_loss_combinations(loss_dict, max_deltas)

                    for isotope in isotopes:
                        for delta in deltas:
                            for ndelta in neutral_delta_combinations:
                                combined_delta = delta + ndelta
                                fragment = _unless_impossible_loss(
                                    ndelta,
                                    partial(frag_one, sub_annot),
                                    ion_type=ion_type,
                                    monoisotopic=monoisotopic,
                                    isotope=isotope,
                                    delta=combined_delta,
                                    calculate_with_composition=calculate_with_composition,
                                    parent_sequence=parent_sequence,
                                    parent_sequence_length=parent_sequence_length,
                                    position=(
                                        start,
                                        end,
                                    ),
                                )
                                if fragment is not None:
                                    yield fragment


def default_fragment_charges(charge_state: int) -> tuple[int, ...]:
    """Return default fragment charge states derived from a precursor charge state.

    For a positive precursor charge ``c``, returns ``1, 2, …, c-1``.
    For a negative precursor charge ``c``, returns ``-1, -2, …, c+1``.
    Falls back to ``(1,)`` when ``charge_state`` is 0 (unannotated) or ±1.
    """
    if charge_state > 1:
        return tuple(range(1, charge_state))
    if charge_state < -1:
        return tuple(range(-1, charge_state, -1))
    return (1,)  # charge 0 (unannotated) or ±1 — last resort


def fragment(
    annot: ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: CHARGE_TYPE | Sequence[CHARGE_TYPE] | None = None,
    *,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (),
    calculate_with_composition: bool = False,
    max_ndeltas: int = 1,
    min_length: int | None = None,
    max_length: int | None = None,
) -> list[Fragment]:
    """Body of :meth:`ProFormaAnnotation.fragment`: every ion of each ion type and charge."""
    ion_types = _as_options(ion_types)
    isotopes = _as_options(isotopes)
    neutral_deltas = _as_options(neutral_deltas)
    if charges is None:
        charges = default_fragment_charges(annot.charge_state)
    charges = _as_options(charges)

    # charge_infos: list[ChargeCarrierInfo] = [ChargeCarrierInfo.from_input(charge) for charge in charges]

    isotope_infos: list[IsotopeInfo] = [IsotopeInfo.from_input(isotope) for isotope in isotopes]

    neutral_deltas_infos: list[NeutralDeltaInfo] = []
    if neutral_deltas_infos is not None:
        for loss in neutral_deltas:
            if loss is None:
                continue
            if isinstance(loss, NeutralDeltaInfo):
                neutral_deltas_infos.append(loss)
            else:
                nd: NeutralDeltaInfo = NEUTRAL_DELTA_LOOKUP[loss]
                neutral_deltas_infos.append(nd)

    delta_infos = [DeltaInfo.from_input(loss) for loss in deltas]

    fragments: list[Fragment] = []
    for charge in charges:
        charged_annot = annot.set_charge(charge, inplace=False)
        sequence = charged_annot.serialize()
        for ion in ion_types:
            fragments.extend(
                list(
                    fragment_ions(
                        charged_annot,
                        ion_type=to_ion_type(ion),
                        monoisotopic=monoisotopic,
                        isotopes=isotope_infos,
                        deltas=delta_infos,
                        neutral_deltas=neutral_deltas_infos,
                        calculate_with_composition=calculate_with_composition,
                        parent_sequence=sequence,
                        parent_sequence_length=len(charged_annot),
                        max_deltas=max_ndeltas,
                        min_length=min_length,
                        max_length=max_length,
                    )
                )
            )
    if not fragments and any(info.data for info in isotope_infos):
        # Each ion the caller's isotopes do not fit is skipped; when none is left, say so.
        plain = fragment(
            annot,
            ion_types,
            charges,
            monoisotopic=monoisotopic,
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            max_ndeltas=max_ndeltas,
            calculate_with_composition=calculate_with_composition,
            min_length=min_length,
            max_length=max_length,
        )
        if plain:
            raise InvalidAdjustmentError(f"isotopes={isotopes!r} do not fit any requested ion: every ion has fewer atoms of the swapped element")
    return fragments


def fast_fragment(
    annot: ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: int | Sequence[int] | None = None,
    *,
    monoisotopic: bool = True,
) -> dict[tuple[IonType, int], list[float]]:
    """Body of :meth:`ProFormaAnnotation.fast_fragment`: prefix/suffix sums of the residue masses."""
    ion_types = _as_options(ion_types)
    if charges is None:
        charges = default_fragment_charges(annot.charge_state)
    charges = _as_options(charges)
    for charge in charges:
        if isinstance(charge, bool) or not isinstance(charge, int) or charge == 0:
            raise PeptacularError("fast_fragment charges must be nonzero integers")
    supported = {IonType.A, IonType.B, IonType.C, IonType.X, IonType.Y, IonType.Z, IonType.PRECURSOR, IonType.NEUTRAL}
    for ion_type_input in ion_types:
        if to_ion_type(ion_type_input) not in supported:
            raise UnsupportedOperationError(f"Ion type {ion_type_input!r} is not supported in fast_fragment(). Use fragment() instead.")

    n = len(annot)
    mass_vec = build_mass_vector(annot, monoisotopic=monoisotopic)
    mod_groups = [annot.nterm_mods, annot.cterm_mods, *annot.internal_mods.values()]
    intrinsic_charge = any(mod.get_charge() for mods in mod_groups for mod in mods)
    intrinsic_charge |= any(mod.get_charge() for mods in annot.map_static_mods_to_indexes().values() for mod in mods)
    needs_fallback = annot.has_isotope_mods or annot.has_labile_mods or intrinsic_charge or not monoisotopic or any(c < 0 for c in charges)
    if needs_fallback:
        fallback: dict[tuple[IonType, int], list[float]] = {}
        for charge in charges:
            for ion_type_input in ion_types:
                ion_type = to_ion_type(ion_type_input)
                ion_info = FRAGMENT_ION_LOOKUP[ion_type]
                if ion_info.is_intact:
                    value = annot.frag(ion_type=ion_type, charge=charge, monoisotopic=monoisotopic).mz
                    fallback[(ion_type, charge)] = [value] * n
                else:
                    fallback[(ion_type, charge)] = [
                        annot.frag(ion_type=ion_type, charge=charge, monoisotopic=monoisotopic, position=position).mz for position in range(1, n + 1)
                    ]
        return fallback
    result: dict[tuple[IonType, int], list[float]] = {}

    # A proton charge carrier weighs PROTON_MASS (monoisotopic), as in fragment().
    proton_offset = _carrier_mass(monoisotopic)
    for charge in charges:
        charge_offset = charge * proton_offset
        for ion_type_input in ion_types:
            ion_type = to_ion_type(ion_type_input)
            ion_info: FragmentIonInfo = FRAGMENT_ION_LOOKUP[ion_type]
            ion_offset = _ion_mass(ion_type, monoisotopic)

            abs_charge = abs(charge)
            if ion_info.is_forward:
                prefix = 0.0
                masses_out: list[float] = []
                for m in mass_vec:
                    prefix += m
                    masses_out.append((prefix + ion_offset + charge_offset) / abs_charge)
                result[(ion_type, charge)] = masses_out

            elif ion_info.is_backward:
                prefix = 0.0
                masses_out = []
                for i in range(n - 1, -1, -1):
                    prefix += mass_vec[i]
                    masses_out.append((prefix + ion_offset + charge_offset) / abs_charge)
                result[(ion_type, charge)] = masses_out

            elif ion_info.is_intact:
                total = sum(mass_vec)
                mz = (total + ion_offset + charge_offset) / abs_charge
                result[(ion_type, charge)] = [mz] * n

            elif ion_info.is_internal:
                result[(ion_type, charge)] = [(mass_vec[i] + ion_offset + charge_offset) / abs_charge for i in range(n)]

    for values in result.values():
        for value in values:
            validate_mass(value)
    return result

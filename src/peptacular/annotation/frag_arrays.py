"""Fragment ions as numpy columns: one array per field, one row per ion.

numpy is optional (``pip install "peptacular[numpy]"``) and imported only when
:func:`fragment_arrays` runs. The values are those of :meth:`ProFormaAnnotation.fragment`, in
the same order. Plain b/y-style series (a, b, c, x, y, z with no isotope swap, formula delta or
neutral loss) are computed with numpy prefix sums in the same arithmetic order as
``fragment()``; everything else is built through ``fragment()`` and copied into the columns.
"""

from collections.abc import Mapping, Sequence
from functools import lru_cache
from typing import TYPE_CHECKING, Any

from tacular import FRAGMENT_ION_LOOKUP, IonType

from ..constants import ELECTRON_MASS
from ..interop._errors import MissingOptionalDependencyError
from ..proforma_components import ChargedFormula
from .cached_comps import DeltaInfo, IsotopeInfo
from .frag import Fragment, _format_counts, proton_binding_offset
from .positions import to_ion_type
from .utils import FRAGMENT_RULES, SATELLITE_TRIM_END, SATELLITE_TRIM_START, _ion_mass, validate_mass

if TYPE_CHECKING:
    import numpy as np

    from .annotation import ProFormaAnnotation

__all__ = ["FRAGMENT_ARRAY_KEYS", "fragment_arrays"]

FRAGMENT_ARRAY_KEYS: tuple[str, ...] = (
    "peptide_index",
    "ion_type",
    "position",
    "end_position",
    "charge_state",
    "mz",
    "mass",
    "isotope",
    "isotope_label",
    "delta_label",
    "delta_mass",
)
"""Keys of the dict :func:`fragment_arrays` returns, in order."""

_INT_KEYS = ("peptide_index", "position", "end_position", "charge_state", "isotope")
_FLOAT_KEYS = ("mz", "mass", "delta_mass")
_STR_KEYS = ("ion_type", "isotope_label", "delta_label")


def _require_numpy() -> Any:
    try:
        import numpy
    except ModuleNotFoundError as exc:
        if exc.name != "numpy":
            raise
        raise MissingOptionalDependencyError('fragment_arrays() needs numpy. Install it with: pip install "peptacular[numpy]"') from exc
    return numpy


@lru_cache(maxsize=4096)
def _extras(isotopes_key: Any, deltas_key: Any, monoisotopic: bool) -> tuple[int, str, str, float]:
    """``(isotope, isotope_label, delta_label, delta_mass)`` for a fragment's isotope and delta fields."""
    isotopes = dict(isotopes_key) if isinstance(isotopes_key, tuple) else isotopes_key
    deltas = dict(deltas_key) if isinstance(deltas_key, tuple) else deltas_key
    probe = Fragment(IonType.B, None, 0.0, monoisotopic, 0, isotopes=isotopes, deltas=deltas)
    raw = probe._isotopes
    if isinstance(raw, int):
        isotope = raw
        isotope_label = _format_counts([("13C", raw)]) if raw else ""
    else:
        isotope = int((raw or {}).get("13C", 0))
        isotope_label = _format_counts((raw or {}).items())
    delta_mass = 0.0
    fragment_deltas = probe.deltas
    for key, count in fragment_deltas.items():
        if isinstance(key, ChargedFormula):
            delta_mass += key.get_mass(monoisotopic=monoisotopic) * count
        else:
            delta_mass += key * count
    return isotope, isotope_label, _format_counts(fragment_deltas.items()), delta_mass


def _key(value: Mapping[Any, int] | int | None) -> Any:
    """Hashable, order-keeping form of an isotope or delta mapping."""
    return tuple(value.items()) if isinstance(value, Mapping) else value


def _fragment_extras(fragment: Fragment) -> tuple[int, str, str, float]:
    return _extras(_key(fragment._isotopes), _key(fragment._deltas), fragment.monoisotopic)


@lru_cache(maxsize=256)
def _charge_info(charge: int, monoisotopic: bool) -> tuple[float, int] | None:
    """``(carrier mass, external charge)`` that ``fragment()``'s fast series uses for an int charge.

    None when that charge takes ``fragment()``'s slow path (a carrier with a negative element
    count, e.g. a deprotonation).
    """
    from .annotation import ProFormaAnnotation

    carriers = ProFormaAnnotation.parse("G").set_charge(charge, inplace=False).charge_adducts
    if any(count < 0 for mod in carriers for count in mod.get_composition().values()):
        return None
    return carriers.get_mass(monoisotopic=monoisotopic) + proton_binding_offset(carriers, monoisotopic), carriers.get_charge()


def _is_fast_ion(ion_type: IonType) -> bool:
    info = FRAGMENT_ION_LOOKUP[ion_type]
    return (
        (info.is_forward or info.is_backward) and ion_type not in FRAGMENT_RULES and ion_type not in SATELLITE_TRIM_END and ion_type not in SATELLITE_TRIM_START
    )


def fragment_arrays(
    annotations: Sequence["ProFormaAnnotation"],
    ion_types: Any,
    charges: Any,
    *,
    monoisotopic: bool,
    isotopes: Any,
    deltas: Any,
    neutral_deltas: Any,
    calculate_with_composition: bool,
    max_ndeltas: int,
    min_length: int | None = None,
    max_length: int | None = None,
) -> dict[str, "np.ndarray"]:
    """Fragment every annotation and return the ions as columns (see :data:`FRAGMENT_ARRAY_KEYS`).

    Takes the arguments of :meth:`ProFormaAnnotation.fragment`; ``peptide_index`` is the
    annotation's index in ``annotations``.
    """
    np = _require_numpy()
    from .annotation import _as_options, get_loss_combinations

    kwargs: dict[str, Any] = {
        "monoisotopic": monoisotopic,
        "isotopes": isotopes,
        "deltas": deltas,
        "neutral_deltas": neutral_deltas,
        "calculate_with_composition": calculate_with_composition,
        "max_ndeltas": max_ndeltas,
        "min_length": min_length,
        "max_length": max_length,
    }

    # Settings shared by every peptide decide whether the numpy path can apply at all.
    fast_ions: list[IonType] = []
    products: list[tuple[float, float, tuple[int, str, str, float]]] = []
    fast = not calculate_with_composition and not any(nd is not None for nd in _as_options(neutral_deltas))
    if fast:
        fast_ions = [to_ion_type(ion) for ion in _as_options(ion_types)]
        fast = all(_is_fast_ion(ion) for ion in fast_ions)
    if fast:
        isotope_infos = [IsotopeInfo.from_input(isotope) for isotope in _as_options(isotopes)]
        delta_infos = [DeltaInfo.from_input(delta) for delta in deltas]
        no_loss = get_loss_combinations({}, max_ndeltas)
        for isotope in isotope_infos:
            for delta in delta_infos:
                for ndelta in no_loss:
                    combined = delta + ndelta
                    if isotope.data or any(isinstance(key, ChargedFormula) for key in combined.deltas):
                        fast = False
                    iso_mass = isotope.get_mass_delta(monoisotopic)
                    delta_mass = combined.get_mass_delta(monoisotopic)
                    extras = _extras(_key(isotope.to_fragment_mapping), _key(combined.to_fragment_mapping), monoisotopic)
                    products.append((iso_mass, delta_mass, extras))
    ion_masses = [
        (FRAGMENT_ION_LOOKUP[ion].ion_type, FRAGMENT_ION_LOOKUP[ion].is_forward, _ion_mass(FRAGMENT_ION_LOOKUP[ion].ion_type, monoisotopic))
        for ion in fast_ions
    ]

    counts = [0] * len(annotations)
    object_rows: dict[int, list[Fragment]] = {}
    buckets: dict[tuple[int, tuple[tuple[float, int], ...]], list[tuple[int, list[float]]]] = {}
    for index, annotation in enumerate(annotations):
        plan = _fast_plan(annotation, charges, monoisotopic, min_length, max_length) if fast else None
        if plan is None:
            fragments = annotation.fragment(ion_types, charges, **kwargs)
            object_rows[index] = fragments
            counts[index] = len(fragments)
            continue
        masses, charge_infos, n_positions = plan
        counts[index] = len(charge_infos) * len(ion_masses) * n_positions * len(products)
        buckets.setdefault((len(masses), charge_infos), []).append((index, masses))

    offsets = np.zeros(len(annotations) + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    total = int(offsets[-1])
    out: dict[str, Any] = {}
    for key in _INT_KEYS:
        out[key] = np.zeros(total, dtype=np.int64)
    for key in _FLOAT_KEYS:
        out[key] = np.zeros(total, dtype=np.float64)
    for key in _STR_KEYS:
        out[key] = np.empty(total, dtype=object)

    if object_rows:
        _fill_objects(np, out, object_rows, offsets)

    lo = 1 if min_length is None else max(1, min_length)
    for (n, charge_infos), members in buckets.items():
        hi = n if max_length is None else min(n, max_length)
        columns = np.arange(lo - 1, hi)
        matrix = np.array([masses for _, masses in members], dtype=np.float64)
        forward = np.cumsum(matrix, axis=1)[:, columns]
        backward = np.cumsum(matrix[:, ::-1], axis=1)[:, columns]
        shape = (len(members), len(charge_infos), len(ion_masses), len(columns), len(products))
        mass = np.empty(shape, dtype=np.float64)
        mz = np.empty(shape, dtype=np.float64)
        for ci, (charge_mass, external_charge) in enumerate(charge_infos):
            electrons = external_charge * ELECTRON_MASS
            for ti, (_, is_forward, ion_mass) in enumerate(ion_masses):
                cumulative = forward if is_forward else backward
                for pi, (iso_mass, delta_mass, _) in enumerate(products):
                    # Same arithmetic order as fragment()'s fast series.
                    value = cumulative + iso_mass
                    value += delta_mass
                    value += charge_mass
                    value += ion_mass
                    value -= electrons
                    mass[:, ci, ti, :, pi] = value
                    mz[:, ci, ti, :, pi] = value / abs(external_charge) if external_charge != 0 else value
        bad = ~np.isfinite(mass) | (mass < 0)
        if bad.any():
            first = members[int(np.argmax(bad.reshape(len(members), -1).any(axis=1)))][0]
            annotations[first].fragment(ion_types, charges, **kwargs)  # raises the real error
            validate_mass(float(mass[bad][0]))

        pattern_shape = shape[1:]
        patterns: dict[str, Any] = {
            "ion_type": np.broadcast_to(np.array([ion.value for ion, _, _ in ion_masses], dtype=object)[None, :, None, None], pattern_shape),
            "position": np.broadcast_to(np.arange(lo, hi + 1, dtype=np.int64)[None, None, :, None], pattern_shape),
            "charge_state": np.broadcast_to(np.array([c for _, c in charge_infos], dtype=np.int64)[:, None, None, None], pattern_shape),
        }
        for key, column in zip(
            ("isotope", "isotope_label", "delta_label", "delta_mass"), zip(*(extras for _, _, extras in products), strict=True), strict=True
        ):
            dtype = object if key.endswith("label") else (np.float64 if key == "delta_mass" else np.int64)
            patterns[key] = np.broadcast_to(np.array(column, dtype=dtype), pattern_shape)

        starts = offsets[[index for index, _ in members]]
        rows = starts[:, None] + np.arange(int(np.prod(pattern_shape)), dtype=np.int64)[None, :]
        flat = rows.ravel()
        out["peptide_index"][flat] = np.repeat(np.array([index for index, _ in members], dtype=np.int64), rows.shape[1])
        out["mass"][flat] = mass.ravel()
        out["mz"][flat] = mz.ravel()
        for key, pattern in patterns.items():
            out[key][flat] = np.broadcast_to(pattern.reshape(1, -1), rows.shape).ravel()
    return {key: out[key] for key in FRAGMENT_ARRAY_KEYS}


def _fast_plan(
    annotation: "ProFormaAnnotation", charges: Any, monoisotopic: bool, min_length: int | None, max_length: int | None
) -> tuple[list[float], tuple[tuple[float, int], ...], int] | None:
    """Residue masses, per-charge carrier info and ion count per series, or None for the object path."""
    from .annotation import _as_options

    n = len(annotation)
    if n == 0:
        return None
    options = annotation._default_fragment_charges(annotation.charge_state) if charges is None else _as_options(charges)
    charge_infos: list[tuple[float, int]] = []
    for charge in options:
        if isinstance(charge, bool) or not isinstance(charge, int):
            return None
        info = _charge_info(charge, monoisotopic)
        if info is None:
            return None
        charge_infos.append(info)
    masses = annotation._series_mass_vector(monoisotopic, False)
    if masses is None:
        return None
    lo = 1 if min_length is None else max(1, min_length)
    hi = n if max_length is None else min(n, max_length)
    return masses, tuple(charge_infos), max(0, hi - lo + 1)


def _fill_objects(np: Any, out: dict[str, Any], object_rows: Mapping[int, list[Fragment]], offsets: Any) -> None:
    index: list[int] = []
    columns: dict[str, list[Any]] = {key: [] for key in FRAGMENT_ARRAY_KEYS}
    for peptide, fragments in object_rows.items():
        start = int(offsets[peptide])
        index.extend(range(start, start + len(fragments)))
        for fragment in fragments:
            position, end_position = fragment.position if isinstance(fragment.position, tuple) else (fragment.position, None)
            isotope, isotope_label, delta_label, delta_mass = _fragment_extras(fragment)
            columns["peptide_index"].append(peptide)
            columns["ion_type"].append(fragment.ion_type.value)
            columns["position"].append(position or 0)
            columns["end_position"].append(end_position or 0)
            columns["charge_state"].append(fragment.charge_state)
            columns["mz"].append(fragment.mz)
            columns["mass"].append(fragment.mass)
            columns["isotope"].append(isotope)
            columns["isotope_label"].append(isotope_label)
            columns["delta_label"].append(delta_label)
            columns["delta_mass"].append(delta_mass)
    rows = np.array(index, dtype=np.int64)
    for key, values in columns.items():
        out[key][rows] = np.array(values, dtype=out[key].dtype)

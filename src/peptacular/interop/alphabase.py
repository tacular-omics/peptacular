"""Conversions for AlphaBase peptide-table columns."""

import warnings
from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Any

from peptacular.annotation import ProFormaAnnotation
from peptacular.proforma_components import ModificationTags, TagName

from ._errors import InteropConversionError, LossyConversionWarning
from ._optional import require_dependency
from ._policy import LossPolicy

if TYPE_CHECKING:
    import pandas as pd

AlphaBaseRow = dict[str, str | int | None]


def _handle_loss(messages: list[str], policy: LossPolicy) -> None:
    if not messages:
        return
    message = "AlphaBase cannot represent: " + "; ".join(messages)
    if policy is LossPolicy.ERROR:
        raise InteropConversionError(message)
    if policy is LossPolicy.WARN:
        warnings.warn(message, LossyConversionWarning, stacklevel=3)


def _alphabase_mod_name(modification: Any, site: str) -> str | None:
    if not isinstance(modification, ModificationTags):
        return None
    for tag in modification.tags:
        if isinstance(tag, TagName) and tag.position_id is None and tag.score is None:
            return f"{tag.name}@{site}"
    return None


def _append_mods(
    output_mods: list[str],
    output_sites: list[int],
    source: Any,
    site_name: str,
    site_index: int,
    valid_mods: Any,
    losses: list[str],
) -> None:
    for mod, count in source.parse_items():
        name = _alphabase_mod_name(mod, site_name)
        if name is None:
            losses.append(f"modification {mod!s} at site {site_index}")
            continue
        if name not in valid_mods:
            losses.append(f"modification {mod!s} at site {site_index} (unregistered AlphaBase key {name!r})")
            continue
        for _ in range(count):
            output_mods.append(name)
            output_sites.append(site_index)


def to_alphabase_row(
    annotation: ProFormaAnnotation,
    *,
    loss_policy: LossPolicy | str = LossPolicy.ERROR,
) -> AlphaBaseRow:
    """Convert an annotation to one AlphaBase precursor-table row.

    Fixed modifications are expanded onto individual residues. Advanced
    ProForma constructs are rejected by default and may only be discarded by
    explicitly selecting ``"warn"`` or ``"drop"``.

    :param annotation: Peptacular annotation to convert.
    :type annotation: ProFormaAnnotation
    :param loss_policy: Behavior for unsupported information.
    :type loss_policy: LossPolicy | str
    :return: Mapping containing ``sequence``, ``mods``, ``mod_sites``, and ``charge``.
    :rtype: AlphaBaseRow
    :raises InteropConversionError: If strict conversion would lose information.
    """
    policy = LossPolicy(loss_policy)
    alphabase_modification = require_dependency("alphabase.constants.modification", "alphabase")
    valid_mods = alphabase_modification.MOD_DF.index
    working = annotation.copy()
    losses: list[str] = []

    if working.compound_name or working.ion_name or working.peptide_name:
        losses.append("annotation names")
    if working.has_isotope_mods:
        losses.append("global isotope modifications")
    if working.has_labile_mods:
        losses.append("labile modifications")
    if working.has_unknown_mods:
        losses.append("unlocalized modifications")
    if working.has_intervals:
        losses.append("ambiguous intervals")
    if working.ambiguous_residues:
        losses.append("ambiguous residues")
    if working.charge_type.value == "adducts":
        losses.append("charge adduct identities")

    if working.has_static_mods:
        try:
            working.condense_static_mods(inplace=True)
        except Exception as exc:
            losses.append(f"fixed modifications ({exc})")

    mods: list[str] = []
    sites: list[int] = []
    _append_mods(mods, sites, working.nterm_mods, "Any_N-term", 0, valid_mods, losses)
    for index, residue_mods in sorted(working.internal_mods.items()):
        _append_mods(mods, sites, residue_mods, working.sequence[index], index + 1, valid_mods, losses)
    _append_mods(mods, sites, working.cterm_mods, "Any_C-term", -1, valid_mods, losses)

    _handle_loss(losses, policy)
    return {
        "sequence": working.sequence,
        "mods": ";".join(mods),
        "mod_sites": ";".join(str(site) for site in sites),
        "charge": working.charge_state or None,
    }


def _from_alphabase_fields(
    sequence: str,
    mods: str = "",
    mod_sites: str = "",
    charge: int | None = None,
) -> ProFormaAnnotation:
    """Build a Peptacular annotation from AlphaBase peptide columns.

    :param sequence: Unmodified peptide sequence.
    :type sequence: str
    :param mods: Semicolon-delimited AlphaBase modification names.
    :type mods: str
    :param mod_sites: Semicolon-delimited AlphaBase modification sites.
    :type mod_sites: str
    :param charge: Optional precursor charge.
    :type charge: int | None
    :return: Peptacular annotation.
    :rtype: ProFormaAnnotation
    :raises InteropConversionError: If columns are inconsistent or invalid.
    """
    mod_values = [] if not mods else mods.split(";")
    site_values = [] if not mod_sites else mod_sites.split(";")
    if len(mod_values) != len(site_values):
        raise InteropConversionError(f"AlphaBase mods and mod_sites must have equal lengths, got {len(mod_values)} and {len(site_values)}")

    alphabase_modification = require_dependency("alphabase.constants.modification", "alphabase")
    valid_mods = alphabase_modification.MOD_DF.index
    annotation = ProFormaAnnotation(sequence=sequence, charge=charge)
    for raw_mod, raw_site in zip(mod_values, site_values, strict=True):
        try:
            site = int(raw_site)
        except ValueError as exc:
            raise InteropConversionError(f"Invalid AlphaBase modification site {raw_site!r}") from exc

        name, separator, target = raw_mod.rpartition("@")
        if not separator or not name or not target:
            raise InteropConversionError(f"Invalid AlphaBase modification {raw_mod!r}; expected 'Name@Target'")
        if raw_mod not in valid_mods:
            raise InteropConversionError(f"AlphaBase modification {raw_mod!r} is not registered")

        if site == 0:
            if target not in {"Any_N-term", "Protein_N-term"}:
                raise InteropConversionError(f"Modification {raw_mod!r} does not match N-terminal site 0")
            annotation.append_nterm_mod(name, inplace=True)
        elif site == -1:
            if target not in {"Any_C-term", "Protein_C-term"}:
                raise InteropConversionError(f"Modification {raw_mod!r} does not match C-terminal site -1")
            annotation.append_cterm_mod(name, inplace=True)
        else:
            index = site - 1
            if index < 0 or index >= len(sequence):
                raise InteropConversionError(f"AlphaBase modification site {site} is outside sequence length {len(sequence)}")
            if target != sequence[index]:
                raise InteropConversionError(f"Modification {raw_mod!r} targets {target!r}, but residue {site} is {sequence[index]!r}")
            annotation.append_internal_mod_at_index(index, name, inplace=True)
    return annotation


def from_alphabase_row(row: Mapping[str, Any]) -> ProFormaAnnotation:
    """Build a Peptacular annotation from one AlphaBase table row.

    :param row: Mapping with AlphaBase peptide columns.
    :type row: Mapping[str, Any]
    :return: Peptacular annotation.
    :rtype: ProFormaAnnotation
    :raises InteropConversionError: If required columns are missing or invalid.
    """
    if "sequence" not in row:
        raise InteropConversionError("AlphaBase row is missing required 'sequence' column")
    sequence = row["sequence"]
    mods = row.get("mods", "")
    mod_sites = row.get("mod_sites", "")
    raw_charge = row.get("charge")
    if not isinstance(sequence, str) or not isinstance(mods, str) or not isinstance(mod_sites, str):
        raise InteropConversionError("AlphaBase sequence, mods, and mod_sites values must be strings")
    try:
        charge = None if raw_charge in (None, "") else int(raw_charge)
    except (TypeError, ValueError) as exc:
        raise InteropConversionError(f"Invalid AlphaBase charge {raw_charge!r}") from exc
    return _from_alphabase_fields(sequence, mods, mod_sites, charge)


def to_alphabase_dataframe(
    annotations: Iterable[ProFormaAnnotation],
    *,
    loss_policy: LossPolicy | str = LossPolicy.ERROR,
) -> "pd.DataFrame":
    """Convert annotations to AlphaBase's native precursor DataFrame.

    AlphaBase's public peptide and precursor representation is a pandas
    DataFrame. The result is passed through
    :func:`alphabase.peptide.precursor.refine_precursor_df` and can be assigned
    directly to ``SpecLibBase.precursor_df``.

    :param annotations: Peptacular annotations to convert.
    :type annotations: Iterable[ProFormaAnnotation]
    :param loss_policy: Behavior for unsupported information.
    :type loss_policy: LossPolicy | str
    :return: Refined pandas DataFrame with AlphaBase columns.
    :rtype: pandas.DataFrame
    :raises InteropConversionError: If charged and uncharged rows are mixed.
    """
    pandas = require_dependency("pandas", "alphabase")
    precursor = require_dependency("alphabase.peptide.precursor", "alphabase")
    rows = [to_alphabase_row(annotation, loss_policy=loss_policy) for annotation in annotations]
    columns = ["sequence", "mods", "mod_sites", "charge"]
    dataframe = pandas.DataFrame(rows, columns=columns)

    if rows:
        has_charge = [row["charge"] is not None for row in rows]
        if any(has_charge) and not all(has_charge):
            raise InteropConversionError("AlphaBase DataFrame conversion cannot mix charged and uncharged annotations")
        if not any(has_charge):
            dataframe.drop(columns="charge", inplace=True)
    else:
        dataframe.drop(columns="charge", inplace=True)

    return precursor.refine_precursor_df(dataframe, ensure_data_validity=True)


def from_alphabase_dataframe(dataframe: "pd.DataFrame") -> list[ProFormaAnnotation]:
    """Convert an AlphaBase precursor DataFrame to Peptacular annotations.

    :param dataframe: AlphaBase pandas DataFrame.
    :type dataframe: pandas.DataFrame
    :return: Annotations in DataFrame row order.
    :rtype: list[ProFormaAnnotation]
    :raises TypeError: If *dataframe* is not a pandas DataFrame.
    :raises InteropConversionError: If required columns are absent.
    """
    pandas = require_dependency("pandas", "alphabase")
    require_dependency("alphabase.peptide.precursor", "alphabase")
    if not isinstance(dataframe, pandas.DataFrame):
        raise TypeError(f"Expected pandas.DataFrame, got {type(dataframe).__name__}")
    missing = {"sequence", "mods", "mod_sites"} - set(dataframe.columns)
    if missing:
        raise InteropConversionError(f"AlphaBase DataFrame is missing required columns: {', '.join(sorted(missing))}")
    return [from_alphabase_row(row) for row in dataframe.to_dict(orient="records")]

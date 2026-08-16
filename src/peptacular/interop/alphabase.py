"""Conversions for AlphaBase peptide-table columns."""

import warnings
from dataclasses import asdict, dataclass
from typing import Any

from peptacular.annotation import ProFormaAnnotation
from peptacular.proforma_components import ModificationTags, TagName

from ._errors import InteropConversionError, LossyConversionWarning
from ._optional import require_dependency
from ._policy import LossPolicy


@dataclass(frozen=True, slots=True)
class AlphaBasePeptide:
    """The AlphaBase columns that describe one peptide precursor."""

    sequence: str
    mods: str = ""
    mod_sites: str = ""
    charge: int | None = None

    def as_dict(self) -> dict[str, str | int | None]:
        """Return values suitable for a precursor DataFrame row.

        :return: Mapping with AlphaBase column names.
        :rtype: dict[str, str | int | None]
        """
        return asdict(self)


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


def to_alphabase(
    annotation: ProFormaAnnotation,
    *,
    loss_policy: LossPolicy | str = LossPolicy.ERROR,
) -> AlphaBasePeptide:
    """Convert an annotation to AlphaBase peptide-table values.

    Fixed modifications are expanded onto individual residues. Advanced
    ProForma constructs are rejected by default and may only be discarded by
    explicitly selecting ``"warn"`` or ``"drop"``.

    :param annotation: Peptacular annotation to convert.
    :type annotation: ProFormaAnnotation
    :param loss_policy: Behavior for unsupported information.
    :type loss_policy: LossPolicy | str
    :return: AlphaBase peptide columns.
    :rtype: AlphaBasePeptide
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
    return AlphaBasePeptide(
        sequence=working.sequence,
        mods=";".join(mods),
        mod_sites=";".join(str(site) for site in sites),
        charge=working.charge_state or None,
    )


def from_alphabase(
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

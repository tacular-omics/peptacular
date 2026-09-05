"""Scientific output schemas, including explicit units and partial field failures."""

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, create_model

from .contracts import Diagnostic, Envelope


class Row(BaseModel):
    # Additional fields preserve upstream lineage through chained operations.
    model_config = ConfigDict(extra="allow", allow_inf_nan=False)
    source_id: str | None = None
    source_key: str
    source_index: int
    row_key: str
    status: Literal["complete", "partial", "error"]
    diagnostics: list[Diagnostic] = Field(default_factory=list)
    proforma: str | None = None
    sequence: str | None = None


class InspectionRow(Row):
    length: int | None = None
    names: dict[str, str | None] | None = None
    modifications: dict[str, Any] | None = None
    encoded_charge: int | list[str] | None = None
    mass_ambiguous_residues: list[str] | None = None


class AnalysisRow(Row):
    length: int | None = None
    neutral_mass_da: float | None = None
    ion_mass_da: float | None = None
    mz: float | None = None
    charge: int | None = None
    external_charge: int | None = None
    intrinsic_charge: int | None = None
    composition: dict[str, float] | None = None
    property: float | None = None
    residue_counts: dict[str, int] | None = None


class FragmentRow(AnalysisRow):
    ion_series: Literal["a", "b", "c", "x", "y", "z", "p"] | None = None
    ordinal: int | None = None
    start: int | None = None
    end: int | None = None
    label: str | None = None
    isotopes: dict[str, int] | None = None
    losses: list[dict[str, Any]] | None = None
    monoisotopic: bool | None = None


class Difference(BaseModel):
    input: float
    reference: float
    delta_input_minus_reference: float


class ComparisonRow(Row):
    reference_id: str | None = None
    reference_proforma: str | None = None
    same_sequence: bool | None = None
    changed_annotation_fields: list[str] | None = None
    neutral_mass_da: Difference | None = None
    mz: Difference | None = None
    composition_delta: dict[str, float] | None = None


class IsotopeRow(Row):
    axis: Literal["neutral_mass_da", "ion_mass_da", "mz", "neutron_offset"] | None = None
    peak_index: int | None = None
    position: float | None = None
    relative_abundance: float | None = None
    neutron_offset: int | None = None
    normalization: str | None = None
    distribution_is_approximate: bool | None = None
    retained_probability: float | None = None


class DigestionRow(Row):
    length: int | None = None
    start: int | None = None
    end: int | None = None
    missed_cleavages: int | None = None
    specificity: Literal["full", "semi", "nonspecific"] | None = None
    enzyme: str | None = None


class EditRow(Row):
    original_proforma: str | None = None
    annotation: dict[str, Any] | None = None


class CandidateRow(Row):
    candidate_index: int | None = None
    modifications: dict[str, Any] | None = None


class MappingRow(Row):
    protein_id: str | None = None
    protein_key: str | None = None
    start: int | None = None
    end: int | None = None
    mapped: bool | None = None
    ambiguous: bool | None = None
    match_count: int | None = None
    coverage_fraction_for_this_peptide: float | None = None


class ConversionRow(Row):
    target: Literal["proforma", "stable_json", "alphabase_row", "pyteomics", "psm_utils"] | None = None
    value: str | dict[str, Any] | None = None


ROWS = {
    "inspect_peptides": InspectionRow,
    "analyze_peptides": AnalysisRow,
    "fragment_peptides": FragmentRow,
    "compare_peptides": ComparisonRow,
    "isotope_envelopes": IsotopeRow,
    "digest_proteins": DigestionRow,
    "edit_peptides": EditRow,
    "enumerate_modifications": CandidateRow,
    "map_peptides": MappingRow,
    "convert_annotations": ConversionRow,
}


# Job submissions and preflights contain execution metadata instead of scientific rows.
class ExecutionRow(BaseModel):
    model_config = ConfigDict(extra="allow")
    mode: Literal["inline", "job"]
    input_records: int


class JobStateRow(BaseModel):
    model_config = ConfigDict(extra="allow")
    state: Literal["queued", "running", "succeeded", "partially_succeeded", "failed", "cancelled", "interrupted", "unavailable"]


OUTPUTS = {
    name: create_model(f"{row.__name__}Envelope", __base__=Envelope, records=(list[row | ExecutionRow | JobStateRow], Field(default_factory=list)))
    for name, row in ROWS.items()
}

"""Typed requests for stateless peptide calculations."""

from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

Count = Annotated[int, Field(strict=True, ge=1, le=5000)]
Index = Annotated[int, Field(strict=True, ge=0, le=1000000)]
Charge = Annotated[int, Field(strict=True, ge=-20, le=20)]
Number = Annotated[float, Field(strict=True, allow_inf_nan=False)]
Text = Annotated[str, Field(min_length=1, max_length=10000)]


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False, strict=True)


class Record(Contract):
    id: Annotated[str, Field(max_length=500)] | None = None
    annotation: Text | dict[str, Any]


Inputs = Annotated[list[Record], Field(min_length=1, max_length=100)]


class Scientific(Contract):
    inputs: Inputs
    max_rows: Count = 1000


class Inspect(Scientific):
    detail: Literal["summary", "full"] = "summary"


class ChargeSelection(Contract):
    charges: (
        Annotated[list[Charge], Field(min_length=1, max_length=10, description="External carrier charges. Total ion charge also includes intrinsic charge.")]
        | None
    ) = None
    charge_policy: Literal["require_agreement", "override"] = "require_agreement"


class ChargeSettings(ChargeSelection):
    monoisotopic: bool = True


class PropertySettings(Contract):
    scale: Annotated[str, Field(min_length=1, max_length=100)] = "hphob_kyte_doolittle"
    aggregation: Literal["avg", "sum"] = "avg"
    missing_residues: Literal["error", "ignore", "zero"] = "error"
    modification_treatment: Literal["ignore"] = "ignore"


class Analyze(Scientific, ChargeSettings):
    measurements: Annotated[
        list[Literal["length", "neutral_mass_da", "ion_mass_da", "mz", "composition", "property", "residue_counts"]],
        Field(min_length=1, max_length=7),
    ] = ["length", "neutral_mass_da"]
    property_settings: PropertySettings = Field(default_factory=PropertySettings)


class Delta(Contract):
    kind: Literal["formula", "mass"]
    value: Text | Number

    @model_validator(mode="after")
    def check_kind(self):
        if (self.kind == "formula") != isinstance(self.value, str):
            raise ValueError("Formula deltas require text, mass deltas require a finite number")
        return self


class Fragment(Scientific, ChargeSettings):
    ion_series: Annotated[list[Literal["a", "b", "c", "x", "y", "z", "p"]], Field(min_length=1, max_length=7)] = ["b", "y"]
    isotope_offsets: Annotated[list[Annotated[int, Field(strict=True, ge=0, le=10)]], Field(min_length=1, max_length=10)] = [0]
    deltas: Annotated[list[Delta], Field(max_length=8)] = []
    include: list[Literal["composition", "sequence", "label"]] = ["label"]
    min_mz: Number | None = None
    max_mz: Number | None = None

    @model_validator(mode="after")
    def check_range(self):
        if self.min_mz is not None and self.max_mz is not None and self.min_mz > self.max_mz:
            raise ValueError("min_mz exceeds max_mz")
        return self


class Compare(Scientific, ChargeSettings):
    reference: Record
    measurements: list[Literal["annotation", "composition", "neutral_mass_da", "mz"]] = ["annotation", "neutral_mass_da"]


class Isotopes(Scientific, ChargeSelection):
    axis: Literal["neutral_mass_da", "ion_mass_da", "mz", "neutron_offset"] = "neutral_mass_da"
    max_peaks: Annotated[int, Field(strict=True, ge=1, le=100)] = 10
    min_relative_abundance: Annotated[float, Field(strict=True, gt=0, le=1)] = 0.001
    resolution: Annotated[int, Field(strict=True, ge=0, le=6)] = 5


class Digest(Scientific):
    enzyme: Annotated[str, Field(min_length=1, max_length=100)] = "trypsin"
    specificity: Literal["full", "semi", "nonspecific"] = "full"
    missed_cleavages: Annotated[int, Field(strict=True, ge=0, le=10)] = 0
    min_length: Count = 1
    max_length: Count = 100

    @model_validator(mode="after")
    def check_lengths(self):
        if self.min_length > self.max_length:
            raise ValueError("min_length exceeds max_length")
        return self


class ModificationEdit(Contract):
    action: Literal["add", "remove", "clear"]
    location: Literal["internal", "nterm", "cterm", "labile", "unknown", "static", "isotope"]
    index: Index | None = None
    modification: Text | None = None

    @model_validator(mode="after")
    def check_location(self):
        if (self.location == "internal") != (self.index is not None):
            raise ValueError("Only internal edits require a zero-based index")
        if (self.action != "clear") != (self.modification is not None):
            raise ValueError("Add/remove require a modification, clear does not")
        return self


class ChargeEdit(Contract):
    action: Literal["charge"]
    charge: Charge | None


class SliceEdit(Contract):
    action: Literal["slice"]
    start: Index
    end: Index


class ExpandEdit(Contract):
    action: Literal["expand_static"]


class Edit(Scientific):
    edits: Annotated[list[ModificationEdit | ChargeEdit | SliceEdit | ExpandEdit], Field(min_length=1, max_length=30)]


class Rule(Contract):
    location: Literal["internal", "nterm", "cterm", "labile"] = "internal"
    residues: Annotated[str, Field(min_length=1, max_length=26, pattern="^[ACDEFGHIKLMNPQRSTVWY]+$")] | None = None
    modification: Annotated[str, Field(min_length=1, max_length=200)]
    variable: bool = True

    @model_validator(mode="after")
    def check_residues(self):
        if (self.location == "internal") != (self.residues is not None):
            raise ValueError("Internal rules require residues, terminal and labile rules do not")
        return self


class Enumerate(Scientific):
    rules: Annotated[list[Rule], Field(min_length=1, max_length=10)]
    max_variable_modifications: Annotated[int, Field(strict=True, ge=0, le=5)] = 2
    max_candidates: Annotated[int, Field(strict=True, ge=1, le=1000)] = 1000


class Map(Scientific):
    proteins: Inputs
    modification_policy: Literal["reject", "ignore"] = "reject"


class Convert(Scientific):
    target: Literal["proforma", "stable_json", "alphabase_row", "pyteomics", "psm_utils"]
    loss_policy: Literal["raise", "warn", "drop"] = "raise"


class FindModifications(Contract):
    query_type: Literal["accession", "name", "mass"]
    query: Text | Number
    vocabularies: Annotated[list[Literal["unimod", "psimod", "xlmod"]], Field(min_length=1, max_length=3)] = ["unimod"]
    name_mode: Literal["exact", "prefix", "contains"] = "contains"
    tolerance: Annotated[float, Field(strict=True, gt=0, le=100)] | None = None
    tolerance_unit: Literal["da", "ppm"] = "da"
    monoisotopic: bool = True
    limit: Annotated[int, Field(strict=True, ge=1, le=500)] = 50
    offset: Index = 0

    @model_validator(mode="after")
    def check_query(self):
        if self.query_type == "mass":
            if isinstance(self.query, str) or self.tolerance is None:
                raise ValueError("Mass lookup requires a number and explicit tolerance")
            if self.query == 0 and self.tolerance_unit == "ppm":
                raise ValueError("Zero mass requires a tolerance in Da")
        elif not isinstance(self.query, str) or self.tolerance is not None:
            raise ValueError("Name/accession lookup requires text and no mass tolerance")
        return self


class GetReference(Contract):
    topic: Literal["capabilities", "enzymes", "ions", "scales", "notation", "conventions", "schemas"] = "capabilities"
    search: Annotated[str, Field(max_length=100)] = ""
    offset: Index = 0
    limit: Annotated[int, Field(strict=True, ge=1, le=100)] = 25


class Diagnostic(Contract):
    code: str
    message: str
    stage: str = "calculate"
    field: str | None = None
    source_key: str | None = None
    recovery: str = "Review the input and applied settings."


class Computation(Contract):
    complete: bool = True
    stop_reason: str | None = None


class Envelope(Contract):
    contract_version: Literal["1.0"] = "1.0"
    request_id: str
    status: Literal["complete", "partial", "error"] = "complete"
    applied_settings: dict[str, Any] = {}
    records: list[dict[str, Any]] = []
    diagnostics: list[Diagnostic] = []
    returned_rows: int = 0
    total_rows: int | None = None
    next_offset: int | None = None
    computation: Computation = Field(default_factory=Computation)


SCIENTIFIC = {
    "inspect_peptides": Inspect,
    "analyze_peptides": Analyze,
    "fragment_peptides": Fragment,
    "compare_peptides": Compare,
    "isotope_envelopes": Isotopes,
    "digest_proteins": Digest,
    "edit_peptides": Edit,
    "enumerate_modifications": Enumerate,
    "map_peptides": Map,
    "convert_annotations": Convert,
}
REQUESTS = {
    **SCIENTIFIC,
    "get_reference": GetReference,
    "find_modifications": FindModifications,
}

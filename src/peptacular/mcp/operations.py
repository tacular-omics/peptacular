"""Scientific adapters for small, bounded in-memory requests."""

import importlib.metadata
import itertools
import json
import warnings
from collections import Counter
from dataclasses import asdict
from typing import Any, Literal

import peptacular as pt
from peptacular.diagnostics import diagnostic_from_exception

from . import contracts as c

RESULT_BYTES = 240000
LIMITS = {
    "records_per_input": 100,
    "request_bytes": 1048576,
    "annotation_characters": 100000,
    "rows": 5000,
    "result_row_bytes": RESULT_BYTES,
    "fragment_combinations_per_charge": 50000,
    "enumeration_residues_per_peptide": 200,
    "modification_candidates_per_peptide": 1000,
    "isotope_residues_per_peptide": 1000,
}


class ServiceError(ValueError):
    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def diagnostic(exc: Exception, field: str | None = None, stage: Literal["parse", "validate", "calculate"] = "calculate") -> dict:
    if isinstance(exc, ServiceError):
        code = exc.code
    elif isinstance(exc, ImportError):
        code = "missing_optional_dependency"
    else:
        code = diagnostic_from_exception(exc, stage).code
    return c.Diagnostic(code=code, message=str(exc)[:1000], stage=stage, field=field).model_dump()


def parse(value):
    a = pt.ProFormaAnnotation.from_dict(value) if isinstance(value, dict) else pt.parse(value)
    a.validate_sequence()
    a.validate_ambiguous_labels()
    return a


def composition(a):
    return {str(k): v for k, v in a.comp().items()}


def charged_annotations(a, settings):
    if settings.charges is None:
        return [a]
    if a.has_charge and settings.charge_policy == "require_agreement" and settings.charges != [a.charge_state]:
        raise ServiceError("charge_conflict", "Requested charges conflict with encoded charge. Select charge_policy='override' explicitly.")
    return [a.set_charge(z, inplace=False) for z in settings.charges]


def require_charge(a):
    if a.frag().charge_state == 0:
        raise ServiceError("missing_charge", "m/z requires a nonzero encoded or requested charge.")


def charge_fields(a):
    # Intrinsic charge may require resolving modification chemistry.
    f = a.frag()
    return {"charge": f.charge_state, "external_charge": f.external_charge, "intrinsic_charge": f.charge_state - f.external_charge}


def inspect_one(a, request):
    data = a.to_dict()
    row = {
        "proforma": a.serialize(),
        "sequence": a.sequence,
        "length": len(a),
        "names": data["names"],
        "modifications": data["modifications"],
        "intervals": data["intervals"],
        "encoded_charge": data["charge"],
        "mass_ambiguous_residues": list(a.mass_ambiguous_residues),
    }
    if request.detail == "full":
        row["annotation"] = data
        row["residue_counts"] = dict(Counter(a.sequence))
    yield row


def analyze_one(a, request):
    for ion in charged_annotations(a, request):
        row = {"proforma": ion.serialize(), "sequence": ion.sequence, "diagnostics": [], "external_charge": ion.charge_state}
        for field in request.measurements:
            try:
                if field == "length":
                    value = len(ion)
                elif field == "residue_counts":
                    value = dict(Counter(ion.sequence))
                elif field == "composition":
                    value = composition(ion)
                elif field == "property":
                    settings = request.property_settings
                    value = pt.calc_property(
                        ion.sequence,
                        scale=settings.scale,
                        aggregation_method=settings.aggregation,
                        missing_aa_handling="skip" if settings.missing_residues == "ignore" else settings.missing_residues,
                        method="sequential",
                    )
                    row["property_settings"] = settings.model_dump()
                elif field == "neutral_mass_da":
                    value = ion.neutral_mass(monoisotopic=request.monoisotopic)
                    row.update(charge_fields(ion))
                else:
                    if field == "mz":
                        require_charge(ion)
                    row.update(charge_fields(ion))
                    value = getattr(ion, "mass" if field == "ion_mass_da" else "mz")(monoisotopic=request.monoisotopic)
                row[field] = value
            except (ValueError, KeyError) as exc:
                row[field] = None
                row["diagnostics"].append(diagnostic(exc, field))
        yield row


def fragment_one(a, request):
    for ion in charged_annotations(a, request):
        require_charge(ion)
        combinations = len(ion) * len(request.ion_series) * len(request.isotope_offsets) * (len(request.deltas) + 1)
        if combinations > 50000:
            raise ServiceError("resource_limit", "Fragment expansion exceeds 50,000 fragments per charge. Reduce sequence length or settings.")
        fragments = ion.fragment(
            ion_types=request.ion_series,
            charges=[ion.charge_state],
            monoisotopic=request.monoisotopic,
            isotopes=request.isotope_offsets,
            deltas=[None, *[delta.value for delta in request.deltas]],
        )
        for f in fragments:
            if f.charge_state == 0:
                yield {
                    "proforma": ion.serialize(),
                    "ion_series": str(f.ion_type),
                    "ordinal": f.position,
                    "charge": 0,
                    "ion_mass_da": f.mass,
                    "mz": None,
                    "diagnostics": [diagnostic(ServiceError("missing_charge", "This fragment has zero total charge and no defined m/z."), "mz")],
                }
                continue
            if request.min_mz is not None and f.mz < request.min_mz:
                continue
            if request.max_mz is not None and f.mz > request.max_mz:
                continue
            ordinal = None if f.ion_type == "p" else f.position
            start = 0 if f.ion_type in ("a", "b", "c", "p") else len(ion) - int(ordinal or 0)
            end = len(ion) if f.ion_type in ("x", "y", "z", "p") else ordinal
            row: dict[str, Any] = {
                "proforma": ion.serialize(),
                "ion_series": str(f.ion_type),
                "ordinal": ordinal,
                "start": start,
                "end": end,
                "charge": f.charge_state,
                "external_charge": f.external_charge,
                "intrinsic_charge": f.charge_state - f.external_charge,
                "mz": f.mz,
                "ion_mass_da": f.mass,
                "neutral_mass_da": f.neutral_mass,
                "monoisotopic": request.monoisotopic,
                "losses": [{"value": str(k), "count": v} for k, v in f.losses.items()],
                "isotopes": {str(k): v for k, v in f.isotopes.items()},
                "diagnostics": [],
            }
            for field in request.include:
                try:
                    if field == "label":
                        value = f.serialize(format="mzpaf", include_sequence=False)
                    elif field == "composition":
                        value = {str(k): v for k, v in f.composition.items()}
                    else:
                        value = str(f.sequence)
                    row[field] = value
                except (ValueError, KeyError) as exc:
                    row[field] = None
                    row["diagnostics"].append(diagnostic(exc, field))
            yield row


def compare_one(a, request):
    reference = parse(request.reference.annotation)
    left = charged_annotations(a, request)
    right = charged_annotations(reference, request)
    for ion, ref in zip(left, right, strict=True):
        row = {"proforma": ion.serialize(), "reference_proforma": ref.serialize(), "reference_id": request.reference.id, "diagnostics": []}
        for field in request.measurements:
            try:
                if field == "annotation":
                    lhs, rhs = ion.to_dict(), ref.to_dict()
                    row["same_sequence"] = ion.sequence == ref.sequence
                    row["changed_annotation_fields"] = [key for key in lhs if lhs[key] != rhs[key]]
                    row["annotation"] = {"input": lhs, "reference": rhs}
                elif field == "composition":
                    lhs, rhs = composition(ion), composition(ref)
                    row["composition_delta"] = {key: lhs.get(key, 0) - rhs.get(key, 0) for key in sorted(lhs.keys() | rhs.keys())}
                else:
                    if field == "mz":
                        require_charge(ion)
                        require_charge(ref)
                    method = "mz" if field == "mz" else "neutral_mass"
                    lhs = getattr(ion, method)(monoisotopic=request.monoisotopic)
                    rhs = getattr(ref, method)(monoisotopic=request.monoisotopic)
                    row[field] = {"input": lhs, "reference": rhs, "delta_input_minus_reference": lhs - rhs}
            except (ValueError, KeyError) as exc:
                row["diagnostics"].append(diagnostic(exc, field))
        yield row


def isotopes_one(a, request):
    if len(a) > 1000:
        raise ServiceError("resource_limit", "Isotope calculations are limited to 1,000 residues per input.")
    for ion in charged_annotations(a, request):
        if request.axis == "mz":
            require_charge(ion)
        effective = ion.set_charge(0, inplace=False) if request.axis == "neutral_mass_da" else ion
        charges = charge_fields(effective)
        peaks = effective.isotopic_distribution(
            max_isotopes=request.max_peaks,
            min_abundance_threshold=request.min_relative_abundance,
            distribution_resolution=request.resolution,
            use_neutron_count=request.axis == "neutron_offset",
        )
        for index, peak in enumerate(peaks):
            yield {
                "proforma": ion.serialize(),
                "peak_index": index,
                "axis": request.axis,
                "position": peak.mass / abs(charges["charge"]) if request.axis == "mz" else peak.mass,
                "relative_abundance": peak.abundance,
                "neutron_offset": peak.neutron_count,
                **charges,
                "normalization": "maximum_retained_peak_equals_one",
                "distribution_is_approximate": True,
                "retained_probability": None,
                "truncation": {"max_peaks_during_convolution": request.max_peaks, "minimum_relative_abundance": request.min_relative_abundance},
            }


def resolve_enzyme(name):
    from tacular import PROTEASE_LOOKUP

    enzyme = PROTEASE_LOOKUP.query_id(name) or PROTEASE_LOOKUP.query_name(name)
    if enzyme is None:
        raise ServiceError("unknown_enzyme", "Unknown enzyme. Use get_reference(topic='enzymes') for supported identifiers.")
    return enzyme


def digest_one(a, request):
    enzyme = resolve_enzyme(request.enzyme)
    if request.specificity == "nonspecific":
        spans = pt.build_non_enzymatic_spans((0, len(a), 0), min_len=request.min_length, max_len=request.max_length)
    else:
        spans = a.digest(
            enzyme=enzyme.regex,
            missed_cleavages=request.missed_cleavages,
            semi=request.specificity == "semi",
            min_len=request.min_length,
            max_len=request.max_length,
        )
    for span in spans:
        peptide = a.slice(span.start, span.end)
        yield {
            "proforma": peptide.serialize(),
            "sequence": peptide.sequence,
            "length": len(peptide),
            "start": span.start,
            "end": span.end,
            "missed_cleavages": None if request.specificity == "nonspecific" else span.missed_cleavages,
            "enzyme": str(enzyme.id),
            "specificity": request.specificity,
        }


def edit_one(a, request):
    changed = a.copy()
    for edit in request.edits:
        if edit.action == "charge":
            changed.set_charge(edit.charge)
        elif edit.action == "slice":
            if not 0 <= edit.start < edit.end <= len(changed):
                raise ServiceError("invalid_coordinates", "Slices require 0 <= start < end <= sequence length.")
            changed = changed.slice(edit.start, edit.end)
        elif edit.action == "expand_static":
            changed.condense_static_mods()
        else:
            if edit.index is not None and edit.index >= len(changed):
                raise ServiceError("invalid_coordinates", "Modification index is outside the sequence.")
            prefix = {"add": "append", "remove": "remove", "clear": "clear"}[edit.action]
            if edit.location == "internal":
                suffix = "mod"
                args = [edit.index]
                method = f"{prefix}_internal_{suffix}_at_index"
            else:
                suffix = "mods" if edit.action == "clear" else "mod"
                args = []
                method = f"{prefix}_{edit.location}_{suffix}"
            if edit.action != "clear":
                args.append(edit.modification)
            getattr(changed, method)(*args)
    changed.validate_annotation()
    yield {"proforma": changed.serialize(), "sequence": changed.sequence, "original_proforma": a.serialize(), "annotation": changed.to_dict()}


def enumerate_one(a, request):
    if len(a) > 200:
        raise ServiceError("resource_limit", "Modification enumeration is limited to 200 residues per peptide.")
    settings = {}
    for rule in request.rules:
        name = f"{rule.location}_{'variable' if rule.variable else 'static'}"
        mapping = settings.setdefault(name, {})
        for residue in rule.residues or [None]:
            mapping.setdefault(residue, []).append(rule.modification)
    candidates = a.modify(**settings, max_variable_mods=request.max_variable_modifications, use_regex=False, inplace=False)
    for index, candidate in enumerate(candidates):
        if index == request.max_candidates:
            raise ServiceError("candidate_limit", "Candidate limit reached. Enumeration is incomplete.")
        yield {
            "proforma": candidate.serialize(),
            "sequence": candidate.sequence,
            "candidate_index": index,
            "modifications": candidate.to_dict()["modifications"],
        }


def convert_one(a, request):
    from peptacular import interop

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if request.target == "proforma":
            value = a.serialize()
        elif request.target == "stable_json":
            value = a.to_dict()
        elif request.target == "alphabase_row":
            value = interop.to_alphabase_row(a, loss_policy="error" if request.loss_policy == "raise" else "warn")
        elif request.target == "pyteomics":
            value = str(interop.to_pyteomics(a))
        else:
            value = str(interop.to_psm_utils(a).proforma)
    yield {
        "proforma": a.serialize(),
        "target": request.target,
        "value": value,
        "diagnostics": [
            c.Diagnostic(code="conversion_warning", message=str(w.message)[:1000]).model_dump()
            for w in caught
            if issubclass(w.category, interop.LossyConversionWarning)
        ],
    }


ADAPTERS = {
    "inspect_peptides": inspect_one,
    "analyze_peptides": analyze_one,
    "fragment_peptides": fragment_one,
    "compare_peptides": compare_one,
    "isotope_envelopes": isotopes_one,
    "digest_proteins": digest_one,
    "edit_peptides": edit_one,
    "enumerate_modifications": enumerate_one,
    "convert_annotations": convert_one,
}


def map_one(a, request, proteins):
    if request.modification_policy == "reject" and a.has_mods():
        raise ServiceError("modification_policy", "Sequence mapping requires modification_policy='ignore' for modified peptides.")
    found = []
    for protein in proteins:
        try:
            target = parse(protein["annotation"])
            if request.modification_policy == "reject" and target.has_mods():
                raise ServiceError("modification_policy", "Modified proteins require modification_policy='ignore'.")
        except (ValueError, KeyError) as exc:
            yield {"protein_key": protein["source_key"], "diagnostics": [diagnostic(exc)]}
            continue
        start = target.sequence.find(a.sequence)
        spans = []
        while start >= 0:
            spans.append((start, start + len(a)))
            if len(found) + len(spans) > request.max_rows:
                raise ServiceError("row_limit", "Mapping expansion exceeds the requested row budget.")
            start = target.sequence.find(a.sequence, start + 1)
        covered = set(itertools.chain.from_iterable(range(start, end) for start, end in spans))
        for start, end in spans:
            found.append(
                {
                    "proforma": a.serialize(),
                    "sequence": a.sequence,
                    "protein_id": protein.get("id"),
                    "protein_key": protein["source_key"],
                    "start": start,
                    "end": end,
                    "protein_length": len(target),
                    "coverage_fraction_for_this_peptide": len(covered) / len(target),
                }
            )
    if not found:
        yield {"proforma": a.serialize(), "sequence": a.sequence, "mapped": False, "match_count": 0}
    for row in found:
        yield {**row, "mapped": True, "match_count": len(found), "ambiguous": len(found) > 1}


def run_operation(name: str, payload: dict, records: list[dict], proteins: list[dict] | None = None) -> dict:
    request = c.SCIENTIFIC[name].model_validate(payload)
    rows, diagnostics = [], []
    complete = True
    stop_reason = None
    byte_count = 0
    for source in records:
        if len(rows) >= request.max_rows:
            complete, stop_reason = False, "row_limit"
            break
        stage: Literal["parse", "validate", "calculate"] = "parse"
        try:
            a = parse(source["annotation"])
            stage = "calculate"
            generated = map_one(a, request, proteins or []) if name == "map_peptides" else ADAPTERS[name](a, request)
            for row in generated:
                if len(rows) >= request.max_rows:
                    raise ServiceError("row_limit", "Requested result row budget reached.")
                row = {
                    **row,
                    "source_id": source.get("id"),
                    "source_key": source["source_key"],
                    "source_index": source["source_index"],
                    "original_annotation": source["annotation"],
                    "row_key": f"row_{len(rows)}",
                }
                row["status"] = "partial" if row.get("diagnostics") else "complete"
                if name == "digest_proteins":
                    row.update(
                        protein_key=source["source_key"],
                        protein_id=source.get("id"),
                        protein_record_index=source["source_index"],
                        protein_start=row["start"],
                        protein_end=row["end"],
                    )
                size = len(json.dumps(row, allow_nan=False).encode())
                if byte_count + size > RESULT_BYTES:
                    raise ServiceError("byte_limit", "Result byte budget reached. Request fewer or smaller fields.")
                byte_count += size
                rows.append(row)
        except (ValueError, KeyError, ImportError) as exc:
            if isinstance(exc, ServiceError) and exc.code in ("row_limit", "candidate_limit", "byte_limit", "resource_limit"):
                complete = False
                stop_reason = exc.code
                diagnostics.append({**diagnostic(exc), "source_key": source["source_key"]})
                if exc.code in ("row_limit", "byte_limit"):
                    break
            else:
                error_row = {
                    "source_id": source.get("id"),
                    "source_key": source["source_key"],
                    "source_index": source["source_index"],
                    "row_key": f"row_{len(rows)}",
                    "status": "error",
                    "diagnostics": [diagnostic(exc, stage=stage)],
                }
                byte_count += len(json.dumps(error_row, allow_nan=False).encode())
                if byte_count > RESULT_BYTES:
                    complete, stop_reason = False, "byte_limit"
                    break
                rows.append(error_row)
    return {
        "records": rows,
        "diagnostics": diagnostics,
        "computation": {"complete": complete, "stop_reason": stop_reason},
    }


def find_modifications(request):
    import tacular

    rows = []
    for vocabulary in request.vocabularies:
        lookup = getattr(tacular, f"{vocabulary.upper()}_LOOKUP")
        for entry in lookup:
            mass = entry.monoisotopic_mass if request.monoisotopic else entry.average_mass
            if request.query_type == "mass":
                tolerance = request.tolerance if request.tolerance_unit == "da" else abs(request.query) * request.tolerance / 1e6
                match = mass is not None and abs(mass - request.query) <= tolerance
            elif request.query_type == "accession":
                match = str(entry.id).casefold() == request.query.split(":")[-1].casefold()
            else:
                query, name = request.query.casefold(), entry.name.casefold()
                match = name == query if request.name_mode == "exact" else name.startswith(query) if request.name_mode == "prefix" else query in name
            if match:
                row = {
                    "vocabulary": vocabulary,
                    "accession": str(entry.id),
                    "name": entry.name,
                    "mass_da": mass,
                    "formula": getattr(entry, "formula", None),
                    "composition": getattr(entry, "dict_composition", None),
                }
                if request.query_type == "mass":
                    row["mass_error_da"] = mass - request.query
                    row["mass_error_ppm"] = (mass - request.query) / abs(request.query) * 1e6 if request.query else None
                rows.append(row)
    rows.sort(key=lambda r: (abs(r.get("mass_error_da", 0)), r["vocabulary"], r["accession"]))
    return rows


def versions():
    result = {"peptacular": pt.__version__}
    for package in ("tacular", "mcp", "pydantic", "pyteomics", "psm-utils", "alphabase"):
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = None
    return result


CONVENTIONS = {
    "coordinates": "start is zero-based, end is exclusive. Fragment ordinal counts residues. Precursor ordinal is null.",
    "charge": "Signed total, external carrier and intrinsic charges are separate. m/z requires nonzero charge. Conflicts require override.",
    "isotopes": "Theoretical approximate distributions. Relative maximum abundance is one. Retained probability is unavailable.",
    "ownership": "Peptacular calculates theoretical sequence properties. Spectacular owns observed spectra and spectrum matching.",
    "inputs": "Pass annotation records directly. IDs and source indexes identify records within this call. No server-side data is retained.",
    "limits": "Inspect computation.complete and stop_reason. For a truncated calculation, narrow the request or split the batch.",
}


def reference_rows(request, limits):
    if request.topic == "capabilities":
        return [{"tools": list(c.REQUESTS), "versions": versions(), "limits": limits, "transport": "local_stdio", "contract_version": "1.0"}]
    if request.topic == "enzymes":
        from tacular import PROTEASE_LOOKUP

        rows = [asdict(entry) for entry in PROTEASE_LOOKUP]
    elif request.topic == "scales":
        rows = [{"id": str(key), "aggregation": ["avg", "sum"], "modification_treatment": "ignore"} for key in pt.PROPERTY_SCALES]
    elif request.topic == "ions":
        rows = [{"ion_series": ion, "kind": "precursor" if ion == "p" else "backbone"} for ion in ("a", "b", "c", "x", "y", "z", "p")]
    elif request.topic == "schemas":
        rows = [{"tool": name, "schema": model.model_json_schema()} for name, model in c.REQUESTS.items()]
    elif request.topic == "notation":
        rows = [
            {
                "description": "Use ProForma 2.1 or versioned Peptacular JSON. Parsing does not guarantee calculation support.",
                "examples": ["PEPTIDE/2", "M[UNIMOD:35]PEPTIDE", "PEPT[+79.9663]IDE", "<[Carbamidomethyl]@C>ACDC"],
            }
        ]
    else:
        rows = [{"topic": key, "description": value} for key, value in CONVENTIONS.items()]
    return [row for row in rows if request.search.casefold() in json.dumps(row).casefold()]

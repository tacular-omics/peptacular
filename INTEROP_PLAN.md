# Optional package interoperability plan

## Goal

Make `ProFormaAnnotation` easy to exchange with commonly used Python
proteomics packages without adding any of those packages to Peptacular's core
runtime dependencies.

The first release should provide typed, tested adapters for Pyteomics,
`psm_utils`, and AlphaBase. Later releases can add pyOpenMS, rustyms, and
`spectrum_utils` after their supported ProForma subsets and stable APIs have
been verified.

## Design principles

1. Keep third-party imports out of `peptacular.__init__` and the core annotation
   implementation. Import an integration dependency only when its adapter is
   called.
2. Use ProForma text as the interchange representation when the other package
   has a public ProForma parser and serializer. Do not manually map its private
   object graph.
3. Never silently discard information. An adapter must either preserve the
   annotation, raise a descriptive `InteropConversionError`, or require the
   caller to explicitly select a lossy conversion policy.
4. Keep adapters usable without installing every integration. Type-only imports
   should be guarded by `TYPE_CHECKING`.
5. Test semantic round trips, not merely object construction. Canonical text is
   preferred where both libraries provide it; otherwise compare parsed
   Peptacular annotations.

## Proposed public API

Create an integration namespace rather than adding third-party methods to the
already-large `ProFormaAnnotation` class:

```text
peptacular.interop
├── _errors.py
├── _policy.py
├── pyteomics.py
├── psm_utils.py
├── alphabase.py
├── pyopenms.py          # phase 2
├── rustyms.py           # phase 2
└── spectrum_utils.py    # phase 2; outbound-first
```

Example usage:

```python
from peptacular.interop.pyteomics import from_pyteomics, to_pyteomics

annotation = peptacular.parse("[Acetyl]-PEM[Oxidation]TIDE/2")
pyteomics_value = to_pyteomics(annotation)
round_tripped = from_pyteomics(pyteomics_value)
```

Use package-local function names because the module already identifies the
target. This keeps imports readable and avoids a growing collection of methods
on the core model.

### Shared policy and errors

```python
class LossPolicy(StrEnum):
    ERROR = "error"
    WARN = "warn"
    DROP = "drop"


class InteropError(Exception): ...
class MissingOptionalDependencyError(InteropError, ImportError): ...
class InteropConversionError(InteropError, ValueError): ...
class LossyConversionWarning(UserWarning): ...
```

`ERROR` must be the default. `WARN` and `DROP` should only exist on adapters
whose target representation is known to be less expressive. The conversion
error should name the unsupported feature and target package.

## Phase 0: complete existing format symmetry

Before adding package adapters, finish or explicitly de-scope the existing
reverse conversions:

- `ProFormaAnnotation.to_ip2()`
- `ProFormaAnnotation.to_diann()`
- `ProFormaAnnotation.to_casanovo()`

These methods are currently public but raise `NotImplementedError`, while the
corresponding `from_*` paths exist. Each implementation needs a documented
supported subset and must reject annotations the target syntax cannot encode.
Add top-level batch helpers matching the existing `convert_*_sequence`
functions only if there is demonstrated batch use.

This work is logically separate from optional package integration and may ship
in its own pull request.

## Phase 1 integrations

### Pyteomics

Public functions:

```python
def to_pyteomics(annotation: ProFormaAnnotation) -> pyteomics.proforma.ProForma
def from_pyteomics(value: pyteomics.proforma.ProForma) -> ProFormaAnnotation
def to_pyteomics_composition(
    composition: Mapping[str, int] | ChargedFormula,
) -> pyteomics.mass.Composition
def from_pyteomics_composition(
    composition: Mapping[str, int],
) -> ChargedFormula
```

Implementation approach:

- Convert annotations with `pyteomics.proforma.ProForma.parse(annotation.serialize())`.
- Convert back with `ProFormaAnnotation.parse(str(value))`.
- Validate feature support before construction so unsupported cross-links or
  other constructs produce a Peptacular error rather than a low-level parser
  traceback.
- Preserve charge and adduct information when supported by both models.
- Keep composition conversion structural because both sides expose elemental
  counts; define electron/proton and isotope-key handling explicitly.

Pyteomics is already used for development cross-validation. Move to a released
version that supports the required ProForma 2.1 behavior before exposing the
adapter; do not depend on the current Git-only source override in a published
extra.

### psm_utils

Public functions:

```python
def to_psm_utils(annotation: ProFormaAnnotation) -> psm_utils.Peptidoform
def from_psm_utils(value: psm_utils.Peptidoform) -> ProFormaAnnotation
```

Implementation approach:

- Construct `Peptidoform(annotation.serialize())`.
- Convert back from its public `proforma` property.
- Document that PSM metadata is not part of the annotation conversion.
- Add one cookbook example showing a `PSMList` read, conversion to Peptacular,
  manipulation, conversion back, and write.
- Mirror the Pyteomics feature restrictions because `psm_utils` uses
  Pyteomics for its peptidoform representation.

This adapter is higher value than adding individual readers for MaxQuant,
mzIdentML, idXML, or similar formats because `psm_utils` already normalizes
those formats into one PSM model.

### AlphaBase

AlphaBase stores peptide information in pandas DataFrames rather than a native
single-peptide object. Match that public model directly and provide a row
mapping only as a lightweight single-item convenience:

```python
def to_alphabase_row(
    annotation: ProFormaAnnotation,
    *,
    loss_policy: LossPolicy = LossPolicy.ERROR,
) -> dict[str, object]

def from_alphabase_row(row: Mapping[str, object]) -> ProFormaAnnotation

def to_alphabase_dataframe(
    annotations: Iterable[ProFormaAnnotation],
    *,
    loss_policy: LossPolicy = LossPolicy.ERROR,
) -> pandas.DataFrame

def from_alphabase_dataframe(
    dataframe: pandas.DataFrame,
) -> list[ProFormaAnnotation]
```

Mapping rules:

- AlphaBase site `0` maps to the peptide N-terminus.
- Site `-1` maps to the peptide C-terminus.
- Sites `1..n` map to zero-indexed Peptacular residues.
- Multiple modifications retain order and may share a site.
- Prefer UniMod accessions or names. Formula-only and mass-only modifications
  require an explicit policy if AlphaBase has no matching modification.
- Reject global, labile, unlocalized, ambiguous, interval, isotope, cross-link,
  branched, chimeric, and non-proton-adduct features by default.
- Protein-terminal modification semantics cannot be inferred from a peptide
  alone; accept optional `is_protein_nterm` and `is_protein_cterm` flags if that
  distinction is added later.

The DataFrame adapter should call AlphaBase's public `refine_precursor_df()` and
return a value that can be assigned directly to `SpecLibBase.precursor_df`.
Pandas remains optional because this module is installed through the AlphaBase
extra. Do not create a Peptacular-owned class that implies it is an AlphaBase
domain object.

## Phase 2 integrations

### pyOpenMS

Target the native `Peptidoform` and `PeptidoformIon` Python bindings once they
are available in a stable OpenMS release. They offer canonical/lossless
ProForma serialization and explicit conversion policies. Do not make legacy
`AASequence` the primary target because it cannot represent the full
Peptacular model.

If users need `AASequence`, provide separately named adapters:

```python
to_pyopenms_aasequence(annotation, *, loss_policy=LossPolicy.ERROR)
from_pyopenms_aasequence(value, *, charge=None)
```

Keep `pyopenms` in its own extra because it is a comparatively large binary
dependency.

### rustyms

Begin with the documented Python `LinearPeptide` type. Gate conversion on a
single linear peptidoform and reject compound, cross-linked, branched, or
chimeric structures. Revisit the richer Rust model when equivalent Python
bindings are stable.

### spectrum_utils

Provide an outbound helper around `spectrum_utils.proforma.parse()` only if
benchmarks or user feedback show value from caching its parsed `Proteoform`
objects. Normal spectrum annotation already accepts a ProForma string, and the
parsed model has no documented public reverse serializer. Do not add an
artificial inbound adapter that relies on private fields.

## Packaging

Add independently installable extras rather than one mandatory interoperability
dependency set:

```toml
[project.optional-dependencies]
pyteomics = ["pyteomics>=<verified-version>"]
psm-utils = ["psm-utils>=<verified-version>"]
alphabase = ["alphabase>=<verified-version>"]
pyopenms = ["pyopenms>=<verified-version>"]
rustyms = ["rustyms>=<verified-version>"]
spectrum-utils = ["spectrum-utils>=<verified-version>"]
interop = [
    "peptacular[pyteomics,psm-utils,alphabase]",
]
```

Verify whether self-referencing extras are accepted by the selected build and
package tooling. If not, repeat the phase 1 dependency strings in `interop`.
Dependency lower bounds must be set from tested public APIs, not current latest
versions. Avoid upper bounds unless an actual incompatibility is known.

## Test strategy

Keep core tests runnable with no extras. Integration tests should use
`pytest.importorskip()` and be grouped with markers such as
`@pytest.mark.interop_pyteomics`.

Use a shared feature matrix:

| Feature | Pyteomics | psm_utils | AlphaBase | pyOpenMS | rustyms | spectrum_utils |
|---|---:|---:|---:|---:|---:|---:|
| Plain sequence | round trip | round trip | round trip | round trip | round trip | outbound |
| Localized name/accession mod | round trip | round trip | round trip | round trip | round trip | outbound |
| N/C-terminal mod | round trip | round trip | round trip | round trip | verify | outbound |
| Multiple mods at one site | verify | verify | round trip | verify | verify | outbound |
| Charge | round trip | round trip | round trip | round trip | round trip | outbound |
| Non-proton adduct | verify | verify | reject | verify | verify | verify |
| Fixed/global mod | verify | verify | reject | verify | verify | verify |
| Ambiguity/interval | verify | verify | reject | verify | verify | verify |
| Glycan | verify | verify | reject | verify | verify | verify |
| Cross-link/branch/chimera | reject if unsupported | reject if unsupported | reject | verify | reject initially | verify outbound |

For every supported row, test both directions and compare:

1. stripped sequence;
2. modification locations and identities;
3. charge and adducts;
4. monoisotopic mass when resolvable; and
5. serialized canonical ProForma where canonicalization rules agree.

For every unsupported row, assert a specific `InteropConversionError` and
message. Add dependency-missing tests in isolated subprocesses or by masking
imports, ensuring the error recommends the correct installation extra.

CI should have one lightweight phase 1 interoperability job plus separate
manual or scheduled jobs for large binary integrations such as pyOpenMS.

## Documentation

- Add an "Interoperability" page with an installation table and one round-trip
  example per integration.
- Clearly label lossless, conditionally lossless, and lossy adapters.
- Link from the existing converter example rather than mixing package-object
  conversion with vendor/search-engine string formats.
- Add a compatibility table generated or checked from the same fixtures used
  by tests so documentation cannot drift from behavior.

## Delivery sequence

1. PR 1: shared interop errors/policy and Pyteomics annotation/composition
   adapters.
2. PR 2: `psm_utils` adapter and PSM I/O cookbook.
3. PR 3: AlphaBase row/DataFrame adapters, loss reporting, and DataFrame
   example.
4. PR 4: implement the existing IP2, DIA-NN, and Casanovo outbound stubs, or
   move this PR earlier if format symmetry is a release priority.
5. PR 5+: add pyOpenMS and rustyms after API/version verification; add
   `spectrum_utils` only with a demonstrated parsed-object use case.

Each PR should update the feature matrix, optional-dependency metadata,
documentation, and tests for only its own integration.

## Decisions required before implementation

1. Whether third-party adapters belong under `peptacular.interop` (recommended)
   or as methods on `ProFormaAnnotation`.
2. Whether the existing outbound format stubs block the first integration
   release.
3. Whether warning/drop loss policies are desirable, or whether Peptacular
   should support only strict conversion initially (strict-only is safer for
   the first release).
4. Minimum supported versions after testing against released packages on
   Python 3.12.
5. Whether composition adapters are part of the first Pyteomics PR or a
   follow-up.

## Definition of done for phase 1

- Core Peptacular installs and imports without any integration dependency.
- Each missing dependency produces a targeted installation message.
- Pyteomics and `psm_utils` round-trip every mutually supported fixture without
  semantic loss.
- AlphaBase round-trips its documented subset and rejects every unsupported
  feature by default.
- Optional dependency metadata installs successfully on all supported Python
  platforms.
- CI runs core tests without extras and the phase 1 adapter suite with extras.
- Public documentation states supported features and loss behavior.

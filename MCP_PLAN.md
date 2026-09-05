# Peptacular MCP design for 3.3.0

The optional MCP integration gives agents typed access to Peptacular's scientific
API without requiring them to write Python scripts. Its main workload is a
single annotation or a small batch. Each call supplies its inputs and receives
its results directly.

This is the current release design. The earlier persistent workspace and custom
execution design has been removed. Calculation speed and the intended workload
do not justify maintaining a database, result handles, or a job manager.

## Scope

| Tool | Purpose |
| --- | --- |
| `get_reference` | Discover capabilities, installed extras, limits, enzymes, scales, notation, and schemas |
| `inspect_peptides` | Explain ProForma or stable JSON annotations without requiring mass calculation |
| `analyze_peptides` | Calculate selected precursor masses, m/z, composition, and sequence properties |
| `fragment_peptides` | Generate bounded theoretical backbone and precursor ions |
| `compare_peptides` | Compare annotations and theoretical values against an explicit reference |
| `isotope_envelopes` | Calculate approximate theoretical isotope distributions with explicit axes |
| `digest_proteins` | Digest inline protein sequences and retain source identifiers and spans |
| `find_modifications` | Find reference candidates by accession, name, or mass tolerance |
| `edit_peptides` | Apply explicit edits atomically per input |
| `enumerate_modifications` | Enumerate bounded structural candidates from residue or terminal rules |
| `map_peptides` | Find repeated and overlapping peptide locations in inline proteins |
| `convert_annotations` | Produce ProForma, stable JSON, or supported optional integration formats |

Spectacular owns observed spectra, matching, scoring, and experimental mass
errors. Cross-link calculations remain on their separate branch. This release
does not add calculation provenance infrastructure or arbitrary code execution.

## Request and response contract

Every tool takes a `request` object. Scientific `inputs` are a list of records
with an `annotation` and optional caller `id`. An annotation is ProForma text or
versioned Peptacular JSON. Mapping proteins use the same record list. Comparison
uses one explicit reference record.

Typed responses carry records, diagnostics, applied settings, a request ID,
returned row counts, and `computation.complete`. Text and structured content
contain the same result. Per-field failures retain successful measurements.
A partial status does not necessarily mean truncation, so clients inspect
completeness and diagnostics separately.

Caller IDs and source indexes associate outputs with inputs. Source keys and row
keys are local to the call. Agents pass actual returned annotations to later
calls and maintain any earlier protein or span association themselves.
Reference lookup supports simple offset paging. Scientific calculations have no
stored continuation or hidden result pages.

## Execution and limits

Scientific adapters invoke the core API directly. The MCP SDK runs synchronous
tools in ordinary threads. There is no custom process pool, scheduler, job API,
progress state, timeout manager, preflight mode, or idempotency cache. Active
threads are not forcibly terminated on client cancellation.

The limits protect against accidental expansion and oversized model responses:

- 100 records per input list, 10,000 characters per text annotation, and a
  combined 100,000-character serialized annotation budget.
- A 1 MiB request budget, including mapping proteins and comparison references.
- 1,000 output rows by default, configurable up to 5,000, with 240,000 bytes of
  serialized scientific row data. Response metadata and the duplicate text
  representation are additional.
- 50,000 eager fragment combinations per charge and 1,000 residues for isotope
  calculations.
- 200 residues and 1,000 candidates per peptide for modification enumeration,
  with at most five variable modifications.

Rows are accumulated only within these output bounds. Reaching a computation
limit sets completeness to false and reports the stop reason. The agent can
narrow the requested output or split a batch. Whole-proteome processing belongs
in the existing streaming and batch Python APIs.

## Scientific behavior

Preserve distinct neutral mass, ion mass, and m/z fields. Charge settings select
external carriers, and intrinsic charge remains explicit in the total. Conflicts
with encoded charge require an explicit override. Require nonzero charge for
m/z. Keep formula changes distinct from numeric mass shifts.

Coordinates are zero-based and end-exclusive. Modification searches return
candidates without assigning identity. Enumeration generates structures without
probabilities. Isotope results disclose their approximation, normalization, and
truncation settings. Optional conversions reject unintended information loss
and report explicitly permitted loss.

## Package structure

- `contracts.py` defines strict request and response envelopes.
- `outputs.py` defines typed scientific result rows.
- `operations.py` adapts the core scientific API and enforces expansion limits.
- `server.py` registers the 12 tools and static reference resources.
- `cli.py` provides stdio startup and the installation check.

Install the optional extra with `pip install -e ".[mcp]"` from this checkout.
Run `peptacular-mcp --check`, then configure the client to launch
`peptacular-mcp`. There are no workspace roots, caches, databases, file export
options, or cleanup commands. Ordinary package imports do not load the MCP SDK.
The [MCP guide](docs/mcp.rst) contains setup instructions and request examples.

## Validation

Automated tests compare scientific adapters with the core API and cover mixed
successes, charge conventions, optional conversions, strict schemas, bounded
results, direct result reuse, side-effect-free imports, and real stdio clients.
A dependency that prints to stdout must not corrupt protocol traffic.
Installed-wheel smoke tests verify the optional entry point and tool discovery.

Manual agent scenarios check useful workflows and whether agents understand
limits and the Spectacular boundary. They are evaluation fixtures, not claims
of measured model success. Future interface additions should follow actual
client usage rather than anticipated infrastructure needs.

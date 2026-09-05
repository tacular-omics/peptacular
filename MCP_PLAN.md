# Peptacular local MCP implementation plan

Status: implemented for the upcoming 3.3.0 release. The detailed sections below preserve the design target. See `docs/mcp.rst` for the implemented interface, tested limits, and operational tradeoffs. Claude Code and similar agents can call Peptacular tools without writing Python scripts.

The local MCP integration is included in the 3.3.0 release scope. Hosted accounts, HTTP deployment, and cloud infrastructure are outside this implementation. Cross-link calculations remain on their existing branch. Spectacular owns spectrum handling and spectra matching, so Peptacular must not introduce a competing implementation.

## 1. Intended experience

The user installs an optional dependency extra and registers a local executable with their agent. The agent discovers Peptacular tools and invokes them with structured arguments. Python runs inside Peptacular's server process and workers, without the agent writing imports, temporary scripts, or parsing terminal output.

The server should support complete tasks such as:

- Explain a ProForma annotation, identify its modifications, and calculate requested precursor properties.
- Compare two possible modification assignments and report their annotation and calculated precursor differences.
- Calculate isotope envelopes and show why a requested composition-dependent operation is unavailable.
- Find known modifications near a specified mass difference without choosing one prematurely.
- Digest a FASTA file, calculate precursors at selected charge states, filter the results, and export a table.
- Locate peptides in proteins and retain all repeated or ambiguous mappings.
- Apply explicit annotation edits or enumerate constrained modification placements.
- Check compatibility with supported proteomics packages and export the appropriate representation.

The server contains deterministic scientific operations and result management. It does not invoke another model. It does not require a model provider API key. Client-specific prompt files can be optional examples rather than a requirement for using the tools.

### Installation and connection

Proposed package interface:

```text
pip install "peptacular[mcp]"
peptacular-mcp --workspace /absolute/project/path
peptacular-mcp --check --workspace /absolute/project/path
```

These commands are implemented in the feature checkout, which is not yet published. `--check` should validate configuration, report installed capabilities and versions, and exit without starting the protocol loop or changing client configuration.

A documented Claude Code registration would have this shape once the executable exists:

```text
claude mcp add --transport stdio --scope local peptacular -- peptacular-mcp --workspace /absolute/project/path
```

Claude Code supports local stdio servers launched as subprocesses. Its registration syntax separates client options from the executable with `--`. Documentation should also show a configuration-file example and how to use an absolute executable path when PATH differs between environments. [Claude Code MCP documentation](https://code.claude.com/docs/en/mcp#option-3-add-a-local-stdio-server).

Provide a version-pinned `uvx` example when an MCP-enabled package version is actually published. Do not document a nonexistent version or install a moving Git branch by default.

## 2. What the tools should cover

The reviewed local interface has 18 tools: 12 for discovery and scientific work, and 6 for reusable data and longer jobs. This is an implementation target, not a quota. Agent evaluations can combine or rename tools when evidence supports doing so.

### Discovery and scientific work

| Tool | Principal arguments | Results and purpose |
| --- | --- | --- |
| `get_reference` | Topic, identifiers or search, cursor | Capabilities, installed extras and versions, limits, enzyme definitions, ion conventions, property scales, notation guidance, and schema references. Defaults to a concise capability summary |
| `inspect_peptides` | Inputs, requested detail | Canonical notation, length, residues, names, modifications, ambiguity, charge encoding, and structural diagnostics |
| `analyze_peptides` | Inputs, requested measurements, charge policy, mass mode, property settings | Neutral mass, ion mass, m/z, composition, selected sequence properties, and per-measurement errors |
| `fragment_peptides` | Inputs, ion series, charges, losses, isotope settings, filters, projection | Typed theoretical fragment rows, labels, spans, applied settings, and optional compositions |
| `compare_peptides` | Explicit pairs or a named reference, requested annotation and precursor comparisons | Sequence, modification, charge, composition, and calculated precursor differences. No spectrum or ion matching |
| `isotope_envelopes` | Inputs, mass or m/z axis, charge settings, resolution, isotope and abundance limits | Peak positions, abundance values, normalization convention, and completeness information |
| `digest_proteins` | Protein inputs, enzyme ID, missed cleavages, specificity, length bounds | Peptides with source records, spans, cleavage information, and per-protein diagnostics |
| `find_modifications` | Accession, name, or mass query, vocabulary selection, explicit tolerance | Matching reference entries, available specificity, mass/composition data, source IDs, and mass errors |
| `edit_peptides` | Inputs and an ordered typed edit list | New annotations, exact applied edits, differences, and failures. Source inputs remain unchanged |
| `enumerate_modifications` | Inputs, static and variable rules, placement constraints, maximum modifications and candidates | Candidate annotations with placement records, counts, and explicit expansion limits |
| `map_peptides` | Peptide inputs, protein inputs, matching policy | Every match, spans, unmapped records, coverage summaries, and ambiguous mappings |
| `convert_annotations` | Inputs, target representation, explicit loss policy | ProForma, stable JSON, supported rows or interchange text, plus compatibility and loss diagnostics |

Every scientific tool should accept small inline batches. Single-record use is a one-record batch, with one stable output shape. Two-input tools should require explicit pairing or a reference selection instead of silently forming a Cartesian product.

A separate diagnostics tool is unnecessary initially. Inspection reports structural problems, and analysis reports the requested operation's actual failures. This avoids making the agent run every expensive calculation twice.

Capabilities belong in `get_reference` instead of a separate discovery tool. Inspection remains separate from analysis because it should work without triggering mass calculations or modification lookup failures.

Properties belong in `analyze_peptides` with an explicit measurement selection. Do not create individual tools for every hydrophobicity scale or residue statistic. Similarly, use one reference tool with a finite topic enum rather than many tiny lookup tools.

### Reusable data and execution

| Tool | Principal arguments | Results and purpose |
| --- | --- | --- |
| `register_dataset` | Inline records or local path, declared format, column mapping where needed | Validated immutable dataset handle, record counts, schema, and ingestion diagnostics |
| `list_workspace` | Object kind, name/status filters, cursor | Registered datasets, stored results, exports, and known jobs with stable IDs, source relationships, and expiry. Lists managed objects rather than browsing arbitrary files |
| `get_job` | Job ID | State, progress counters, elapsed time, errors, and available result handles |
| `cancel_job` | Job ID | Cancellation request status followed by a verifiable terminal state |
| `query_result` | Result ID, projection, typed filters, sorting, aggregation, cursor | Bounded pages or summaries from stored results |
| `export_result` | Result ID, format, optional destination relative to output root, overwrite flag | Managed export path, resource URI, record counts, and format details |

These supporting tools let an agent perform useful file-based work without generating helper scripts. They also avoid sending an entire proteome or fragment library into the conversation.

## 3. Important operation semantics

### Inspect and analyze

Inspection should provide concise, structured explanations of notation. It should distinguish an accession, named modification, formula, numeric shift, ambiguous location, and cross-link notation. Parsing something does not imply every calculation is supported.

Analysis should evaluate only requested measurements. The default can return sequence length, modification summary, and neutral monoisotopic mass. Expensive isotope and fragmentation work remains explicit.

For a peptide with a numeric mass shift, a requested mass calculation can succeed while composition fails. Return both outcomes. Never discard a successful measurement because another requested measurement is unavailable.

Return theoretical precursor values for declared charge states. Comparing measured precursor values with predictions belongs in Spectacular, rather than adding observed-data inputs or experimental mass-error analysis here.

An unknown modification should preserve the accession or supplied name and provide actionable diagnostics. Reference lookup can suggest candidates, but the server does not silently rewrite the annotation to one candidate.

### Fragmentation and comparison

Use the ordinary fragment engine as the correctness reference. The fast path is an internal optimization available only when its semantics and output agree with the requested operation. The agent does not need separate fast and slow fragmentation tools.

Return fields such as source record ID, ion series, ordinal or span, signed charge, external charge, m/z, ion mass, neutral mass, loss description, and isotope adjustment. Composition and fragment sequence are optional projections because they can require additional work. If a requested label format cannot represent a mass-only loss, retain the numerical fragment and report the label limitation.

Initially advertise tested backbone series and precursor outputs. Unsupported or insufficiently reviewed residue-specific and cross-link fragmentation should not appear as supported simply because a Python enum lists it.

Comparison has explicit modes: annotation structure, composition where available, and calculated precursor values. Do not invent a sequence alignment when two different sequences are supplied. Position-wise differences apply directly to equal sequences, and a separate alignment algorithm would require its own contract and tests.

Peptacular can calculate the theoretical ions for either annotation. Spectrum matching, tolerance-based peak/ion correspondence, and interpretation of distinguishing spectral evidence belong in Spectacular. Keep `compare_peptides` useful for annotation and calculated-value differences without implementing a matching engine.

### Boundary with Spectacular

The package ownership rule is: Peptacular describes sequences and calculates their theoretical properties. Spectacular handles spectrum data and spectra matching.

| Responsibility | Owner |
| --- | --- |
| ProForma parsing, modifications, sequence transformations, and annotation comparison | Peptacular |
| Theoretical precursor values, elemental composition, and isotope envelopes | Peptacular |
| Theoretical fragment ion generation and labels | Peptacular |
| Protein digestion and peptide-to-protein sequence mapping | Peptacular |
| Reference modification lookup, including a user-supplied mass-difference search | Peptacular and its existing tacular reference dependency |
| Spectrum objects, observed peak input, and spectral data handling | Spectacular |
| Observed-to-theoretical matching and experimental mass errors | Spectacular |
| Spectrum comparison, scoring, or spectral localization evidence | Spectacular, only where its actual implementations support them |

Theoretical isotope envelopes and fragment ion tables remain Peptacular outputs. They do not require Peptacular to own a spectrum object, measured intensities, or a spectra-matching implementation. Sequence mapping remains here because peptide-to-protein matching is a sequence operation.

Do not expose `match_fragments` in this MCP server or add a Peptacular matching module. Do not move existing Spectacular algorithms into the MCP adapters. Before implementing any cross-package bridge, inspect Spectacular's actual public interfaces and supported data formats. This plan establishes ownership from the user's package description, not a claim that every possible spectral feature is already implemented there.

### Cross-package exchange

Make theoretical outputs easy for Spectacular to consume. Preserve source record IDs, canonical ProForma, ion series, ordinal or span, signed and external charge, theoretical m/z, mass mode, loss/isotope settings, and package/contract versions. Add composition only when requested and available.

Prefer Spectacular's existing interchange representation if it can preserve these fields. Otherwise, define a small versioned exchange contract in one place and test both producers and consumers. Avoid maintaining two nearly identical serializers or a circular dependency between the packages. Peptacular's normal calculation path must remain independent of Spectacular.

For a future combined agent workflow, the agent can obtain theoretical results from Peptacular and ask a Spectacular integration to match supplied spectrum data against them. Small results can be passed as structured data. Larger results should use a supported exported file in a configured shared local directory or an explicit bridge.

Opaque result IDs are local to their owning server. A Spectacular tool cannot automatically dereference a Peptacular result ID. File/resource access and any reference transfer must be explicit and supported by the receiving integration. A new server-to-server protocol is not needed for the initial Peptacular MCP release.

### Isotopes and properties

Isotope output must state whether positions are neutral masses, ion masses, m/z values, or neutron offsets. Report the normalization rule, applied abundance threshold, and truncation information. Only report retained probability when the backend exposes enough information to calculate it correctly.

Each sequence property names its scale, aggregation, missing-residue policy, and treatment of modifications. Hydrophobicity and secondary-structure tendencies are model-derived sequence properties, not measurements of a particular modified peptide.

### Digestion and mapping

Resolve enzymes by known identifiers. Unknown enzyme names should be rejected or offered as lookup suggestions, rather than silently interpreted as regular expressions. Expert custom cleavage rules are a separately gated extension.

Expose full, semi-specific, and nonspecific digestion with mandatory expansion bounds for the latter modes. For a resulting peptide, retain the protein record key, original FASTA identifier, original record index, start, end, and missed-cleavage count. Do not infer counts from display names.

Preserve duplicate FASTA identifiers. Use an internal record key based on dataset and record ordinal so duplicate names remain distinguishable. Preserve all source relationships when identical peptides are generated from multiple proteins.

Mapping must retain overlapping and repeated matches. Exact sequence matching is the default. Any modification-insensitive or residue-equivalence policy is explicit. Do not implement an I/L-equivalence mode by pretending the core API already provides one.

### Edits and enumeration

Initially support explicit localized and terminal modification edits, fixed modification expansion, charge edits, removing selected modification types, and slicing with declared coordinates. Every edit is validated before being applied to a copy. For each record, edits are atomic: a failure does not return a partially edited annotation as a success.

Enumeration is separate from deterministic editing because its expansion and result semantics differ. Require finite candidate and modification limits. Return whether all candidates were enumerated. A candidate is a structural possibility under the supplied rules, not a probability-ranked biological claim.

### Modification lookup and conversion

Modification lookup must retain all qualifying entries under a deterministic sort, with vocabulary and accession included. Mass tolerances are required for mass searches. Define signed mass errors and how ppm is computed. Near zero, require an absolute tolerance instead of dividing by zero.

Use available tacular fields without inventing missing residue specificity or chemical composition. Exact lookup and substring/prefix search should be clearly distinguished from any later fuzzy matching.

Conversion returns portable data, not Python object handles. JSON and ProForma are always available. Optional adapters add AlphaBase row conversion and compatibility-checked representations for Pyteomics and psm_utils. List precisely what each target produces.

The AlphaBase row adapter is a good default because it preserves input associations. If native DataFrame refinement is requested, carry a source-record column through sorting and test it. Do not recover associations by matching sequence strings, which may be duplicated.

Retain the existing default rejection of lossy conversion. An explicit warning/drop policy produces an itemized loss report. Known target-parser charge changes must remain errors rather than successful conversions.

## 4. Shared contracts

### Input records

Use a root JSON object for every tool. Nested discriminated input forms select inline records, a dataset, or a result column. Accept exactly one form. A record contains a caller ID when provided, plus ProForma text or supported versioned annotation JSON.

The parser should reject unknown fields, non-finite numeric inputs, booleans supplied as counts, unsupported enum values, and contradictory settings. Do not use an unrestricted `kwargs` object or arbitrary Python function selector.

Dataset references specify the source column and expected record type. A numeric result column must not accidentally become a peptide input. Inputs and results remain immutable across calls.

### Scientific conventions

- Use explicit `neutral_mass_da`, `ion_mass_da`, and `mz` fields, rather than an ambiguous generic `mass` output.
- Preserve signed charge, and distinguish external carriers, intrinsic charge, and total charge using backend definitions.
- Use encoded charge when present. If m/z is requested without a charge, report `missing_charge`. Do not silently assume charge 2.
- When requested charge settings conflict with the input annotation, require an explicit override mode and echo the applied settings.
- Distinguish a formula delta from a numeric mass delta. Only a formula implies elemental changes.
- Machine-readable spans use zero-based, end-exclusive coordinates. Human-facing residue positions use separately named one-based fields.
- Precursor fragments have no residue ordinal.
- Preserve original notation separately from canonical notation.
- Keep full calculation precision. Display rounding is a presentation setting.
- Echo property scales, ion settings, and tolerance units.

### Result envelope

Use one envelope shape across tools:

```json
{
  "contract_version": "1.0",
  "request_id": "opaque-request-id",
  "status": "complete",
  "applied_settings": {},
  "records": [],
  "diagnostics": [],
  "result_id": null,
  "page": {
    "returned_rows": 0,
    "total_rows": 0,
    "next_cursor": null
  },
  "computation": {
    "complete": true,
    "stop_reason": null
  }
}
```

The example describes the proposed shape, not a real analysis result. Job creation returns the same envelope with a job reference and the appropriate status. Each record has its own status, source identifiers, and per-measurement values or diagnostics.

Pagination is separate from computation completeness. A first page can represent a fully completed calculation. A capped calculation can have no next page while still being incomplete. Unknown totals are null, not estimates presented as exact counts.

Use the existing diagnostic code as the core category, with service fields for record ID, operation, argument path, and a documented recovery action. Report source positions only when the parser supplies reliable positions. Do not extract them through fragile regular expressions over exception messages.

Differentiate invalid requests, expected scientific limitations, partial record failures, expired references, resource limits, cancellation, and unexpected internal failures. Unexpected errors include a request ID for logs rather than a traceback in model output.

Publish output schemas and return structured content. Include a bounded textual representation for clients that primarily consume text. An initial inline response target of approximately 4,000 tokens should be measured and tuned, with full results available by query or export. [MCP tools specification](https://modelcontextprotocol.io/specification/2026-07-28/server/tools).

Record package versions and applied settings. A separate calculation-provenance framework is not a prerequisite.

## 5. Local data and workflow behavior

### File support

Initial ingestion formats:

- FASTA and gzip-compressed FASTA.
- CSV and TSV with an explicit sequence/ProForma column mapping.
- JSONL records with a declared annotation field.
- Stable ProForma JSON documents and inline record arrays.

File registration should stream into a managed snapshot with a stable dataset ID. Record the original path and file metadata for display, but run subsequent jobs against the snapshot so changing the original file does not silently change an existing dataset.

Use configured workspace/read roots and an explicit output root. Resolve symlinks and reject paths outside those roots. A client-provided path or dataset name is data, not permission to expand the roots. The server can operate on inline inputs with file ingestion disabled.

Support ordinary file encoding and gzip settings through the current FASTA API. Bound decompressed bytes and total residues, not just compressed file size. Preserve per-record errors and distinguish a malformed file structure from a malformed sequence.

### Querying results

`query_result` supports a finite expression grammar with validated column names and operators. Initial operators should cover equality, membership, numeric ranges, null/error status, sorting, count, min, max, and grouped counts. Do not accept SQL, Python expressions, or unbounded regular expressions.

Use deterministic ordering with a stable record key as a tie breaker. Cursors bind to a result snapshot and query shape. Querying a later page must not rerun the scientific calculation.

### Tool composition instead of a workflow language

Do not implement `run_workflow` or a custom step-list language initially. Claude Code already orchestrates tool calls. The useful requirement is that tools can consume each other's result references without copying whole tables through conversation.

Example user task:

> Digest this FASTA with trypsin, allow one missed cleavage, calculate precursors at charges 2 and 3, keep peptides in a specified length and m/z range, and save a CSV with the source protein and span.

The agent calls `register_dataset`, passes the dataset ID to `digest_proteins`, passes the digestion result to `analyze_peptides`, filters through `query_result`, and calls `export_result`. It uses `get_job` only for calls that return a job handle. Each scientific operation remains a separately testable request.

Queries should return an immutable view reference when their selected rows will feed another operation or export. The view binds to its original result snapshot and validated filter. It is not just a displayed preview. Every downstream tool validates the referenced record type and required columns before starting work.

All stages preserve source relationships. The agent can inspect a stage, correct a setting, or reuse a completed result without rerunning earlier stages. Use a preflight option on expensive individual tools to explain parameters and bounds. Add a server-side workflow abstraction only if measured agent use shows repeated orchestration errors or material overhead.

### Workspace discovery and recovery

Include `list_workspace` because conversations lose context and users return to previous analyses. It should find datasets, results, views, exports, and known jobs by readable label, original filename, operation, status, and creation time. Return identifiers and lineage metadata, not entire data contents.

A lost result ID must not force the agent to re-register a FASTA or repeat a long calculation. Limit listings to the configured workspace, paginate them, and identify expired or unavailable objects where metadata remains. Live job control remains specific to the owning server instance.

### Persistence and cleanup

Start with SQLite metadata and managed result files in a standard local application cache directory, namespaced by workspace. Use immutable typed result tables with normalized query columns and serialized annotation payloads. This avoids adding pandas or a distributed database to the MCP extra.

Dataset IDs, result IDs, and job IDs are opaque and distinct. Every call supplies its references explicitly. Do not rely on hidden state such as a current peptide or active dataset.

A proposed default retention is 24 hours, configurable by the user. Cache cleanup only removes service-owned data. Exported user files are outside cache cleanup. Defer the model-facing `release_data` tool. Automatic retention and a scoped `peptacular-mcp cache clean` CLI command are sufficient initially. Cleanup must skip live inputs and must not delete original files or user exports.

After a server restart, completed retained results remain readable. Do not implement checkpoints or automatic job resumption initially. Record graceful interruptions, and report earlier jobs as unavailable when final state could not be established after a crash. Do not claim background jobs continue after the client shuts down its stdio server.

Concurrent agent calls must be independent and explicitly reference their inputs. Each local server instance owns its live jobs and worker budget. Persist completed immutable results and small manifests in the workspace cache, using transactional writes and instance namespaces. Defer shared job ownership, heartbeats, takeover, and a multi-process coordinator. One server instance must never cancel or rewrite another instance's live job state.

On client disconnect or stdin EOF, stop accepting requests, cancel owned work, close workers, flush metadata, and exit. A reconnect can retrieve completed stored results, but a process crash must not be reported as a successful cancellation or completion.

## 6. Execution model and resource budgets

Small requests execute directly within a bounded worker context. Larger work returns an ordinary job handle. Every scientific tool uses the same execution machinery and has one implementation of its operation.

Use `execution.mode` with `auto`, `inline`, and `job` options. In auto mode, preflight routes oversized work to a job before computation begins. Defer transparent promotion of an already running inline call. Inline calls still have an enforced deadline and report an explicit timeout or bounded partial result. Never duplicate work while switching execution paths.

Keep computation off the async protocol event loop. Use server-owned spawn-compatible worker processes with a global concurrency cap. Do not expose Python multiprocessing options directly to the model or allow each request to create its own full-sized pool.

Jobs have queued, running, succeeded, partially succeeded, failed, cancelled, and interrupted states. Track records consumed, generated, and failed independently. Cancellation is a request until the worker actually stops. Check between bounded work units and terminate an isolated worker if a calculation cannot cooperate within its deadline. Replace terminated workers safely.

Preflight combinatorial operations, and avoid materializing unbounded fragment or modification lists. Output pagination alone does not limit computation. If an existing backend path is eager, reject excessive requests or add a bounded backend path before exposing it for large jobs.

Initial benchmark targets:

| Budget | Proposed default |
| --- | --- |
| Inline input | 100 records and 1 MiB |
| Direct calculation | Roughly 5 seconds before job handling is needed |
| Inline result preview | 25 to 50 rows with a bounded text representation |
| Result page | At most 500 rows and 256 KiB |
| Active compute workers | Configurable, initially capped at 2 and available CPUs |
| Modification candidates | 1,000 unless the user explicitly requests a higher configured limit |
| Stored data | Configurable byte quota enforced during ingestion and execution |
| Cached result retention | 24 hours |

These are initial limits to test, not measured package performance promises. Additional per-operation limits cover sequence length, total residues, isotope expansion, fragment combinations, and nonspecific digestion.

Use idempotency keys for scientific job creation and file export. A repeated key with the same request returns the existing operation. Reusing it for different input is an error. Cache identities include the dataset snapshot, settings, output projection, and backend versions.

## 7. Codebase organization

Recommended layout:

```text
src/peptacular/
    mcp/
        __init__.py
        __main__.py
        cli.py
        server.py
        config.py
        capabilities.py
        contracts/
            common.py
            inputs.py
            analysis.py
            fragmentation.py
            transformations.py
            datasets.py
            jobs.py
        adapters/
            inspection.py
            analysis.py
            fragmentation.py
            comparison.py
            isotopes.py
            digestion.py
            modifications.py
            mapping.py
            conversion.py
        serializers.py
        execution/
            runner.py
            workers.py
            budgets.py
            jobs.py
        storage/
            datasets.py
            results.py
            queries.py
            exports.py
        resources.py
        operations.py
```

Responsibilities and boundaries:

- `cli` owns startup, configuration checks, logging setup, and graceful shutdown.
- `server` registers named MCP tools and resources. It delegates rather than calculating.
- `contracts` defines finite, typed requests and portable responses, with unknown fields forbidden.
- `adapters` call the public annotation APIs where possible. They consume iterators deliberately and never expose live Python objects.
- `serializers` projects Fragment, ElementInfo-keyed compositions, isotope values, spans, and diagnostics into explicit JSON fields.
- `execution` owns workload estimation, processes, cancellation, and limits.
- `storage` owns snapshots, stable IDs, result queries, retention, and atomic exports.
- `operations` registers the finite set of internal execution handlers. It is not an exposed arbitrary-function dispatcher or workflow language.
- `resources` provides short notation/convention documentation and references to stored artifacts.

The core package does not import MCP modules from `peptacular.__init__`. Normal `import peptacular` must not start workers, create cache directories, import the MCP SDK, or probe optional integrations.

Use the official MCP Python SDK and a tested version range. The currently documented stable SDK line is v2, so implementation should use its supported interfaces and test against the actual Claude Code client rather than copying old server examples. Declare directly imported validation/configuration packages as direct optional dependencies, even if the SDK also depends on them. [Python SDK](https://github.com/modelcontextprotocol/python-sdk).

The entry point can exist in the base wheel, but must exit with an actionable installation message when the MCP extra is absent. Users can combine extras such as `peptacular[mcp,alphabase]` when those features are implemented.

Pure operations should remain usable in unit tests without starting a transport. If a new comparison algorithm is valuable to ordinary Python users, place that calculation in a core module with tests, and have its MCP adapter delegate to it. Do not bury new scientific mathematics in transport handlers.

## 8. Discovery, resources, and agent ergonomics

Every tool description should explain when to use it, the meaning of ambiguous scientific parameters, its most important limits, and how to obtain omitted results. Include realistic schema examples rather than a large generic manual in every result.

Use domain vocabulary in names and descriptions so a client can discover a tool from requests about precursor mass, tryptic digestion, modification localization, or isotope envelopes.

Provide resources for supported notation, charge conventions, enzymes, ion series, property scales, JSON schema, and stored result artifacts. Essential information must also be reachable through tools because resource interfaces vary across clients. [MCP resources specification](https://modelcontextprotocol.io/specification/2026-07-28/server/resources).

Resource content should be local and versioned with the installed package. Reference entries and FASTA headers are returned as data. The server must not treat instructions embedded in those fields as executable workflow steps.

Measure the total discovered tool-schema size. A planning target is at most approximately 12,000 tokens for the full interface, with short descriptions and shared conventions. Do not shrink schemas by replacing useful typed fields with unrestricted strings. Configure optional capability profiles at startup if evaluation shows that a smaller advertised set materially helps a client.

Mark tool effects truthfully. Inline inspection and calculation are read-only. Dataset registration, job creation, export, and cache deletion have storage effects. Protocol hints do not replace the configured file boundaries.

Logging goes to stderr in stdio mode. Stdout is reserved for protocol traffic. Capture dependency warnings and accidental worker prints so they cannot corrupt the connection. Default logs contain request IDs, operation names, sizes, durations, and failure codes rather than complete research sequences.

## 9. Backend prerequisites and known gaps

The recently added JSON representation and diagnostics are useful foundations, but they do not serialize all calculation result types automatically.

A concrete existing issue was reproduced during planning: `batch('digest', ..., method='sequential')` returns a BatchResult containing a generator. Pickling that result raises `TypeError: cannot pickle 'generator' object`, so the process executor cannot return it as-is. Resolve this before depending on the batch API for service workers, preferably before publishing 3.3.0. Bounded iterator consumption and per-record error reporting need tests.

Do not automatically call `diagnose()` before running the same expensive calculation. Execute once and capture its value or diagnostic. Define what diagnosis means for iterator-producing operations so returning an iterator is not mistaken for validating every generated record.

Other readiness checks:

- Verify that each advertised fragment mode has a supported serializer and tested scientific behavior.
- Confirm that labels unavailable for numeric mass shifts become field diagnostics rather than destroying valid numerical results.
- Verify isotope truncation and abundance semantics before promising probability coverage.
- Preserve source IDs across optional table conversions and any native sorting.
- Keep cross-link calculations excluded until their branch is independently reviewed and integrated.
- Do not advertise existing outbound IP2, DIA-NN, or Casanovo stubs as implemented formats.

## 10. Tests and acceptance criteria

Suggested test layout:

```text
tests/test_mcp/
    test_contracts.py
    test_serializers.py
    test_analysis.py
    test_comparison.py
    test_digestion.py
    test_mapping.py
    test_modifications.py
    test_conversion.py
    test_tool_composition.py
    test_workspace.py
    test_results.py
    test_limits.py
    test_jobs.py
    test_stdio.py
    test_installation.py
    fixtures/
evals/mcp/
    tasks.jsonl
    expected/
    README.md
```

Required verification:

1. Golden calculation cases cover modifications, labels, positive/negative charge, intrinsic charge, losses, and unavailable composition. Results agree with core tests and independent reference cases where available.
2. Every returned result validates against its published schema, including mixed failures and limit/cancellation responses.
3. No supported result contains a generator, enum-keyed JSON object, NaN, infinity, or a live third-party Python object.
4. Tool errors identify the affected record and requested measurement. Successful independent measurements are retained.
5. Mapping, digestion, and conversion preserve duplicates, original indexes, source relationships, and coordinate conventions.
6. Pagination has no missing or duplicate rows, and a page request does not rerun calculations.
7. Candidate and fragment limits are enforced during computation, including eager backend paths and decompressed input limits.
8. Cancellation stops the owned worker and frees resources. Restarted clients receive accurate interrupted or unavailable job status, completed retained results remain retrievable, and one instance cannot take over another instance's live job.
9. Paths and symlinks cannot escape configured read/write roots. Cache cleanup cannot delete original input files or unrelated exports.
10. The base wheel, MCP extra alone, and MCP plus each supported integration install independently. Test Python 3.12 through 3.14 and the supported desktop operating systems.
11. Real stdio tests exercise discovery, tool calls, resources, failures, and stdout cleanliness. Include a Claude Code smoke test before declaring compatibility.
12. Optional features are advertised accurately when their packages are absent or fail to import.
13. Theoretical exports preserve the fields needed by Spectacular without adding observed-spectrum inputs or matching behavior. When a bridge is implemented, test its round trip against Spectacular's actual supported representation.
14. An agent can recover dataset/result IDs through workspace listings and pass a filtered result view directly to another tool without transferring all rows into context.

Agent evaluation tasks should include ambiguity, partial errors, file-based workflows, and recovery from a wrong argument. Compare against direct Python use and a small wrapper interface. Measure task success, numerical correctness, corrective calls, latency, and context consumed. Manual scripting should not be necessary for workflows explicitly covered by the tools.

Useful adversarial cases include a pair of different annotations with identical calculated precursor masses, repeated FASTA IDs, an unlocalized modification, an invalid charge, an expired result cursor, a huge compressed input, and a conversion target that changes annotation data.

## 11. Implementation sequence

| Stage | Deliverable | Exit condition |
| --- | --- | --- |
| A. Backend readiness and contracts | Batch digestion fix, serializers, shared schemas, operation conventions | Golden cases and serialization/error tests pass |
| B. Local connection and analysis | Optional extra, entry point, reference/discovery, inspection, precursor analysis, fragments, isotopes, annotation comparison, lookup | A real agent completes inline analytical tasks without writing Python |
| C. Transformations and interoperability | Editing, bounded enumeration, mapping, and supported conversions | No silent input mutation or representation loss, and source IDs remain stable |
| D. Files and complete workflows | Dataset snapshots, workspace discovery, digestion, result views, queries, jobs, exports, retention | An agent completes a FASTA-to-filtered-table workflow with cancellation and partial-error recovery |
| E. Release hardening and evaluation | Cross-platform installation tests, protocol tests, limits, user docs, measured agent trials | The full advertised local interface meets the acceptance criteria |

Stages are implementation milestones, not an argument for shipping an unfinished five-tool product. The intended functional release includes the useful analytical tools and file/result workflow support. Keep the current 3.3.0 release separate and select the MCP release version when this work is ready.

## 12. Deliberate exclusions and later options

Do not initially expose random peptide generation, unrestricted permutations, codon enumeration, arbitrary Python execution, shell access, or database SQL. These add interface surface without serving the main analytical workflows.

Do not add hosted authentication, cloud storage, a browser dashboard, a distributed scheduler, or an LLM inside the server for this local use case.

Potential later additions should have specific demand:

- A server-side workflow abstraction, only if actual agent use demonstrates a benefit over composing tools through result references.
- Model-facing data deletion or resumable job infrastructure, only if retention and instance-owned jobs prove insufficient.
- Local artifact plots for isotope envelopes, coverage, or fragment ladders, accompanied by the underlying values.
- Additional interchange formats once their core conversion implementations are real and tested.
- Cross-link operations after that independent feature is reviewed and merged.

Spectrum scoring, database identification, false-discovery estimation, intensity prediction, and probability-based localization are new scientific capabilities. An MCP interface should not imply they exist merely because it can generate and compare theoretical fragments.

## 13. Scope review decisions

The priority is avoiding scripts for real analytical work while keeping the service maintainable.

| Decision | Reason |
| --- | --- |
| Assign spectrum matching to Spectacular | Avoids duplicating an existing package responsibility in Peptacular or its MCP layer |
| Add `list_workspace` | Lets an agent recover handles after context loss or return to existing analyses |
| Keep `analyze_peptides` theoretical | Observed precursor matching and experimental mass errors belong in Spectacular |
| Make filtered result views reusable inputs | Enables full workflows through ordinary tool composition |
| Fold capability discovery into `get_reference` | Avoids two overlapping metadata tools |
| Defer `run_workflow` | A custom workflow language adds contracts and failure modes before a demonstrated need |
| Defer `release_data` | Retention and explicit local CLI cleanup cover the initial requirement |
| Defer job resumption, takeover, and shared scheduling | Local instance-owned jobs plus retained results are sufficient for the initial target |
| Keep inspection separate from analysis | Understanding notation should not require calculation support |
| Keep editing and enumeration separate | Deterministic edits and potentially large candidate generation need different limits and semantics |

Optional adapters remain useful when they produce portable rows or validate target compatibility. Do not expose a call whose only result would be a Python object the MCP client cannot use.

Known theoretical limitations should remain visible. Equal precursor masses or a plausible reference modification mass do not by themselves establish an experimental identification. Spectacular owns spectral matching and any supported evidence-based interpretation.

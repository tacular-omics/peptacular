Local MCP tools for AI agents
=============================

Peptacular provides an optional local MCP server with 18 tools. An agent can
inspect annotations, calculate theoretical peptide properties, digest protein
files, and export results without writing Python scripts. The server does not
use an AI model or require a provider API key.

Spectacular owns observed spectra, spectrum matching, scoring, and experimental
mass errors. Peptacular's MCP tools produce theoretical values and portable
files that a separate Spectacular integration can consume explicitly.

Install and connect
-------------------

The MCP extra is part of the current unreleased source. From a checkout that
contains this feature, install it into a virtual environment::

    pip install -e ".[mcp]"
    peptacular-mcp --workspace /absolute/project --check

For optional conversions, combine extras, for example ``.[mcp,alphabase]``.
A base installation imports no MCP SDK, starts no workers, and creates no cache.
The console entry point explains how to install the extra if it is missing.

Register the executable with Claude Code::

    claude mcp add --transport stdio --scope local peptacular -- peptacular-mcp --workspace /absolute/project

Use an absolute path to the environment's executable if the client's PATH does
not include it. An equivalent client configuration is:

.. code-block:: json

    {
      "mcpServers": {
        "peptacular": {
          "command": "/absolute/venv/bin/peptacular-mcp",
          "args": ["--workspace", "/absolute/project"]
        }
      }
    }

On Windows the executable is under the environment's ``Scripts`` directory.
This configuration launches a local stdio process. It does not configure HTTP,
open a listening port, or install anything into the client automatically.

Discover and call tools
-----------------------

Every tool takes one root object with a ``request`` field. Unknown fields and
invalid numeric values are rejected. ``get_reference`` accepts an empty request:

.. code-block:: json

    {"request": {}}

Its topics include capabilities, enzymes, ions, scales, notation, conventions,
and request schemas. Scientific tools accept inline records or references:

.. code-block:: json

    {
      "request": {
        "inputs": {
          "kind": "inline",
          "records": [
            {"id": "sample-a", "annotation": "M[Oxidation]PEPTIDE/2"},
            {"id": "sample-b", "annotation": "PEP[+15.5]TIDE/2"}
          ]
        },
        "measurements": ["neutral_mass_da", "mz", "composition"]
      }
    }

Pass that request to ``analyze_peptides``. A numeric mass shift can have a valid
mass and m/z while its composition is unavailable. The result retains both
successful measurements and field diagnostics.

The tools are:

* ``get_reference`` and ``find_modifications`` for discovery and reference lookup.
* ``inspect_peptides``, ``analyze_peptides``, ``fragment_peptides``,
  ``compare_peptides``, and ``isotope_envelopes`` for annotation and theoretical calculations.
* ``digest_proteins``, ``map_peptides``, ``edit_peptides``,
  ``enumerate_modifications``, and ``convert_annotations`` for bounded transformations.
* ``register_dataset``, ``list_workspace``, ``query_result``, and ``export_result``
  for reusable snapshots and files.
* ``get_job`` and ``cancel_job`` for longer calculations.

Structured responses publish output schemas and use contract version ``1.0``.
Each response includes a request ID, status, records, diagnostics, page metadata,
and ``computation.complete``. Scientific results also have a persistent result ID.
A bounded text representation supports clients that display text primarily.

Rows preserve caller IDs alongside internal source keys. Duplicate caller IDs
are valid. Derived rows use distinct row keys, and downstream inputs identify
the source result and row. Digestion preserves original protein keys and spans
through later calculations. Earlier computed values are not silently copied
into a new calculation.

Scientific conventions
----------------------

``neutral_mass_da``, ``ion_mass_da``, and ``mz`` are separate quantities. m/z
requires a nonzero encoded or requested charge. Requested charges that conflict
with encoded charge require ``charge_policy: "override"``. Charge output names
total, external, and intrinsic charge separately where chemistry is available.
The ``charges`` argument controls external carriers, and total ion charge also
includes intrinsic charge. Neutral mass follows the core API's removal of
external carriers. Intrinsic charge encoded in a modification remains explicit.

Machine coordinates are zero-based and end-exclusive. Fragment ordinal counts
residues from the appropriate terminus. Precursor ordinal is null. Peptide
comparison uses one explicit reference, with numeric differences defined as
input minus reference. It does not infer an alignment or compare spectra.

Isotope distributions are theoretical approximations. Choose neutral mass, ion
mass, m/z, or neutron offset explicitly. Abundances are normalized to a maximum
retained peak of one. The backend prunes during convolution and applies a relative
abundance threshold. Retained probability is unavailable, and the server does
not claim that all isotopic probability was calculated.

Sequence properties use an explicit named scale and average or sum aggregation.
They ignore modifications by an explicit contract setting. These are predictions
from residue scales, not experimentally measured properties of modified peptides.

Modification searches return candidates without assigning an identity. Mass
searches require an explicit tolerance. Da errors are reference minus query, and
ppm uses the absolute query mass as the denominator. Zero mass requires Da.
Edits are atomic per input. Enumeration reports a candidate cap as incomplete.

Only tested a/b/c/x/y/z backbone series and precursor ions are exposed. A label or
composition failure does not discard a valid numeric fragment. Formula and
numeric deltas are distinct. Cross-link calculations and unfinished outbound
annotation formats remain excluded.

Reusable file workflow
----------------------

An agent can complete a FASTA-to-precursor-table task with these calls:

1. ``register_dataset`` with ``path: "proteins.fasta.gz"`` and ``format: "fasta"``.
2. ``digest_proteins`` with ``inputs`` referencing the returned dataset ID,
   ``enzyme: "trypsin"``, desired length bounds, and missed cleavages.
3. ``analyze_peptides`` with ``inputs`` referencing the digestion result ID,
   ``charges: [2, 3]``, and ``measurements: ["mz", "length"]``.
4. ``query_result`` with finite filters and ``create_view: true``.
5. ``export_result`` with the view ID, ``format: "csv"``, and a relative destination.

A reference input looks like:

.. code-block:: json

    {
      "kind": "reference",
      "reference_id": "result_ID_RETURNED_BY_THE_SERVER",
      "column": "proforma"
    }

Column values must be annotation text or stable JSON on every selected row.
Filter out failed rows before reusing a mixed result. A view contains the entire
selection, including rows outside its displayed first page. Views and exports do
not rerun scientific calculations. Query cursors bind to a result and query
shape. Repeat the filters, projection, and sort when requesting another page.

Queries support equality, inequality, numeric comparisons, membership, null
checks, sorting, count, min, max, and grouped counts. They accept no SQL or code.
Projection is a list of top-level column names. Request an export when an entire
row exceeds the page byte limit.

Input formats are FASTA, gzip-compressed input, CSV, TSV, JSONL, and Peptacular
stable JSON. Table inputs require an annotation column, with an optional ID
column. FASTA identifiers and full headers are retained separately. Duplicate
identifiers keep distinct internal record keys. Files are snapshotted at
registration, so later source edits do not alter the registered data.

Exports support JSON, JSONL, CSV, and plain-sequence FASTA. JSON exports include
result metadata. CSV serializes nested cells as JSON and prefixes formula-like
text cells with an apostrophe for spreadsheet safety. Numeric cells keep their
numeric values. FASTA rejects annotated sequences unless a plain sequence column
is selected. User destinations require explicit overwrite, with atomic file
publication. JSON and JSONL are preferable for exact portable structured data.

Execution and limits
--------------------

``execution.mode`` accepts ``auto``, ``inline``, or ``job``. Auto routes larger
requests, isotope calculations, and enumeration to jobs before calculation.
``preflight: true`` explains the chosen mode and basic input budgets without
performing science. Jobs execute the same adapters as inline calls.

Default limits include:

* At most 100 inline records, 16 MiB ingested bytes, 50,000 registered records,
  and 2 million annotation characters in a resolved input.
  The complete tool argument object is limited to 1 MiB.
* At most two active spawned workers, 32 pending requests, five seconds for
  inline calculation, and a configurable per-job deadline up to 300 seconds.
* At most 50,000 output rows, 32 MB worker result data, and 1,000 modification
  candidates per peptide. Enumeration accepts peptides up to 200 residues.
* At most 50,000 eager fragment combinations per charge and 1,000 residues
  for an isotope calculation. Mapping is capped at 100,000 peptide/protein pairs.
* At most 500 rows and 256 KB per query page. Scientific previews target 16 KB.
* A 256 MiB workspace storage quota and 24-hour retention.

The character budget conservatively includes annotation syntax. Byte, record,
sequence, expansion, and elapsed-time limits complement each other. Increasing
the page limit does not increase a computation budget. Some bounded datasets
and result queries are materialized in memory, so this local server is intended
for bounded analyses rather than whole-proteome fragment libraries.

Job progress reports input records consumed and rows generated. Final progress
also includes failed rows. A deadline or cancellation stops the isolated worker.
A cancelled computation does not publish a success result. Completed retained
results survive a restart. Jobs from another server instance report unavailable
live state, and one instance cannot cancel another instance's work. Client
disconnect stops this instance's jobs. There is no background daemon or resume.

Scientific requests and exports accept idempotency keys. Reusing a key with a
different request is rejected. Keys are scoped to the workspace and operation
family. Keep referenced snapshots until a retry has completed.

Local files and cleanup
-----------------------

``--workspace`` defines the default input root. Add ``--read-root`` explicitly
for additional directories, and choose ``--output-root`` when exports belong
elsewhere. Paths are resolved against those roots, including symlinks. Export
destinations are relative to the output root.

The default cache is an application cache directory namespaced by a hash of the
workspace path. ``--cache`` can override it, but a cache cannot be reused for a
different workspace. SQLite stores immutable snapshots and metadata, and worker
scratch files belong to their creating instance. Completed managed exports count
toward the storage quota. Temporary worker and export staging files can require
additional disk space within their operation limits.

Clean expired service-owned objects with::

    peptacular-mcp --workspace /absolute/project cache clean

Cleanup runs at server startup as well. It never deletes source files or exports
written to user destinations. It does not take over another instance's active
jobs or remove unverified scratch files from a crashed instance. Such orphaned
scratch files may require manual cleanup after all instances have stopped.

``--check`` validates configuration and reports versions and limits without
creating a cache or starting protocol traffic. It does not prove that future
file writes will succeed. Logs use stderr. Stdout is reserved for MCP traffic.

Resources and testing
---------------------

Resources expose conventions, individual request schemas, bounded result pages,
and export metadata at ``peptacular://`` URIs. The same essential information is
available through tools for clients with limited resource support. Export
resource URIs describe files, they do not transfer an arbitrary file's contents.
A different MCP server cannot dereference these opaque workspace IDs directly.

The implementation targets official MCP Python SDK 2.1.x, with a deliberately
narrow dependency range. Tests cover SDK discovery, structured schemas, real
stdio subprocess calls, worker cancellation and deadlines, optional imports,
scientific parity, file boundaries, and multi-step workflows. The repository
also includes agent scenario fixtures for manual evaluation. Passing protocol
tests is distinct from measuring Claude Code's task success with a live model.

With SDK 2.1.1, the 18 input schemas occupy approximately 37 KB of compact JSON.
Full discovery metadata, including typed output schemas, is approximately 104 KB.
This exceeds the initial full-interface token-budget target in the design plan.
The implementation keeps explicit scientific types. Client-specific capability
profiles remain a possible follow-up if agent evaluations show context pressure.

Local MCP tools for AI agents
=============================

Peptacular provides an optional local MCP server with 12 tools. An agent can
inspect annotations, calculate theoretical peptide properties, digest protein
sequences, and transform annotations without writing Python scripts. Each call
accepts records and returns results directly. The server needs no AI model or
provider API key.

Spectacular owns observed spectra, spectrum matching, scoring, and experimental
mass errors. Peptacular produces theoretical values and portable annotations
that a separate Spectacular integration can consume explicitly.

Install and connect
-------------------

Install the MCP extra into a virtual environment::

    pip install "peptacular[mcp]"
    peptacular-mcp --check

For optional conversions, combine extras, for example ``peptacular[mcp,alphabase]``.
A base installation imports no MCP SDK. The console entry point explains how to
install the extra if it is missing.

Register the executable with Claude Code::

    claude mcp add --transport stdio --scope local peptacular -- peptacular-mcp

Use an absolute path to the environment's executable if the client's PATH does
not include it. An equivalent client configuration is:

.. code-block:: json

    {
      "mcpServers": {
        "peptacular": {
          "command": "/absolute/venv/bin/peptacular-mcp",
          "args": []
        }
      }
    }

On Windows the executable is under the environment's ``Scripts`` directory.
This launches a local stdio process. No workspace configuration is required.
``--check`` reports installed versions, tools, and limits without starting
protocol traffic. Logs use stderr, and stdout is reserved for MCP traffic.

Discover and call tools
-----------------------

Every tool takes one root object with a ``request`` field. Unknown fields and
invalid numeric values are rejected. ``get_reference`` accepts an empty request:

.. code-block:: json

    {"request": {}}

Its topics include capabilities, enzymes, ions, scales, notation, conventions,
and request schemas. Scientific inputs are lists of annotation records. Each
record has ProForma text or versioned Peptacular JSON and an optional caller ID.
For example, pass this request to ``analyze_peptides``:

.. code-block:: json

    {
      "request": {
        "inputs": [
          {"id": "sample-a", "annotation": "M[Oxidation]PEPTIDE/2"},
          {"id": "sample-b", "annotation": "PEP[+15.5]TIDE/2"}
        ],
        "measurements": ["neutral_mass_da", "mz", "composition"]
      }
    }

A numeric mass shift can have a valid mass and m/z while its composition is
unavailable. The result retains successful measurements and field diagnostics.

The tools are:

* ``get_reference`` and ``find_modifications`` for discovery and reference lookup.
* ``inspect_peptides``, ``analyze_peptides``, ``fragment_peptides``,
  ``compare_peptides``, and ``isotope_envelopes`` for annotation and theoretical calculations.
* ``digest_proteins``, ``map_peptides``, ``edit_peptides``,
  ``enumerate_modifications``, and ``convert_annotations`` for bounded transformations.

Responses include typed records, diagnostics, applied settings, a request ID,
and contract version ``1.0``. Both structured content and text contain the same
records. ``returned_rows`` counts the delivered rows. ``computation.complete``
distinguishes a finished calculation from one stopped at a limit, with
``stop_reason`` explaining an early stop. A ``partial`` status can also mean a
measurement failed while other measurements succeeded. Inspect the diagnostics
and completeness separately.

Scientific results are returned entirely within the call's output budget.
``total_rows`` is null when calculation stopped early because the full count is
unknown. Nothing is retained for later retrieval. Only reference lookups use
``offset``, ``limit``, and ``next_offset`` for paging through reference data.

Rows preserve caller IDs alongside source indexes. Duplicate caller IDs are
valid. Source keys and derived row keys are local to one call. An agent can
reuse returned ``proforma`` or ``annotation`` values as inputs to another tool.
It must keep any association with earlier rows itself.

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

A short multi-step workflow
---------------------------

To digest a protein and calculate precursor m/z, first call ``digest_proteins``:

.. code-block:: json

    {
      "request": {
        "inputs": [{"id": "protein-a", "annotation": "MPEPTIDERAAK"}],
        "enzyme": "trypsin",
        "min_length": 3,
        "max_length": 30
      }
    }

The returned rows contain peptide ProForma strings, protein identifiers, and
zero-based end-exclusive spans. Pass the selected peptide strings to
``analyze_peptides``, with distinct IDs that associate them with the earlier rows:

.. code-block:: json

    {
      "request": {
        "inputs": [
          {"id": "protein-a:0-9", "annotation": "MPEPTIDER"},
          {"id": "protein-a:9-12", "annotation": "AAK"}
        ],
        "charges": [2, 3],
        "measurements": ["mz", "length"]
      }
    }

For file ingestion, filtering, and saving tables, use the client's existing
file tools or Peptacular's Python APIs. Large FASTA workflows can use
``iter_fasta`` and ``iter_batch`` from the :doc:`streaming` guide. The MCP
interface is intended for individual annotations and small batches.

Execution and limits
--------------------

Tools call the existing scientific API synchronously using the SDK's normal
thread handling. The server has no custom worker processes, job queue, stored
results, cache, or cleanup command. Repeating a request recalculates its result.
There is no hard computation deadline, and client cancellation does not
forcibly stop an active calculation thread.

Limits include:

* At most 100 records per input list and 10,000 characters per text annotation.
  The request is limited to 1 MiB, with a combined 100,000-character budget for
  serialized annotations, including comparison references and mapping proteins.
* ``max_rows`` defaults to 1,000 and can be increased to 5,000. Scientific row
  data is limited to 240,000 serialized bytes, excluding response metadata and
  the second content representation.
* At most 50,000 eager fragment combinations per charge and 1,000 residues per
  isotope calculation.
* At most 1,000 modification candidates per peptide. Enumeration accepts
  peptides up to 200 residues and at most five variable modifications.

The character budget conservatively includes annotation syntax. Expansion
limits also constrain work that occurs before output rows become available.
If a result stops at a row, byte, candidate, or resource limit, narrow the
request or split its inputs into smaller calls. There is no hidden next page of
scientific results. Splitting a batch cannot make an oversized individual
annotation or expansion fit, so reduce its settings instead.

Resources and testing
---------------------

Resources expose conventions at ``peptacular://conventions`` and request schemas
at ``peptacular://schemas/{tool}``. The same information is available through
``get_reference`` for clients with limited resource support.

The implementation targets official MCP Python SDK 2.1.x, with a deliberately
narrow dependency range. Tests cover discovery, structured schemas, real stdio
subprocess calls, protocol-safe dependency output, optional imports, scientific
parity, request and result bounds, and direct multi-step workflows. Manual agent
scenario fixtures are separate from these automated protocol tests.

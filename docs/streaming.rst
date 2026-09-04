Streaming, batch results, and diagnostics
=========================================

Use ``iter_fasta`` to read one protein at a time and ``iter_batch`` to process a
bounded number of peptide inputs. Existing ``parse_fasta`` and sequence functions
keep their list-returning interfaces. Existing batch calls such as ``mass(list)``
continue to raise exceptions on errors.

**Streaming FASTA**

.. testcode::

   import io
   import peptacular as pt

   stream = io.StringIO(">protein_a\nPEPTIDE\n>protein_b\nMKR\n")
   for protein in pt.iter_fasta(stream):
       print(protein.header, protein.sequence)
   assert not stream.closed

.. testoutput::

   protein_a PEPTIDE
   protein_b MKR

Paths ending in ``.gz`` are decompressed automatically. For example,
``pt.iter_fasta("proteins.fasta.gz")`` yields the same records as an uncompressed
file. File encodings are sampled at the start. Pass ``encoding="latin-1"`` or
another explicit encoding when a file's later contents require it. Binary streams
default to UTF-8, with an optional byte order mark. Text streams are already decoded.

The iterator uses memory proportional to the largest protein record. Collecting
all records or all downstream results into a list still requires memory for that
list. Empty records are skipped, matching ``parse_fasta``. Invalid structure raises
``ValueError`` with a line number where applicable. An error in a later record can
occur after earlier records have already been processed.

Files opened by the iterator close when it is exhausted, closed, or raises an error.
Caller-owned streams remain open. Use ``contextlib.closing`` when stopping early:

.. code-block:: python

   from contextlib import closing

   with closing(pt.iter_fasta("proteins.fasta.gz")) as proteins:
       first = next(proteins)

Binary streams can be buffered ahead of the last returned record. Their position
after early termination is not guaranteed to point to the next FASTA header.

**Collecting errors without losing successful results**

``batch`` returns a list of ``BatchResult`` objects. Each has an input index,
original input, value, optional diagnostic, and an ``ok`` property. Results remain
in input order. ``iter_batch`` has the same behavior but yields bounded batches.

.. testcode::

   inputs = ["PEPTIDE", "PEP[UnknownModification]TIDE", "PEPTIDES"]
   results = pt.batch("mass", inputs, errors="collect")
   for result in results:
       if result.ok:
           print(result.index, round(result.value, 3))
       else:
           print(result.index, result.error.code)

.. testoutput::

   0 799.36
   1 unresolved_modification
   2 886.392

The default is ``errors="raise"``. Collection applies to expected sequence errors,
including invalid notation, unresolved modifications, and unavailable compositions.
Bad operation names, misspelled keywords, failures reading the input iterator,
unexpected programming errors, serialization failures, and crashed workers still
raise. It does not silently suppress every exception.

Supported operations are ``parse``, ``mass``, ``mz``, ``comp``, ``fragment``,
``fast_fragment``, ``digest``, and ``isotopic_distribution``. Extra keyword arguments
are passed to the corresponding annotation method. A digest result contains the
annotation method's spans. It is not the functional ``digest`` method's sequence/span
pairs. Parsing keeps its existing ``validate=False`` default, with ``validate=True``
available for full annotation validation.

.. testcode::

   charged = pt.batch("mz", ["PEPTIDE"], charge=2)
   assert charged[0].value == pt.mz("PEPTIDE", charge=2)
   assert pt.batch("parse", ["PEP[CustomName]TIDE"])[0].ok

An iterator can connect FASTA processing with calculation without accumulating the
whole database:

.. code-block:: python

   from contextlib import closing

   if __name__ == "__main__":
       with closing(pt.iter_fasta("proteins.fasta.gz")) as proteins:
           sequences = (protein.sequence for protein in proteins)
           with closing(pt.iter_batch("mass", sequences, errors="collect", batch_size=256)) as results:
               for result in results:
                   save_result(result)

``save_result`` represents your own output function. Batch memory is bounded by
``batch_size`` inputs and their results, but one result can be large, such as a
protein's fragment list. The input iterator remains owned by the caller. Close both
iterators explicitly when breaking out of processing early.

Automatic execution uses a conservative 1,000-item threshold. Smaller batches run
sequentially. Explicit ``method="process"`` or ``method="thread"`` overrides this
choice, as does setting ``n_workers``. Small ``batch_size`` values can therefore keep
a large stream sequential. Benchmark your operation before selecting a backend.
Pools are reused across batches and close on completion or failure. Closing waits
for already running work. Process execution requires a guarded script entry point
on spawn platforms. ``start_method="spawn"`` selects a local process context without
changing application-wide multiprocessing settings.

**Understanding a failure**

``diagnose`` executes one requested operation, returning ``None`` on success or a
``Diagnostic`` on an expected input error. The result contains ``code``, ``stage``,
``message``, and ``exception_type``. Stages distinguish parsing, validation, and
calculation. Parser positions and other details remain in the original message.
Diagnostics have the same computational cost as executing the requested operation.

.. testcode::

   assert pt.diagnose("PEP[+42]TIDE", "mass") is None
   error = pt.diagnose("PEP[+42]TIDE", "comp")
   print(error.code, error.stage)
   assert pt.diagnose("PEP[CustomName]TIDE", "parse") is None

.. testoutput::

   unavailable_composition calculate

Codes are ``invalid_input``, ``invalid_notation``, ``invalid_annotation``,
``unresolved_modification``, ``unresolved_reference``, ``unavailable_composition``,
``invalid_adjustment``, ``unsupported_operation``, and ``calculation_error``.
The last code preserves errors without a more specific category. ``invalid_input``
applies to the wrong input type. Human-readable messages can evolve independently
of these codes. Existing exception-based code can catch the new
``CompositionError``, ``UnknownModificationError``, ``InvalidAdjustmentError``, and
``UnsupportedOperationError`` classes, all of which inherit from ``ValueError``.

**Calculation behavior tightened in this release**

Global isotope substitutions and elemental adjustments now share the same order
across mass and composition calculations. Explicit isotope counts and elemental
losses are validated against the complete neutral ion composition, including its
terminal atoms. Mass-only deltas remain additive, without an inferred composition.
Impossible counts and non-finite or negative calculated masses raise errors.

``fast_fragment`` supports a, b, c, x, y, z, p, and n. It uses regular calculations
when isotope labels, intrinsic charges, labile modifications, negative charge, or
average masses require them. Its dictionary keys retain the requested external
proton charge, while intrinsic charge contributes to the calculated m/z. Other ion
series now raise a clear error directing callers to ``fragment``. Charges must be
nonzero integers. These changes deliberately reject cases that previously produced
misleading numbers.

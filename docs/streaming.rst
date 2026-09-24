Streaming, batch results, and diagnostics
=========================================

peptacular does not read files. Use `fastatacular <https://github.com/tacular-omics/fastatacular>`_
(``pip install fastatacular``) or pefftacular to read FASTA or PEFF, and pass the entries
straight in: every sequence function, ``batch``, ``iter_batch`` and ``diagnose`` accept any
object with a ``sequence`` string attribute (the :class:`~peptacular.sequence.util.HasSequence` protocol).
``iter_batch`` processes a bounded number of inputs at a time. Existing batch calls such as
``mass(list)`` continue to raise exceptions on errors.

**Streaming FASTA**

.. testcode::

   import io

   from fastatacular import FastaReader

   import peptacular as pt

   stream = io.StringIO(">protein_a\nPEPTIDE\n>protein_b\nMKR\n")
   with FastaReader(stream) as proteins:
       for protein in proteins:
           print(protein.identifier, protein.sequence, round(pt.mass(protein), 3))

.. testoutput::

   protein_a PEPTIDE 799.36
   protein_b MKR 433.247

``FastaReader`` reads one entry at a time and decompresses ``.gz``, ``.bz2`` and ``.xz``
paths. ``read_fasta`` returns a list. See the fastatacular documentation for header
fields, encodings and error reporting. Unlike the removed ``pt.iter_fasta``, fastatacular
keeps residue case as written and raises on an entry with no sequence; peptacular 4.x
uppercased sequences and silently skipped empty entries.

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
are passed to the corresponding annotation method. The ``digest`` operation calls
``ProFormaAnnotation.digest_spans``, so a result contains spans, not the functional
``pt.digest`` sequence/span pairs. Parsing keeps its existing ``validate=False`` default, with ``validate=True``
available for full annotation validation.

.. testcode::

   charged = pt.batch("mz", ["PEPTIDE"], charge=2)
   assert charged[0].value == pt.mz("PEPTIDE", charge=2)
   assert pt.batch("parse", ["PEP[CustomName]TIDE"])[0].ok

An iterator can connect FASTA reading with calculation without accumulating the
whole database:

.. code-block:: python

   from contextlib import closing

   from fastatacular import FastaReader

   if __name__ == "__main__":
       with FastaReader("proteins.fasta.gz") as proteins:
           with closing(pt.iter_batch("mass", proteins, errors="collect", batch_size=256)) as results:
               for result in results:
                   save_result(result.input.accession, result)

``save_result`` represents your own output function. Batch memory is bounded by
``batch_size`` inputs and their results, but one result can be large, such as a
protein's fragment list. The input iterator remains owned by the caller. Close both
iterators explicitly when breaking out of processing early.

Automatic execution uses a conservative 1,000-item threshold. Smaller batches run
sequentially. Explicit ``method="process"`` or ``method="thread"`` overrides this
choice, as does setting ``n_workers``. Small ``batch_size`` values can therefore keep
a large stream sequential. Benchmark your operation before selecting a backend.
The same rule applies to the list form of every functional call (``pt.digest``,
``pt.mass``, ``pt.fragment``, ...): a list of 1,000 or more items with neither
``n_workers`` nor ``method`` fans out to a process pool with one worker per available
CPU (``os.process_cpu_count()``), or a thread pool on free-threaded Python. On a
shared machine, pass ``n_workers`` to cap it or ``method="sequential"`` to stay in
the calling process.
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
of these codes.

**Exceptions**

Code that raises instead of collecting can catch these classes, all importable from
``peptacular`` (``pt.ProFormaFormatError`` and so on):

- ``PeptacularError``: the base class of all of them (4.2). It subclasses ``ValueError``,
  so existing ``except ValueError`` handlers keep working.
- ``ProFormaFormatError``: the string is not valid ProForma. Raised by ``parse`` and
  ``parse_chimeric``, and by calculations that parse a modification, glycan, isotope
  label or adduct lazily (``pt.mass("<113C>PEPTIDE")``).
- ``UnknownModificationError``: a modification name or accession does not resolve.
- ``CompositionError``: a composition or mass is not available (delta-mass
  modifications with a composition request, or an empty sequence).
- ``InvalidAdjustmentError``: impossible isotope or delta counts.
- ``InvalidPositionError`` (4.2): a slice index or fragment position is outside the
  sequence.
- ``UnsupportedOperationError``: the operation does not support this input, for
  example an unknown ion type (the message lists the valid ones).

.. testcode::

   try:
       pt.parse("PEP[TIDE")
   except pt.PeptacularError as e:
       print(type(e).__name__, isinstance(e, ValueError))

.. testoutput::

   ProFormaFormatError True

**Calculation behavior tightened in 3.3.0**

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

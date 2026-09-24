Tables with pandas or polars
============================

peptacular does not ship pandas or polars, not even as an optional extra. What it gives you
is rows: :func:`peptacular.digest_records` and :func:`peptacular.fragment_records` return a
list of plain ``dict`` objects, one per peptide or ion, holding only strings, numbers,
booleans and ``None``. Install whichever DataFrame library you use and pass the list to it.

Digest rows
-----------

.. testcode::

   import peptacular as pt

   rows = pt.digest_records("MKVLATSAGERTIDEK", "trypsin", missed_cleavages=1)
   print(rows[1])
   print(pt.DIGEST_RECORD_KEYS)

.. testoutput::

   {'peptide': 'MKVLATSAGER', 'stripped_sequence': 'MKVLATSAGER', 'start': 0, 'end': 11, 'missed_cleavages': 1, 'semi': False, 'accession': None}
   ('peptide', 'stripped_sequence', 'start', 'end', 'missed_cleavages', 'semi', 'accession')

- ``peptide`` is ProForma and keeps the protein's modifications; ``stripped_sequence`` is the
  residues only.
- ``start`` and ``end`` are 0-based and half-open: ``protein[start:end]``.
- ``semi`` is True when one end is not a cleavage site or a protein terminus (only with
  ``semi=True``).
- ``accession`` is copied from the input's ``accession`` attribute (a fastatacular entry), or
  else its ``db_unique_id`` (a pefftacular entry), so each entry labels its own rows. Pass a
  list of entries, or the generator a FASTA or PEFF reader returns, to get one flat table for
  a whole proteome.

Fragment rows
-------------

:func:`~peptacular.fragment_records` takes the list :func:`~peptacular.fragment` returns, or
the list of lists it returns for several peptides (the rows come out flat, and
``parent_sequence`` says which peptide each ion belongs to). The keys are named after the
:class:`~peptacular.Fragment` constructor arguments:

.. testcode::

   frags = pt.fragment("PEPTIDE/2", ion_types=("b", "y"), charges=(1, 2))
   rows = pt.fragment_records(frags)
   print(len(rows), pt.FRAGMENT_RECORD_KEYS)
   print({k: rows[1][k] for k in ("ion_type", "position", "charge_state", "sequence", "mzpaf")})

.. testoutput::

   28 ('ion_type', 'position', 'end_position', 'charge_state', 'mz', 'mass', 'neutral_mass', 'monoisotopic', 'deltas', 'isotopes', 'sequence', 'parent_sequence', 'mzpaf')
   {'ion_type': 'b', 'position': 2, 'charge_state': 1, 'sequence': 'PE', 'mzpaf': 'b2{PE}'}

- ``position`` is the ion number, or the start of an internal ion, whose end is in
  ``end_position``. Both are int or None (None for precursor ions; ``end_position`` is None
  for every ion that is not internal), so each column has one type.
- ``deltas`` holds the neutral losses and gains, ``isotopes`` the isotope swaps, both as
  strings (``""`` for none). A formula counts as a loss: water loss is ``"H-2O-1"``, a
  water gain ``"H-2O-1^-1"``; a mass is written with its sign (``"-17.0^2"``).
- ``sequence`` and ``parent_sequence`` leave out the ``/charge`` suffix; the charge is in
  ``charge_state``.
- ``mzpaf`` is the :meth:`~peptacular.Fragment.to_mzpaf` label (``"b3{PEP}-H2O"``,
  ``"b3{PEP}-34.0"``), or None for an ion type mzPAF cannot write.

pandas
------

.. testcode::
   :skipif: pd is None

   import pandas as pd

   peptides = pd.DataFrame(pt.digest_records("MKVLATSAGERTIDEK", "trypsin", missed_cleavages=1))
   print(peptides[["peptide", "start", "end", "missed_cleavages"]].to_string(index=False))

   ions = pd.DataFrame(pt.fragment_records(frags))
   print(ions.groupby("charge_state").size().to_dict())

.. testoutput::
   :skipif: pd is None

          peptide  start  end  missed_cleavages
               MK      0    2                 0
      MKVLATSAGER      0   11                 1
        VLATSAGER      2   11                 0
   VLATSAGERTIDEK      2   16                 1
            TIDEK     11   16                 0
   {1: 14, 2: 14}

polars
------

.. testcode::
   :skipif: pl is None

   import polars as pl

   ions = pl.DataFrame(pt.fragment_records(pt.fragment("PEPTIDE", ion_types=("y",), charges=(1,))))
   print(ions.filter(pl.col("position") <= 3)["mz"].round(3).to_list())

.. testoutput::
   :skipif: pl is None

   [148.06, 263.087, 376.171]


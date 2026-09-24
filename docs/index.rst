Peptacular Documentation
========================

peptacular parses, edits and analyzes peptide and protein sequences written in
**ProForma 2.1** notation. It calculates mass, m/z and elemental composition,
generates fragment ions and isotope envelopes, digests proteins, and predicts
physicochemical properties. List inputs are processed in parallel automatically.

New to peptacular? Start with :doc:`installation` and :doc:`quickstart`, then browse
:doc:`features` for a tour of what's available.

.. code-block:: python

   import peptacular as pt

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   peptide.mass()          # 851.357...
   peptide.mz()            # 425.679...

Related packages
----------------

* `tacular <https://tacular.readthedocs.io/>`_ provides the amino acid, element and
  modification (UNIMOD, PSI-MOD, RESID, XLMOD, GNOme, UniProt-PTM) lookups that
  peptacular uses.
* `paftacular <https://paftacular.readthedocs.io/>`_ parses and serializes mzPAF
  peak annotations, and can use peptacular for sequence-aware mass calculations.

User Guide
==========

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   features
   streaming
   proforma_compliance
   examples
   interoperability
   masses
   mass_calculation
   json_serialization
   mcp

API Reference
=============

.. toctree::
   :maxdepth: 2

   api

Project
=======

.. toctree::
   :maxdepth: 1

   migration
   changelog
   citation

Indices and Tables
===================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

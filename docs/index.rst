Peptacular Documentation
========================

A Python package for peptide sequence analysis built around **ProForma 2.1** notation.

Peptacular provides tools for parsing ProForma sequences, calculating masses, generating fragments,
predicting isotopic patterns, and analyzing physicochemical properties. Built with performance in mind,
it supports parallel processing and lazy loading for efficient batch operations.

New to Peptacular? Start with :doc:`quickstart`, then browse :doc:`features` for a tour of what's available.

.. code-block:: python

   import peptacular as pt

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   peptide.mass()          # 851.357...
   peptide.mz()            # 425.679...

User Guide
==========

.. toctree::
   :maxdepth: 2

   quickstart
   features
   proforma_compliance
   examples
   masses
   mass_calculation

API Reference
=============

.. toctree::
   :maxdepth: 2

   api

Indices and Tables
===================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
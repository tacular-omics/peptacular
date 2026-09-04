API Reference
=============

.. automodule:: peptacular
   :members:
   :undoc-members:
   :show-inheritance:

Core
----

Peptacular contains a functional and object-oriented API for working with peptides and proteins. Everything
can be accessed through the peptacular namespace (``import peptacular as pt``), but for clarity the API is broken down into sections below.


Sequence
~~~~~~~~

Processing
^^^^^^^^^^

.. automodule:: peptacular.sequence.digestion
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.transformations
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.mod_builder
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.subseqs
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.basic
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.combinatoric
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.converters
   :members:
   :undoc-members:
   :show-inheritance:

Mass/Comp/Isotope/Fragment
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: peptacular.sequence.isotope
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.mass_funcs
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.fragmentation
   :members:
   :undoc-members:
   :show-inheritance:

Property
^^^^^^^^

.. automodule:: peptacular.sequence.properties
   :members:
   :undoc-members:
   :show-inheritance:

Data Classes
~~~~~~~~~~~~

Some of the core data classes used throughout Peptacular.

.. autoclass:: peptacular.annotation.mod.Interval
.. autoclass:: peptacular.annotation.mod.Mod
.. autoclass:: peptacular.annotation.mod.Mods
.. autoclass:: peptacular.spans.Span
.. autoclass:: peptacular.isotope.IsotopicData
.. autoclass:: peptacular.isotope.IsotopeLookup

Optional interoperability
-------------------------

Third-party dependencies are imported only when an adapter is called. See
:doc:`interoperability` for installation and compatibility details.

.. automodule:: peptacular.interop.pyteomics
   :members:
   :show-inheritance:

.. automodule:: peptacular.interop.psm_utils
   :members:
   :show-inheritance:

.. automodule:: peptacular.interop.alphabase
   :members:
   :show-inheritance:

Data Modules
------------

Peptacular relies on the python tacular package for data on amino acids, modifications, elements, and other
chemical entities, but are accessible from the peptacular namespace.

see https://tacular.readthedocs.io

.. automodule:: peptacular.batch
   :members:

.. automodule:: peptacular.diagnostics
   :members:

.. autofunction:: peptacular.iter_fasta

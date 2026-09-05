Examples
========

This section provides practical examples demonstrating key features of Peptacular.

ProForma Notation
-----------------

Basic usage of ProForma notation for representing modified peptides.

.. literalinclude:: ../examples/proforma_notation.py
   :language: python
   :linenos:

Mass, m/z, and Composition
---------------------------

Calculate masses, m/z ratios, and elemental compositions.

.. literalinclude:: ../examples/mass_mz_comp.py
   :language: python
   :linenos:

Annotation
----------

Parse and work with ProForma annotations.

.. literalinclude:: ../examples/annotation.py
   :language: python
   :linenos:

Annotation with Modification Objects
-------------------------------------

Advanced annotation examples using modification objects.

.. literalinclude:: ../examples/annotation_mod_objects.py
   :language: python
   :linenos:

Digestion
---------

Digest protein sequences with various proteases.

.. literalinclude:: ../examples/digest.py
   :language: python
   :linenos:

FASTA Digestion
---------------

Parse protein sequences from a FASTA file and digest each one.

.. literalinclude:: ../examples/digest_fasta_file.py
   :language: python
   :linenos:

Fragmentation
-------------

Generate theoretical fragment ions for peptides.

.. literalinclude:: ../examples/fragment.py
   :language: python
   :linenos:

Isotope Calculations
--------------------

Calculate isotopic distributions for peptides.

.. literalinclude:: ../examples/isotope.py
   :language: python
   :linenos:

Physiochemical Properties
-------------------------

Calculate various properties like pI, hydrophobicity, etc.

.. literalinclude:: ../examples/physiochemical_properties.py
   :language: python
   :linenos:

Converters
----------

Convert sequences from other tools (IP2, DIANN, Casanovo, MS2PIP) to ProForma format.

.. literalinclude:: ../examples/converters.py
   :language: python
   :linenos:

Parallel Processing
--------------------

Batch-process peptides across mass, composition, and property calculations, automatically
parallelized for list inputs.

.. literalinclude:: ../examples/parallel.py
   :language: python
   :linenos:


Streaming and Diagnostics
-------------------------

.. literalinclude:: ../examples/streaming.py
   :language: python
   :linenos:

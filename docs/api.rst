API Reference
=============

.. automodule:: peptacular
   :members:
   :undoc-members:
   :show-inheritance:
   :ignore-module-all:

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

Localization
^^^^^^^^^^^^

.. automodule:: peptacular.sequence.localization
   :members:
   :show-inheritance:

Records
^^^^^^^

.. automodule:: peptacular.sequence.records
   :members:
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

Annotation classes
~~~~~~~~~~~~~~~~~~

:class:`~peptacular.annotation.ProFormaAnnotation` is the main object returned by
``pt.parse``. Its ``prop`` attribute is an
:class:`~peptacular.property.AnnotationProperties` (see `Properties and scales`_).

.. automodule:: peptacular.annotation
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: AnnotationProperties, ProFormaAnnotation

.. autoclass:: peptacular.annotation.ProFormaAnnotation
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

ProForma components
~~~~~~~~~~~~~~~~~~~

The structured component model behind parsing and ProForma JSON.

.. automodule:: peptacular.proforma_components
   :members:
   :show-inheritance:

Chemistry
~~~~~~~~~

.. automodule:: peptacular.chem
   :members:

Isotope distributions
~~~~~~~~~~~~~~~~~~~~~

Low-level BRAIN isotope functions. The sequence-level wrappers are in
:mod:`peptacular.sequence.isotope` above.

.. automodule:: peptacular.isotope
   :members:
   :show-inheritance:

Digestion and spans
~~~~~~~~~~~~~~~~~~~

.. automodule:: peptacular.digestion
   :members:
   :show-inheritance:

.. automodule:: peptacular.spans
   :members:
   :show-inheritance:

Input protocol
~~~~~~~~~~~~~~

Sequence functions accept a ProForma ``str``, a :class:`~peptacular.ProFormaAnnotation`,
or any object with a ``sequence`` string attribute, such as a fastatacular
``SequenceEntry``. peptacular itself does not read FASTA files.

.. autoclass:: peptacular.sequence.util.HasSequence

Properties and scales
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: peptacular.property
   :members:
   :show-inheritance:

Constants and utilities
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: peptacular.constants
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.utils
   :members:
   :undoc-members:

.. automodule:: peptacular.sequence.parallel
   :members:

Constants, type aliases and data tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are available from the top-level ``peptacular`` namespace, except the
``proforma_components`` type aliases, which are importable from that module.

**Physical constants** (Da)

.. py:data:: peptacular.constants.PROTON_MASS
   :type: float
   :value: 1.007276466621

   Proton mass (CODATA 2018). Re-exported from ``tacular.constants``.

.. py:data:: peptacular.constants.HYDROGEN_BINDING_MASS
   :type: float
   :value: 1.4300064909988919e-08

   ``PROTON_MASS - (HYDROGEN_MASS - ELECTRON_MASS)``, the hydrogen 1s binding energy as a
   mass. A monoisotopic proton charge adds ``PROTON_MASS``; a mass summed from an ion's
   composition (one H atom per proton) adds this term per proton. See :doc:`mass_calculation`.

.. py:data:: peptacular.constants.ELECTRON_MASS
   :type: float
   :value: 0.000548579909065

   Electron mass (CODATA 2018). Re-exported from ``tacular.constants``.

.. py:data:: peptacular.constants.NEUTRON_MASS
   :type: float
   :value: 1.00866491595

   Neutron mass (CODATA 2018). Re-exported from ``tacular.constants``.

.. py:data:: peptacular.constants.C13_NEUTRON_MASS
   :type: float
   :value: 1.00335483507

   Mass difference between carbon-13 and carbon-12 (``tacular.constants.C13_C12_MASS_DIFF``).

.. py:data:: peptacular.constants.PEPTIDE_AVERAGINE_NEUTRON_MASS
   :type: float
   :value: 1.002856

   Average isotope spacing used for peptide averagine estimates.

.. py:data:: peptacular.isotope.AVERAGINE_RATIOS
   :type: dict[ElementInfo, float]

   Elemental ratios (C, H, N, O, S per Da) of the peptide averagine model.

**Literal and union type aliases**

.. py:data:: peptacular.constants.ModTypeLiteral

   ``Literal['nterm', 'cterm', 'isotope', 'static', 'labile', 'unknown', 'interval', 'internal', 'charge']``

.. py:data:: peptacular.constants.ParallelMethodLiteral

   ``Literal['process', 'thread', 'sequential']``

.. py:data:: peptacular.batch.BatchOperation

   ``Literal['parse', 'mass', 'mz', 'comp', 'fragment', 'fast_fragment', 'digest', 'isotopic_distribution']``

.. py:data:: peptacular.property.AggregationMethodLiteral

   ``Literal['sum', 'avg']``

.. py:data:: peptacular.property.MissingAAHandlingLiteral

   ``Literal['zero', 'avg', 'min', 'max', 'median', 'error', 'skip']``

.. py:data:: peptacular.property.WeightingMethodsLiteral

   ``Literal['uniform', 'linear', 'exponential', 'gaussian', 'sigmoid', 'cosine', 'sinusoidal']``

.. py:data:: peptacular.proforma_components.SEQUENCE_TYPE

   ``SequenceElement | SequenceRegion``.

.. py:data:: peptacular.proforma_components.GLOBAL_CHARGE_TYPE

   ``int | tuple[GlobalChargeCarrier, ...]``.

.. py:data:: peptacular.proforma_components.MODIFICATION_AMBIGUOUS_TYPE

   ``ModificationAmbiguousPrimary | ModificationAmbiguousSecondary``.

.. py:data:: peptacular.proforma_components.MODIFICATION_TYPE

   Any modification component: ``MODIFICATION_AMBIGUOUS_TYPE | ModificationCrossLinker | ModificationTags``.

.. py:data:: peptacular.proforma_components.MODIFICATION_TAG_TYPE

   Union of the modification tag component classes (``TagAccession``, ``ChargedFormula``, ``GlycanTag``, ``TagInfo``, ``TagMass``, ``TagName``, ``TagCustom``, ``PositionScore``, ``PositionTag``, ``LimitTag``, ``ComkpTag``, ``ComupTag``).

**Property scale tables**

Each maps a scale enum member to its per-residue value table. See :class:`~peptacular.property.PropertyScale` and the scale enums above.

.. py:data:: peptacular.property.PROPERTY_SCALES
   :type: dict[str, dict[str, float]]

   All scales.

.. py:data:: peptacular.property.HYDROPHOBICITY_SCALES
   :type: dict[str, dict[str, float]]

   Hydrophobicity scales.

.. py:data:: peptacular.property.HYDROPHILICITY_SCALES
   :type: dict[str, dict[str, float]]

   Hydrophilicity scales.

.. py:data:: peptacular.property.SURFACE_ACCESSIBILITY_SCALES
   :type: dict[str, dict[str, float]]

   Surface accessibility scales.

.. py:data:: peptacular.property.HPLC_SCALES
   :type: dict[str, dict[str, float]]

   Hplc scales.

.. py:data:: peptacular.property.FLEXIBILITY_SCALES
   :type: dict[str, dict[str, float]]

   Flexibility scales.

.. py:data:: peptacular.property.POLARITY_SCALES
   :type: dict[str, dict[str, float]]

   Polarity scales.

.. py:data:: peptacular.property.COMPOSITION_SCALES
   :type: dict[str, dict[str, float]]

   Composition scales.

.. py:data:: peptacular.property.PHYSICAL_PROPERTY_SCALES
   :type: dict[str, dict[str, float]]

   Physical property scales.

**ProForma JSON**

.. py:data:: peptacular.proforma_json.PROFORMA_JSON_SCHEMA_VERSION
   :type: str
   :value: '1.0'

   Current ProForma JSON schema version.

.. py:data:: peptacular.proforma_json.PROFORMA_JSON_SCHEMA_ID
   :type: str

   ``$schema`` URL written into every document. The schema file is published at that address.

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

Batch processing and diagnostics
--------------------------------

See :doc:`streaming` for a guide.

.. automodule:: peptacular.batch
   :members:

.. automodule:: peptacular.diagnostics
   :members:

Reference data from tacular
---------------------------

Amino acid, element, modification (UNIMOD, PSI-MOD, RESID, XLMOD, GNOme,
UniProt-PTM), ion type, neutral delta, protease and reference molecule lookups
come from `tacular <https://tacular.readthedocs.io/>`_. Import them from
``tacular`` (for example ``from tacular import UNIMOD_LOOKUP``). Only the enums
``pt.IonType``, ``pt.NeutralDelta`` and ``pt.Protease`` are re-exported by
peptacular, because its own signatures take them. They are documented in the
`tacular API reference <https://tacular.readthedocs.io/en/latest/api/index.html>`_.

ProForma JSON
-------------

.. automodule:: peptacular.proforma_json
   :members:

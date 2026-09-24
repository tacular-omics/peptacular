Migrating to 5.0
================

peptacular 5.0 requires tacular 2.0 (``tacular>=2.0,<3``). It removes names that belonged to tacular or were internal, unifies the
digestion API and replaces bare ``ValueError``/``IndexError`` with typed errors. Every
typed error still subclasses ``ValueError``, so ``except ValueError`` keeps working.

Renamed or removed names
------------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - 4.x
     - 5.0
   * - ``pt.AA_LOOKUP``, ``pt.ELEMENT_LOOKUP``, ``pt.UNIMOD_LOOKUP``,
       ``pt.PSIMOD_LOOKUP``, ``pt.PROTEASE_LOOKUP``, ``pt.FRAGMENT_ION_LOOKUP`` and the
       other tacular lookups, ``*Info``/``*Lookup`` classes and literal types
       (``pt.ElementInfo``, ``pt.FragmentIonInfo``, ``pt.IonTypeProperty``, ...)
     - Import from tacular: ``from tacular import AA_LOOKUP``. Only ``pt.IonType``,
       ``pt.NeutralDelta`` and ``pt.Protease`` are still re-exported. The tacular lookups
       themselves changed in tacular 2.0; see its migration guide.
   * - ``pt.Proteases``
     - ``pt.Protease`` (tacular 2.0 renamed the enum)
   * - ``pt.PROTON_MASS``, ``pt.ELECTRON_MASS``, ``pt.NEUTRON_MASS``
     - Same names, now re-exported from ``tacular.constants`` (CODATA 2018). The values
       moved in the 10th decimal place or beyond; see the changelog.
   * - ``pt.parse_composition``
     - ``tacular.parse_composition``
   * - ``pt.parallelMethod``, ``pt.parallelMethodLiteral``
     - ``pt.ParallelMethod``, ``pt.ParallelMethodLiteral``
   * - ``pt.regex_utils``, ``pt.get_regex_match_indices``, ``pt.get_regex_match_range``
     - Removed (private ``peptacular._regex_utils``). Use ``re.finditer``.
   * - ``pt.CV_TO_NAME_PREFIX``, ``pt.CV_TO_ACCESSION_PREFIX``, ``pt.CV_TO_MASS_PREFIX``
     - Removed (private).
   * - ``pt.ReadableProtocol``, ``pt.SupportsStr``, ``pt.handle_number_and_intern_mod``,
       ``pt.Any``, ``pt.SEQUENCE_TYPE``, ``pt.MODIFICATION_*_TYPE``,
       ``pt.GLOBAL_CHARGE_TYPE``, ``pt.ModLocation``, ``pt.MassPropertyMixin``
     - Removed. Annotate with ``str``, ``ProFormaAnnotation`` or ``typing`` types.
   * - ``FLIXIBILITY_SCALES`` (deprecated misspelt alias)
     - ``FLEXIBILITY_SCALES``
   * - ``peptacular.isotope.isotopic_distribution``
     - ``peptacular.isotope.brain_isotopic_distribution`` (``pt.isotopic_distribution``
       is unchanged)
   * - ``pt.isotopic_distribution(annotations=...)``,
       ``pt.estimate_isotopic_distribution(annotations=...)``
     - ``sequence=...``

Digestion
---------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - 4.x
     - 5.0
   * - ``pt.digest(seq, enzyme_regex="trypsin")``, ``pt.cleavage_sites(seq, enzyme_regex=...)``
     - ``pt.digest(seq, enzyme="trypsin")``
   * - ``pt.digest(seq, "([KR])")`` (a regex string)
     - ``pt.digest(seq, re.compile("([KR])"))``. A string is only a protease name; an
       unknown name raises ``UnknownEnzymeError``.
   * - ``EnzymeConfig(enzyme_regex=...)``
     - ``EnzymeConfig(enzyme=...)``. ``EnzymeConfig`` is now frozen.
   * - ``annot.digest(...)``, ``annot.simple_digest(...)``,
       ``annot.sequential_digest(...)``
     - ``annot.digest_spans(...)``, ``annot.simple_digest_spans(...)``,
       ``annot.sequential_digest_spans(...)``. The functional ``pt.digest`` names are
       unchanged.

Behaviour changes
-----------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - 4.x
     - 5.0
   * - ``pt.mass(seqs, 4)`` (positional ``n_workers``, ``chunksize``, ``method``)
     - Keyword-only: ``pt.mass(seqs, n_workers=4)``
   * - ``hash(annotation)``, annotations as dict keys or set members
     - ``ProFormaAnnotation`` is mutable and unhashable. Key on ``annot.serialize()``.
   * - Bare ``ValueError`` from parsers and internals
     - ``ProFormaFormatError`` for bad ProForma, ``PeptacularError`` otherwise
   * - ``IndexError`` for positions outside the sequence
     - ``InvalidPositionError`` (a ``PeptacularError`` and still an ``IndexError``)

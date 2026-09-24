Migrating to 5.0
================

peptacular 5.0 requires tacular 2.0 (``tacular>=2.0,<3``). It removes names that belonged to tacular or were internal, unifies the
digestion API and replaces bare ``ValueError``/``IndexError`` with typed errors. Every
typed error still subclasses ``ValueError``, so ``except ValueError`` keeps working for bad input
values. Arguments of the wrong *type* now raise ``TypeError`` in a few places (listed under
Behaviour changes), which ``except ValueError`` does not catch.

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

Every removed ``pt`` name
~~~~~~~~~~~~~~~~~~~~~~~~~~

82 names left the top-level namespace: 52 that tacular owns and 30 that were internal, renamed or moved to
fastatacular. Search this table for the name your code uses.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - ``pt.<name>`` in 4.x
     - 5.0
   * - ``pt.AA_LOOKUP``
     - ``from tacular import AA_LOOKUP``
   * - ``pt.AALookup``
     - ``from tacular import AALookup``
   * - ``pt.AMINO_ACID_INFOS``
     - ``from tacular import AMINO_ACID_INFOS``
   * - ``pt.AminoAcid``
     - ``from tacular import AminoAcid``
   * - ``pt.AminoAcidInfo``
     - ``from tacular import AminoAcidInfo``
   * - ``pt.Any``
     - ``from typing import Any``
   * - ``pt.CV_TO_ACCESSION_PREFIX``
     - Removed (private).
   * - ``pt.CV_TO_MASS_PREFIX``
     - Removed (private).
   * - ``pt.CV_TO_NAME_PREFIX``
     - Removed (private).
   * - ``pt.Element``
     - ``from tacular import Element``
   * - ``pt.ELEMENT_LOOKUP``
     - ``from tacular import ELEMENT_LOOKUP``
   * - ``pt.ElementInfo``
     - ``from tacular import ElementInfo``
   * - ``pt.ElementLookup``
     - ``from tacular import ElementLookup``
   * - ``pt.fasta``
     - fastatacular; see :ref:`fasta-migration`
   * - ``pt.FASTA_INPUT_TYPE``
     - fastatacular; see :ref:`fasta-migration`
   * - ``pt.FastaFormatError``
     - fastatacular; see :ref:`fasta-migration`
   * - ``pt.FastaSequence``
     - fastatacular; see :ref:`fasta-migration`
   * - ``pt.FRAGMENT_ION_LOOKUP``
     - ``from tacular import FRAGMENT_ION_LOOKUP``
   * - ``pt.FragmentIonInfo``
     - ``from tacular import FragmentIonInfo``
   * - ``pt.FragmentIonLookup``
     - ``from tacular import FragmentIonLookup``
   * - ``pt.get_regex_match_indices``
     - Removed. Use ``re.finditer``.
   * - ``pt.get_regex_match_range``
     - Removed. Use ``re.finditer``.
   * - ``pt.GLOBAL_CHARGE_TYPE``
     - Removed. Annotate with ``str``, ``ProFormaAnnotation`` or ``typing`` types.
   * - ``pt.GNO_LOOKUP``
     - ``from tacular import GNO_LOOKUP``
   * - ``pt.GnoInfo``
     - ``from tacular import GnoInfo``
   * - ``pt.GnoLookup``
     - ``from tacular import GnoLookup``
   * - ``pt.handle_number_and_intern_mod``
     - Removed (internal).
   * - ``pt.IonTypeLiteral``
     - ``from tacular import IonTypeLiteral``
   * - ``pt.IonTypeProperty``
     - ``from tacular import IonTypeProperty``
   * - ``pt.iter_fasta``
     - fastatacular; see :ref:`fasta-migration`
   * - ``pt.MassPropertyMixin``
     - Removed (internal).
   * - ``pt.MODIFICATION_AMBIGUOUS_TYPE``
     - Removed. Annotate with ``str``, ``ProFormaAnnotation`` or ``typing`` types.
   * - ``pt.MODIFICATION_TAG_TYPE``
     - Removed. Annotate with ``str``, ``ProFormaAnnotation`` or ``typing`` types.
   * - ``pt.MODIFICATION_TYPE``
     - Removed. Annotate with ``str``, ``ProFormaAnnotation`` or ``typing`` types.
   * - ``pt.ModLocation``
     - ``from tacular import ModLocation``
   * - ``pt.Monosaccharide``
     - ``from tacular import Monosaccharide``
   * - ``pt.MONOSACCHARIDE_LOOKUP``
     - ``from tacular import MONOSACCHARIDE_LOOKUP``
   * - ``pt.MonosaccharideInfo``
     - ``from tacular import MonosaccharideInfo``
   * - ``pt.MonosaccharideLookup``
     - ``from tacular import MonosaccharideLookup``
   * - ``pt.NEUTRAL_DELTA_DICT``
     - ``from tacular import NEUTRAL_DELTA_DICT``
   * - ``pt.NEUTRAL_DELTA_LOOKUP``
     - ``from tacular import NEUTRAL_DELTA_LOOKUP``
   * - ``pt.NeutralDeltaInfo``
     - ``from tacular import NeutralDeltaInfo``
   * - ``pt.NeutralDeltaLiteral``
     - ``from tacular import NeutralDeltaLiteral``
   * - ``pt.NeutralDeltaLookup``
     - ``from tacular import NeutralDeltaLookup``
   * - ``pt.OboEntity``
     - ``from tacular import OboEntity``
   * - ``pt.OntologyLookup``
     - ``from tacular import OntologyLookup``
   * - ``pt.ORDERED_AMINO_ACIDS``
     - ``from tacular import ORDERED_AMINO_ACIDS``
   * - ``pt.parallelMethod``
     - ``pt.ParallelMethod``
   * - ``pt.parallelMethodLiteral``
     - ``pt.ParallelMethodLiteral``
   * - ``pt.parse_composition``
     - ``from tacular import parse_composition``
   * - ``pt.parse_fasta``
     - fastatacular; see :ref:`fasta-migration`
   * - ``pt.parse_fasta_text``
     - fastatacular; see :ref:`fasta-migration`
   * - ``pt.PROTEASE_LITERALS``
     - ``from tacular import ProteaseLiteral``
   * - ``pt.PROTEASE_LOOKUP``
     - ``from tacular import PROTEASE_LOOKUP``
   * - ``pt.ProteaseInfo``
     - ``from tacular import ProteaseInfo``
   * - ``pt.ProteaseLookup``
     - ``from tacular import ProteaseLookup``
   * - ``pt.Proteases``
     - ``pt.Protease`` (or ``from tacular import Protease``)
   * - ``pt.PROTEASES_DICT``
     - ``from tacular import PROTEASE_DICT``
   * - ``pt.PSIMOD_LOOKUP``
     - ``from tacular import PSIMOD_LOOKUP``
   * - ``pt.PsimodInfo``
     - ``from tacular import PsimodInfo``
   * - ``pt.PsimodLookup``
     - ``from tacular import PsimodLookup``
   * - ``pt.ReadableProtocol``
     - Removed (internal).
   * - ``pt.REFMOL_LOOKUP``
     - ``from tacular import REFMOL_LOOKUP``
   * - ``pt.RefMolID``
     - ``from tacular import RefMolID``
   * - ``pt.RefMolInfo``
     - ``from tacular import RefMolInfo``
   * - ``pt.RefMolLiteral``
     - ``from tacular import RefMolLiteral``
   * - ``pt.RefMolLookup``
     - ``from tacular import RefMolLookup``
   * - ``pt.regex_utils``
     - Removed (private ``peptacular._regex_utils``).
   * - ``pt.RESID_LOOKUP``
     - ``from tacular import RESID_LOOKUP``
   * - ``pt.ResidInfo``
     - ``from tacular import ResidInfo``
   * - ``pt.ResidLookup``
     - ``from tacular import ResidLookup``
   * - ``pt.SEQUENCE_TYPE``
     - Removed. Annotate with ``str``, ``ProFormaAnnotation`` or ``typing`` types.
   * - ``pt.SupportsStr``
     - Removed (internal).
   * - ``pt.UNIMOD_LOOKUP``
     - ``from tacular import UNIMOD_LOOKUP``
   * - ``pt.UnimodInfo``
     - ``from tacular import UnimodInfo``
   * - ``pt.UnimodLookup``
     - ``from tacular import UnimodLookup``
   * - ``pt.UNIPROT_PTM_LOOKUP``
     - ``from tacular import UNIPROT_PTM_LOOKUP``
   * - ``pt.UniprotPtmInfo``
     - ``from tacular import UniprotPtmInfo``
   * - ``pt.UniprotPtmLookup``
     - ``from tacular import UniprotPtmLookup``
   * - ``pt.XLMOD_LOOKUP``
     - ``from tacular import XLMOD_LOOKUP``
   * - ``pt.XlModInfo``
     - ``from tacular import XlmodInfo``
   * - ``pt.XlModLookup``
     - ``from tacular import XlmodLookup``

Names added to ``pt``: ``HasSequence``, ``ParallelMethod``, ``ParallelMethodLiteral``, ``PeptacularKeyError``,
``Protease``, ``UnknownElementError`` and ``UnknownEnzymeError``.

Renamed misspelt names
----------------------

In ``peptacular.property.data`` (and the matching ``property`` enums):

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - 4.x
     - 5.0
   * - ``AMIGUOUS_AMINO_ACID_MAP``
     - ``AMBIGUOUS_AMINO_ACID_MAP``
   * - ``surface_accessiblility_janin``
     - ``surface_accessibility_janin``
   * - ``hphob_agros`` (enum member ``AGROS``)
     - ``hphob_argos`` (``ARGOS``)
   * - ``hphob_adoberin`` (enum member ``ADOBERIN``)
     - ``hphob_aboderin`` (``ABODERIN``)
   * - ``FLIXIBILITY_SCALES`` (deprecated alias)
     - ``FLEXIBILITY_SCALES``

.. _fasta-migration:

FASTA reading moved to fastatacular
-----------------------------------

peptacular no longer reads files. Install ``fastatacular`` and pass its entries straight to
peptacular: any object with a ``.sequence`` string (the ``HasSequence`` protocol) is accepted
by the sequence functions, ``batch``, ``iter_batch`` and ``diagnose``.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - 4.x
     - 5.0
   * - ``pt.parse_fasta(path)``
     - ``fastatacular.read_fasta(path)``
   * - ``pt.iter_fasta(path)``
     - ``with fastatacular.FastaReader(path) as reader: for entry in reader: ...``
   * - ``pt.parse_fasta_text(text)``
     - ``fastatacular.read_fasta(io.StringIO(text))``
   * - ``pt.FastaSequence`` (``.header``, ``.sequence``)
     - ``fastatacular.SequenceEntry`` (``.raw_header`` or ``.identifier``, ``.sequence``)
   * - ``pt.FastaFormatError``
     - ``fastatacular.FastaParseError``
   * - ``pt.FASTA_INPUT_TYPE``, ``pt.fasta``
     - Removed.

Two behaviour differences: peptacular uppercased sequences and silently dropped empty entries;
fastatacular keeps the case as written and raises on an empty entry.

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
   * - ``EnzymeConfig("trypsin", 1, True)``, ``config.semi_enzymatic``
     - ``EnzymeConfig("trypsin", missed_cleavages=1, semi=True)``, ``config.semi``
       (the same name as ``digest(..., semi=)``). Options are keyword-only.
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
   * - ``except ValueError`` around ``get_mod_type`` or a ``mods=`` argument of the wrong type
     - ``TypeError`` now (a string naming no mod type still raises ``PeptacularError``)
   * - ``pt.digest(seq, enzyme="")`` for a nonspecific digest
     - ``UnknownEnzymeError``. Use ``enzyme="unspecific"`` or ``pt.nonspecific_digest``.
   * - ``calculate_composition=True``
     - ``calculate_with_composition=True``
   * - ``brain_isotopic_distribution(chemical_formula, ..., charge_state=2)``
     - ``brain_isotopic_distribution(formula, *, ..., charge=2)``
   * - ``annot.coverage(annotations=...)``, ``percent_coverage``, ``modification_coverage``
     - ``subsequences=...``
   * - ``annot[3]``
     - ``UnsupportedOperationError``. Use ``annot[3:4]`` or ``annot.stripped_sequence[3]``.
   * - ``peptacular.utils.get_mods``
     - Removed. Use ``annot.get_mods()`` or ``pt.get_mods``.
   * - ``hash(batch_result)``
     - ``BatchResult`` is compared by value and is not hashable when it holds an annotation, list or dict.
   * - ``fragment.mz = ...`` (mutating a ``Fragment``)
     - ``Fragment`` is immutable and raises ``dataclasses.FrozenInstanceError``. Use ``fragment.replace(mass=...)`` (constructor names; ``mz`` is derived from ``mass``).
   * - ``fragment.losses``, ``fragment.asdict()["losses"]``, MCP fragment ``losses``
     - ``fragment.deltas``, ``asdict()["deltas"]``, MCP ``deltas``.
   * - ``Fragment(ion_type, position, mass, mono, charge, adducts, ...)``
     - Options after ``charge_state`` are keyword-only.
   * - Fragments compared by identity
     - Compared and hashed by value.
   * - ``fast_fragment`` m/z values
     - Move by -1.4e-8 Da per charge (charge carrier is now H - e, as in ``fragment()``). ``fragment()`` and ``mass()`` are unchanged.
   * - ``ProFormaAnnotation("PEPTIDE", None, None, ...)``, ``Interval(1, 3, True, mods)``
     - Every option after ``sequence`` (and after ``start``, ``end``) is keyword-only:
       ``Interval(1, 3, ambiguous=True, mods=mods)``.
   * - ``interval.set_mods(m, True)``, ``append_mod(m, True, False)``, ``extend_mods(m, True)``
     - ``validate=`` and ``inplace=`` are keyword-only. ``append_mod`` returns the interval
       (the copy when ``inplace=False``; 4.x returned None).
   * - ``frag.to_mzpaf(False)``, ``frag.serialize("mzpaf")``
     - ``frag.to_mzpaf(include_sequence=False)``, ``frag.serialize(format="mzpaf")``
   * - ``annot.prop.calc_property(scale, "avg", ...)`` and ``property_windows`` /
       ``property_partitions`` options
     - Keyword-only after ``scale``: ``calc_property(scale, normalize=True)``
   * - ``mod.get_mass(False)``, ``formula.get_mass(False)`` (any component)
     - ``get_mass(monoisotopic=False)``, as in tacular 2.0. Also keyword-only:
       ``ChargedFormula.from_string``/``serialize``/``from_composition`` and
       ``FormulaElement.from_string`` options, ``ModificationTags.validate(all_tags=)``.
   * - ``PeptidoformIon.get_mass()`` / ``get_composition()`` (always raised
       ``NotImplementedError``)
     - ``get_mass`` raises ``UnsupportedOperationError`` with a hint; ``get_composition``
       is removed. Use ``pt.parse`` / ``pt.parse_chimeric`` and the annotation methods.
   * - mzPAF neutral-loss labels ``-H3CON``, ``-H2CO2``
     - Canonical names: ``-HCONH2``, ``-HCOOH``. Other formulas are written in Hill order
       (``+NaS``, not ``+SNa``). Every other 4.2.0 label, including ``-NH3``, is unchanged.
   * - ``fragment(..., ion_types="by")`` (a string of letters meant b and y)
     - A string is one ion type. Write ``ion_types=("b", "y")``.
   * - ``fragment(..., neutral_deltas=["H3PO4"])`` on a peptide with an unmodified S/T/Y
       raised ``InvalidAdjustmentError``
     - The impossible loss is skipped for that ion; possible losses are still produced.
       A delta passed explicitly in ``deltas=`` still raises.
   * - ``pt.parse(b"PEPTIDE")``, ``pt.parse(None)`` (``TypeError: ... has no len()``)
     - ``TypeError`` naming the accepted inputs. An object with a str ``sequence``
       (a FASTA entry) is now accepted.

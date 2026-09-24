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

MCP server
----------

The ``peptacular-mcp`` tools use the library's argument and record names. The
response ``contract_version`` is ``2.0``. The old names are rejected, since
unknown arguments are errors.

.. list-table::
   :header-rows: 1

   * - 4.x MCP name
     - 5.0 MCP name
   * - ``fragment_peptides`` ``ion_series``
     - ``ion_types``
   * - ``fragment_peptides`` ``isotope_offsets``
     - ``isotopes``
   * - ``fragment_peptides`` ``include: ["label"]``, row ``label``
     - ``include: ["mzpaf"]``, row ``mzpaf``
   * - ``enumerate_modifications`` ``max_variable_modifications``
     - ``max_variable_mods``
   * - fragment rows ``ion_series``, ``ordinal``, ``charge``
     - ``ion_type``, ``position``, ``charge_state``
   * - fragment rows ``ion_mass_da``, ``neutral_mass_da``
     - ``mass``, ``neutral_mass``
   * - ``get_reference(topic="ions")`` rows ``ion_series``
     - ``ion_type``
   * - fragment formula delta ``"-H3PO4"`` (failed the call)
     - A loss, as ``"H3PO4"``. ``"+HPO3"`` is a gain. Ions that cannot lose it are skipped.

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
   * - Monoisotopic ``mass()``, ``mz()``, ``fragment()`` values of charged ions
     - A proton charge carrier is now CODATA ``PROTON_MASS`` (was H - e), as mzPAF 4.4.1, pyteomics and OpenMS use. Mass moves by +1.43e-8 Da per charge and m/z by +1.43e-8 Da (-1.43e-8 Da for deprotonated ions). ``fast_fragment`` already used ``PROTON_MASS`` and is unchanged; it now matches ``fragment()`` to 1e-9 Da. Average masses are unchanged.
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
       (``+NaS``, not ``+SNa``). ``-NH3`` is unchanged. The other 5.0 label changes are
       in the rows below.
   * - ``fragment(..., ion_types="by")`` (a string of letters meant b and y)
     - A string is one ion type. Write ``ion_types=("b", "y")``.
   * - ``fragment(..., neutral_deltas=["H3PO4"])`` on a peptide with an unmodified S/T/Y
       raised ``InvalidAdjustmentError``
     - The impossible loss is skipped for that ion; possible losses are still produced.
       A delta passed explicitly in ``deltas=`` still raises.
   * - ``pt.parse(b"PEPTIDE")``, ``pt.parse(None)`` (``TypeError: ... has no len()``)
     - ``TypeError`` naming the accepted inputs. An object with a str ``sequence``
       (a FASTA entry) is now accepted.
   * - ``fragment.to_mzpaf()`` with a gain (``deltas={"H2O": -1}``) or a numeric delta
     - A gain is written ``+H2O`` (4.x wrote ``-H2O``). A numeric delta is written as a signed mass (``b2-34.0``) instead of raising.
   * - ``fragment.to_mzpaf()`` at a negative charge: ``y3{IDE}^1`` for z=-1, ``^2`` for z=-2
     - Signed: ``y3{IDE}^-1``, ``^-2``, so the label reads back to the same m/z.
       ``to_mzpaf(signed_charge=False)`` restores the magnitude-only form. That form is
       mzPAF 1.0.1 section 4.8 and is only valid next to negative-mode spectrum metadata.
       At z=-1 it has no charge suffix, so on its own it reads as +1.
   * - ``fragment.to_mzpaf()`` for an immonium ion wrote only the residue's own mod
       (``IP`` for ``[Acetyl]-PEP``, ``<[Oxidation]@P>PEP`` or ``<13C>PEP``)
     - A terminal mod or a global fixed mod on the residue is written as the immonium
       mod: ``IP[Acetyl]``, ``IP[Oxidation]``. A global isotope label is written as isotope
       shifts, one per labelled atom of the final ion: ``IP+4i13C``, ``IP+6i2H^-1`` for
       ``<D>P`` at charge -1, ``IK-NH3+i15N``. Two or more mods raise ``PeptacularError``. The mod tag is the
       plain name: ``P[U:Oxidation]`` gives ``IP[Oxidation]`` (4.x wrote ``IP[U:Oxidation]``).
   * - ``fragment.to_mzpaf()`` for ax/bx internal ions: ``-2H``, ``+CO-2H``
     - Hill order, like every other delta: ``-H2``, ``+CO-H2``.
   * - ``charge="H:z-1"`` (hydride): ``is_protonated`` True, mzPAF ``y3{IDE}^-1``
     - A hydride is an adduct, not a proton: ``is_protonated`` is False, and the mzPAF is
       ``y3{IDE}[M+H]^-1``. The mass is unchanged.
   * - d/da/db/w/wa/wb ion of a modified residue (``PEPV[Oxidation]K`` d4 gave ``d4{PEPV[Oxidation]}``)
     - Not defined: ``frag()`` raises ``PeptacularError`` (explicit or global fixed mod),
       ``fragment()`` leaves the ion out. v ions still drop the mod.
   * - ``fragment.to_mzpaf()`` of an uncharged fragment: ``b3{PEP}`` (reads as +1)
     - Raises ``PeptacularError``; build the ion with ``charge=1``.
   * - Full-length d/da/db ions ignored the C-terminal mod, v/w/wa/wb ions the N-terminal mod
     - Every full-length ion type carries both terminal mods, as a/b/c/x/y/z already did.
       Only those ions' masses change.

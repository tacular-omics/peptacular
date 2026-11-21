Features
========

ProForma 2.1 Parsing
~~~~~~~~~~~~~~~~~~~~

Parse, serialize, and modify ProForma sequences.

- Full support for Base-ProForma and Level 2-ProForma and select other features
- Named modifications from Unimod / PSI-MOD / RESID / XLMOD / GNOme databases
- Delta mass and formula-based modifications
- Terminal, labile, and ambiguous modifications
- Charge states and adducts
- Global isotope labels and fixed modifications

Mass Calculations
~~~~~~~~~~~~~~~~~

Calculate monoisotopic and average masses with comprehensive modification support.

- Monoisotopic and average mass calculations
- Support for isotope labeling and neutral losses
- Charge state and adduct handling
- Precursor and fragment ion masses
- Elemental composition tracking

Isotopic Distributions
~~~~~~~~~~~~~~~~~~~~~~

Generate aggregated theoretical isotopic patterns with exact center masses and abundance filtering.

- High-performance isotopic pattern generation
- Adaptive envelope sizing and configurable abundance thresholds
- Support for modifications (must have composition data)

Protein Digestion
~~~~~~~~~~~~~~~~~

Simulate enzymatic digestion with various proteases and cleavage rules.

- Built-in protease database (trypsin, chymotrypsin, pepsin, etc.)
- Missed cleavage support
- Semi-specific and non-specific digestion
- Custom protease definitions

Fragment Generation
~~~~~~~~~~~~~~~~~~~

Create theoretical fragment ions for MS/MS analysis.

- b, y, a, x, c, z ion types
- Internal fragments (ax, by, etc.) & Immonium ions
- Multiple charge states / adducts
- Neutral losses (H2O, NH3, custom)
- Isotopes

**Fast Fragment**

For high-throughput workflows, ``fast_fragment()`` uses a prefix/suffix-sum
algorithm to compute fragment m/z values directly, without constructing
:class:`~peptacular.Fragment` objects.  It is faster than ``fragment()`` but
does not support neutral losses, isotope shifts, adduct charges, or custom
deltas.  Returns a ``dict`` mapping ``(IonType, charge)`` to a list of m/z
values ordered from fragment position 1 to N.

.. code-block:: python

   import peptacular as pt

   # OOP method
   mz_map = pt.parse("PEPTIDE").fast_fragment(ion_types=["b", "y"], charges=[1, 2])
   b1_mzs = mz_map[(pt.IonType.B, 1)]  # list of 7 floats

   # Functional API — supports batch list inputs (auto-parallelised)
   results = pt.fast_fragment(["PEPTIDE", "ACDEFGHIK"], ion_types=["y"], charges=[1])

**Fragment columns (numpy)**

``fragment_arrays()`` returns every ion ``fragment()`` would, with the same values, as a
``dict`` of numpy arrays (one row per ion, with ``peptide_index``, ``ion_type``,
``position``, ``charge_state``, ``mz``, ``mass`` and delta/isotope columns) that
``polars.DataFrame`` or ``pyarrow.table`` take as is. Unlike ``fast_fragment()`` it keeps
every ``fragment()`` option. Needs ``pip install "peptacular[numpy]"``; see
:doc:`records`.

Property Calculations
~~~~~~~~~~~~~~~~~~~~~

Calculate various physicochemical properties of peptides.

- Hydrophobicity (multiple scales: Kyte-Doolittle, Eisenberg, etc.)
- Isoelectric point (pI)
- Aromaticity
- Aliphatic index
- Instability index
- Molar extinction coefficient
- GRAVY score
- Secondary structure predictions
- Custom property scales

Format Conversion
~~~~~~~~~~~~~~~~~

Import and export sequences from popular proteomics tools.

**Supported Formats:**

- IP2 (flanking amino acids notation)
- DIANN (modification notation)
- MaxQuant (modified sequence notation)
- Skyline (modified sequence notation)
- Casanovo (de novo sequencing output)
- MS2PIP (fragment prediction input/output)
- MSFragger (peptide and assigned modifications pairs)

Parallel Processing
~~~~~~~~~~~~~~~~~~~

Automatic multiprocessing for batch operations.

- Automatic parallelization for list inputs
- No code changes needed
- Applies to all major functions (mass, digest, fragment, etc.)
- Efficient memory usage, cached modifications, and frozen dataclasses

Design Principles
-----------------

Lazy Loading / Caching
~~~~~~~~~~~~~~~~~~~~~~~

Modifications are stored as strings and only parsed when needed for calculations. 
This minimizes overhead and improves performance for operations that don't require 
mass and composition data.

.. code-block:: python

   import peptacular as pt
   
   # Modifications are not parsed until needed
   peptide = pt.parse("PEM[Oxidation]TIDE")
   
   # First call parses the modifications and caches the return values of 'Oxidation'
   mass1 = peptide.mass()
   
   # Second call uses cached values directly (No re-parsing)
   mass2 = peptide.mass()


**Validation**

By default, parsing does not validate inputs. This allows users to potentially create annotations
with invalid modifications/sequences. If you need strict validation:

.. code-block:: python

   import peptacular as pt
   
   # Strict validation
   a = pt.parse('PEPTIDE', validate=True)
   a.static_mods = 'INVALID_MOD'  # raises ProFormaFormatError (a PeptacularError)

   a.validate = False  # Disable validation
   a.static_mods = 'INVALID_MOD'  # Works without error

   # get validation status
   print(a.validate)

``pt.validate()`` is the whole-string check, and it is stricter than ``parse(validate=True)``:
on top of every field it reads the rules that span the whole annotation, among them the
``#label`` rules of ProForma 2.1 sections 7.6.2, 7.9 and 9.2. A label must be defined by
attaching the modification to one position, one label may name one modification only, no two
labels may differ only in case, and a labile modification's label may not be reused for another
modification. Referencing a labile label in the sequence, which section 7.9 is for, stays valid,
and a label written as its own element of a pipe-separated bracket, ``T[Phospho|#g1]``, counts as
that modification's label.

.. testcode::

   import peptacular as pt

   print(pt.validate("EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK"))  # a defined group
   print(pt.validate("PEPT[#g1]IDES[#g1]"))                            # #g1 never defined
   print(pt.validate("PEPT[Phospho#g1]S[#g1]IDES[Formyl#G1]K[#G1]"))   # #g1 against #G1
   print(pt.validate("{HexNAc#g1}PEPN[#g1]ITE"))                       # section 7.9
   print(pt.validate("{Glycan:Hex#g1}PEPT[Phospho#g1]IDES[#g1]"))      # the labile label reused

.. testoutput::

   True
   False
   False
   True
   False

A ``True`` is not a promise that a mass can be calculated. ``B``, ``Z``, ``J`` and a bare ``X``
are valid ProForma residues and an empty sequence is valid, while ``pt.mass()`` raises for a bare
``X``, for ``B`` or ``Z``, and for an empty sequence.
   

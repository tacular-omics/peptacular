Quick Start
===========

Installation
------------

Install Peptacular using pip:

.. code-block:: bash

   pip install peptacular

Basic Usage
-----------

Import as a namespace:

.. testcode::

   import peptacular as pt

.. testsetup:: *

   import peptacular as pt

Parse/Serialize
~~~~~~~~~~~~~~~

`pt.parse` will return a ProFormaAnnotation object representing the sequence with modifications.

.. testcode::
   
   # Parse a ProForma sequence
   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   print(peptide.charge_state)
   print(peptide.serialize())

.. testoutput::

   2
   PEM[Oxidation]TIDE/2


Parse and Modify Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A number of methods exist to modify the ProFormaAnnotation object after parsing. You can set, append, extend, and remove modifications.

There are a lot of accepted input values for these methods, see the API documentation for details. 
But generally, strings representing modification names or identifiers, and float / int values are accepted.

append / extend methods will add to existing modifications (increasing the count if the modification already exists at that position),
while set methods will overwrite existing modifications at that position.

Annotations have the following modification types: static, isotope, labile, nterm, cterm, internal, interval, charge (charge state or adducts).

.. testcode::

   # Parse ProForma sequences
   simple = pt.parse("PEPTIDE")
   modified = pt.parse("[Acetyl]-PEP[Oxidation]TIDE/2")
   labeled = pt.parse("<13C>PEPTIDE")

   # Programmatic modification
   annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
   annot.set_nterm_mods({"Acetyl": 1})
   annot.set_internal_mods_at_index(2, {"Oxidation": 1})
   annot.set_charge(2)
   print(annot.serialize())

.. testoutput::

   [Acetyl]-PEP[Oxidation]TIDE/2


Fragment Generation
~~~~~~~~~~~~~~~~~~~

The fragment method requires that you provide a list of ion types and charge states to generate. It will return a list of FragmentAnnotation objects.
If no charge states are provided, it will default to the charge state of the peptide, or charge state 1 if the peptide does not have a charge state assigned.

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Fragment the peptide
   for frag in peptide.fragment(ion_types=['b', 'y'], charges=[1, 2])[:5]:
       print(f"{frag.ion_type}{frag.position}+{frag.charge_state}: {frag.mz:.3f}")

.. testoutput::

   b1+1: 98.060
   b2+1: 227.103
   b3+1: 374.138
   b4+1: 475.186
   b5+1: 588.270

Protein Digestion
~~~~~~~~~~~~~~~~~

The digest method will return spans representing the (start, end, missed_cleavages) of each peptide generated from the digestion.
You can readily create peptide annotations by slicing the protein annotation with the spans.

A number of proteases are available in pt.Proteases, but you can also provide your own cleavage rules via a regex pattern.

.. testcode::

   # Digest a protein
   protein = pt.parse("PEM[Oxidation]TRPEPTIDEKPEPTIDEIDE/2")
   for span in protein.digest(pt.Proteases.TRYPSIN, missed_cleavages=1):
       print(f"  {protein[span].serialize()}")

.. testoutput::
   :hide:

   ...

There is an additional digest option which is based on amino acids rather than regex patterns.

.. testcode::

   # Digest a protein
   protein = pt.parse("PEM[Oxidation]TRPEPTIDEKPEPTIDEIDE/2")
   for span in protein.simple_digest(cleave_on="KR", restrict_after="P", missed_cleavages=1):
       print(f"  {protein[span].serialize()}")

.. testoutput::
   :hide:

   ...


Mass / MZ / Composition
~~~~~~~~~~~~~~~~~~~~~~~

There are a number of methods available to calculate masses, m/z values, and elemental compositions of peptides. 

This is where modifications are parsed and applied to the calculations. So previously parsed ProFormaAnnotation objects just required valid syntax, 
but for mass/comp calculations the modifications must also be recognized. For mass calculations to work, all modifications must be mass resolvable. For compositions,
all modifications must be composition resolvable (This means delta mass modifications eg. [+15.99] will not work for compositions).

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Calculate masses with different parameters
   precursor_mass = peptide.mass(ion_type="p")
   neutral_mass = peptide.neutral_mass(ion_type="p")
   mz_2plus = peptide.mz(charge=2)

   # With isotopes and losses
   mass_c13 = peptide.mass(isotopes=1)
   mass_h2o_loss = peptide.mass(deltas={'H-2O-1': 1})

   # Adduct masses
   mass_na = peptide.mass(charge='Na:z+1')

Isotopic Patterns
~~~~~~~~~~~~~~~~~

Peptacular generates aggregated nominal isotope envelopes for peptides and fragments with the BRAIN algorithm. Each result contains a nominal neutron offset, relative abundance, and exact probability-weighted center mass. This relies on parsing a valid composition for the peptide and its modifications.
If the composition cannot be resolved, use estimate_isotopic_distribution to predict the envelope from peptide averagine and the molecular mass.

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Generate isotopic distribution
   dist = peptide.isotopic_distribution(
       max_isotopes=5,
       min_abundance_threshold=0.01
   )

   # Omit max_isotopes to size the envelope adaptively
   adaptive_dist = peptide.isotopic_distribution(
       min_abundance_threshold=0.001
   )

   # For fragments with modifications
   dist_fragment = peptide.isotopic_distribution(
       ion_type='y',
       charge=2,
       isotopes=1,
       deltas={'H-2O-1': 1}
   )


Physicochemical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~
These are based solely on the amino acid sequence, modifications are not considered in these calculations.

Properties are accessible via the `prop` property. See the API documentation for a full list of available properties and methods. There are also options to calculate custom properties via the calc_property method.


.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Simple properties
   hydro = peptide.prop.hydrophobicity
   pi = peptide.prop.pi
   arom = peptide.prop.aromaticity

   # Secondary structure prediction
   ss = peptide.prop.secondary_structure()
   print(f"Alpha helix: {ss['alpha_helix']:.1f}%")
   print(f"Beta sheet: {ss['beta_sheet']:.1f}%")

.. testoutput::
   :hide:

   ...

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Custom property calculations
   custom = peptide.prop.calc_property(
       scale=pt.HydrophobicityScale.KYTE_DOOLITTLE,
       aggregation_method='avg'
   )

Format Conversion
~~~~~~~~~~~~~~~~~

Peptacular supports conversion from and to a number of popular peptide sequence formats.

.. testcode::

   # From other formats
   ip2_seq = pt.convert_ip2_sequence('K.PEPTIDE.K')
   diann_seq = pt.convert_diann_sequence('_PEM[Carbamidomethyl]TIDE_')
   casanovo_seq = pt.convert_casanovo_sequence('+43.006PEPTIDE')

   # To MS2PIP format
   peptide = pt.parse("PEM[Oxidation]TIDE")
   unmod_seq, mod_str = peptide.to_ms2_pip()

Batch Processing
~~~~~~~~~~~~~~~~

Most functions in Peptacular support batch processing of multiple sequences via lists. Small automatic batches run sequentially,
while larger batches use parallel execution. The sequence API functions support strings and annotation objects as input.

.. testcode::

   # Process multiple sequences in parallel
   sequences = ['PEPTIDE', 'PEPTIDES', 'PEPTIDEST']

   # Automatic parallelization for lists
   masses = pt.mass(sequences)
   compositions = pt.comp(sequences)
   annotations = pt.parse(sequences)

   # Properties
   hydrophobicity = pt.hydrophobicity(sequences)
   pi_values = pt.pi(sequences)

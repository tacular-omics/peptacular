=====================================
Mass Calculation
=====================================

Peptacular calculates the molecular mass and isotopic patterns for amino acid
(AA) sequences. This page presents the mathematical framework that underlies
these calculations.

Base Mass
---------

The base mass :math:`M_{base}` of an AA sequence containing modifications is
calculated as the sum of all constituent molecular components:

.. math::

   M_{base} = \sum_{i=1}^{n} m_{AA_i} + M_{N} + M_{C} + M_{S} + M_{I} + M_{R} + M_{U} + \mathbb{1}_{\textrm{precursor} \, \lor \, \textrm{neutral}} \cdot M_{L}

where:

- :math:`n` is the sequence length
- :math:`m_{AA_i}` is the mass of amino acid at position :math:`i`
- :math:`M_{N}` is the total mass of N-terminal modifications
- :math:`M_{C}` is the total mass of C-terminal modifications
- :math:`M_{S}` is the total mass of static/fixed modifications (applied to sequence)
- :math:`M_{I}` is the total mass of position-specific modifications
- :math:`M_{R}` is the total mass of modifications within defined sequence intervals
- :math:`M_{U}` is the total mass of modifications with unknown positions
- :math:`M_{L}` is the total mass of labile modifications
- :math:`\mathbb{1}_{\textrm{precursor} \, \lor \, \textrm{neutral}}` is an indicator function: 1 for precursor and neutral ion types, 0 otherwise

Labile modifications are included only for precursor and neutral ion types, since they are
lost during fragmentation.

**Modification masses.** A named modification (a Unimod, PSI-MOD, RESID, XLMOD or GNO
name or accession, e.g. ``[Oxidation]`` or ``[UNIMOD:35]``) contributes its listed
database mass (tacular's ``monoisotopic_mass`` / ``average_mass``, e.g. 15.994915 for
Oxidation). The elemental composition is used only when there is no listed mass
(``Formula:`` mods, glycans, unlisted entries). Because listed masses are rounded (Unimod
gives 6 decimals), the mass of a named mod's ``composition`` can differ from ``mass()`` by
up to about 1e-6 Da. Loss, isotope and deprotonated ions use the same listed masses, so a
loss or isotope peak is exactly the plain ion plus its delta; only
``calculate_with_composition=True`` and global isotope labels (``<13C>``) sum compositions.
ProForma text is never rewritten: names stay names and numbers stay as written.

Neutral Mass
------------

The neutral mass :math:`M_{\textrm{neutral}}` is calculated by combining the base
mass with ion-type adjustments, isotope modifications, and neutral deltas:

.. math::

   M_{\textrm{neutral}} = M_{\textrm{base}} + M_{\textrm{ion}} + M_{\textrm{isotope}} + M_{\textrm{ndelta}}

where:

- :math:`M_{\textrm{base}}` is the peptide base mass from the previous section
- :math:`M_{\textrm{ion}}` is the ion-type-specific mass offset
- :math:`M_{\textrm{isotope}}` is the mass shift from a specific isotopic species
- :math:`M_{\textrm{ndelta}}` is the mass change from neutral losses/gains

Mass-to-charge Ratio
--------------------

The mass-to-charge (*m/z*) ratio is calculated by incorporating charge carriers
and electron mass corrections to the neutral mass:

.. math::

   \frac{m}{z} = \frac{M_{\textrm{neutral}} + M_{\textrm{adduct}} - z \cdot m_e}{|z|}

where:

- :math:`M_{\textrm{neutral}}` is the neutral fragment mass
- :math:`M_{\textrm{adduct}}` is the total mass of charge carriers
- :math:`z` is the total charge state (adduct charge plus any charge contributed by modifications), which may be negative for negative-mode ions
- :math:`m_e = 0.000548579909065` Da (electron mass, CODATA 2018)

The denominator uses :math:`|z|` so that negative charge states still produce a
positive *m/z* value.

**The proton charge carrier.** A default (protonated) monoisotopic charge adds the
CODATA ``PROTON_MASS`` (1.007276467 Da), as mzPAF (section 4.4.1), pyteomics and OpenMS
do. An ion's elemental composition counts one hydrogen atom per charge, and
``HYDROGEN_MASS - ELECTRON_MASS`` (1.007276452 Da) is lighter than a proton by the
hydrogen 1s binding energy, ``HYDROGEN_BINDING_MASS`` (1.43e-8 Da). A mass summed from a
charged composition therefore adds ``HYDROGEN_BINDING_MASS`` once per net proton (a
deprotonation subtracts it), so the composition path and the mass path give the same
number. ``mass()``, ``mz()``, ``fragment()`` and ``fast_fragment()`` all agree. Other
adducts (``Na:z+1``, a hydride ``H:z-1``) are their atoms minus or plus electrons, with
no binding term. Average masses use the average hydrogen mass minus one electron.

Negative charge
---------------

**A negative charge is deprotonation.** ``PEPTIDE/-2`` is the ion
:math:`[M-2H]^{2-}`: two protons are removed, not added.

This departs from the notation's own wording, deliberately. ProForma 2.1 section 11.5
reads:

.. code-block:: text

   The charge of a peptidoform ion can be specified after the definition. This specifies
   the charge of the peptidoform ion, excluding any local charges of elemental formulas
   (section 11.1). Any positive charges are defined as addition of protons H+, any
   negative charges are defined as addition of electrons e-.

peptacular follows that sentence for positive charges, which are protons, and parts with
it for negative ones, because the ion that negative-mode electrospray produces is the
deprotonated peptide, not the peptide with extra electrons on it. Reading a negative
charge as added electrons would make ``pt.mass("PEPTIDE/-2")`` heavier than the neutral
peptide, which is not an ion any negative-mode spectrum contains.

**The size of the difference is one neutral hydrogen atom per charge**, which is
``PROTON_MASS + ELECTRON_MASS``, about 1.007825 Da, so 2.01565 Da at charge 2:

.. testcode::

   import peptacular as pt

   print(round(pt.mass("PEPTIDE/-2"), 6))                              # what peptacular returns
   print(round(pt.parse("PEPTIDE").neutral_mass(), 6))                 # the neutral peptide
   print(round(pt.parse("PEPTIDE").neutral_mass() + 2 * pt.ELECTRON_MASS, 6))  # the specification's reading

.. testoutput::

   797.345411
   799.359964
   799.361061

**Getting the specification's reading.** There is no charge-carrier spelling for it. A
charge carrier is a charged formula of elements (section 11.1) and an electron is not an
element, so ``PEPTIDE/[e:z-1]`` and ``PEPTIDE/[e-1:z-1]`` both raise
``ProFormaFormatError``, and so does an empty or zero-count formula. A carrier can only
state the deprotonation peptacular already performs, and then gives the same number as
``/-2``:

.. testcode::

   print(round(pt.mass("PEPTIDE/[H-1:z-1^2]"), 6))

   try:
       pt.mass("PEPTIDE/[e:z-1]")
   except pt.ProFormaFormatError as exc:
       print(type(exc).__name__)

.. testoutput::

   797.345411
   ProFormaFormatError

If you need :math:`[M+ne]^{n-}`, add ``n * pt.ELECTRON_MASS`` to ``neutral_mass()``
yourself, as the third line above does. Nothing in peptacular computes it.

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

   M_{base} = \sum_{i=1}^{n} m_{AA_i} + M_{N} + M_{C} + M_{S} + M_{I} + M_{R} + M_{U} + \mathbb{1}_{\textrm{precursor}} \cdot M_{L}

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
- :math:`\mathbb{1}_{\textrm{precursor}}` is an indicator function: 1 for precursor ions, 0 for fragment ions

Labile modifications are included only for precursor and neutral ion types.

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

   \frac{m}{z} = \frac{M_{\textrm{neutral}} + M_{\textrm{adduct}} - z \cdot m_e}{z}

where:

- :math:`M_{\textrm{neutral}}` is the neutral fragment mass
- :math:`M_{\textrm{adduct}}` is the total mass of charge carriers
- :math:`z` is the total charge state
- :math:`m_e = 0.0005485799` Da (electron mass)

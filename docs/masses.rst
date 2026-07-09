=====================================
Ion Types
=====================================

Peptacular supports the fragment ion types below, passed as the ``ion_types`` argument
to methods like ``fragment()`` and ``mass()``. Each ion type is defined by a delta
formula that is added to the neutral peptide backbone mass — see :doc:`mass_calculation`
for how these deltas fit into the overall mass equation.

Standard Fragment Ions
-----------------------

These are the most commonly used fragment ions in peptide mass spectrometry:

.. list-table::
   :header-rows: 1
   :widths: 15 15 25 45

   * - Ion Type
     - String Value
     - Delta Formula
     - Description
   * - **Precursor**
     - ``"p"``
     - H₂O
     - Intact peptide with water addition
   * - **Neutral**
     - ``"n"``
     - (none)
     - Neutral peptide backbone (no modifications)
   * - **a-ion**
     - ``"a"``
     - C⁻¹O⁻¹
     - N-terminal fragment
   * - **b-ion**
     - ``"b"``
     - (none)
     - N-terminal fragment
   * - **c-ion**
     - ``"c"``
     - H₃N
     - N-terminal fragment
   * - **x-ion**
     - ``"x"``
     - CO₂
     - C-terminal fragment
   * - **y-ion**
     - ``"y"``
     - H₂O
     - C-terminal fragment
   * - **z-ion**
     - ``"z"``
     - H⁻¹N⁻¹O
     - C-terminal fragment
   * - **z• radical**
     - ``"z."``
     - N⁻¹O
     - Radical z-ion
   * - **z+H**
     - ``"z+H"``
     - HN⁻¹O
     - Protonated z-ion (additional H from c-ion)
   * - **c-H**
     - ``"c-H"``
     - H₂N
     - Dehydrogenated c-ion (loss of H from c-ion)
   * - **Immonium**
     - ``"i"``
     - C⁻¹O⁻¹
     - Amino acid immonium ions (residue-specific)

Side-Chain Fragment Ions
--------------------------

These ions involve side-chain cleavage and are residue-specific:

.. list-table::
   :header-rows: 1
   :widths: 20 20 25 35

   * - Ion Type
     - String Value
     - Delta Formula
     - Description
   * - **d-ion**
     - ``"d"``
     - C₂H₃N
     - General side-chain loss from a-ion
   * - **d-ion (Val)**
     - ``"d-valine"``
     - C₃H₆N
     - Valine-specific d-ion
   * - **da-ion**
     - ``"da"``
     - 
     - General da-ion (side-chain loss)
   * - **da-ion (Thr)**
     - ``"da-threonine"``
     - C₂H₄NO
     - Threonine-specific da-ion
   * - **da-ion (Ile)**
     - ``"da-isoleucine"``
     - C₄H₈N
     - Isoleucine-specific da-ion
   * - **db-ion**
     - ``"db"``
     - 
     - General db-ion (side-chain loss)
   * - **db-ion (Thr)**
     - ``"db-threonine"``
     - C₃H₆N
     - Threonine-specific db-ion
   * - **db-ion (Ile)**
     - ``"db-isoleucine"``
     - C₃H₆N
     - Isoleucine-specific db-ion
   * - **v-ion**
     - ``"v"``
     - C₂H₂NO
     - Side-chain retention (complement of d-ion)
   * - **w-ion**
     - ``"w"``
     - C₃H₃O
     - General side-chain retention from x-ion
   * - **w-ion (Val)**
     - ``"w-valine"``
     - C₄H₆O
     - Valine-specific w-ion
   * - **wa-ion**
     - ``"wa"``
     - 
     - General wa-ion (side-chain retention)
   * - **wa-ion (Thr)**
     - ``"wa-threonine"``
     - C₃H₄O₂
     - Threonine-specific wa-ion
   * - **wa-ion (Ile)**
     - ``"wa-isoleucine"``
     - C₅H₈O
     - Isoleucine-specific wa-ion
   * - **wb-ion**
     - ``"wb"``
     - 
     - General wb-ion (side-chain retention)
   * - **wb-ion (Thr)**
     - ``"wb-threonine"``
     - C₄H₆O
     - Threonine-specific wb-ion
   * - **wb-ion (Ile)**
     - ``"wb-isoleucine"``
     - C₄H₆O
     - Isoleucine-specific wb-ion

Internal Fragment Ions
-----------------------

Internal fragments result from cleavage at both N- and C-termini. Formulas are based on the mzPAF specification.

.. list-table::
   :header-rows: 1
   :widths: 20 20 25 35

   * - Ion Type
     - String Value
     - Delta Formula
     - Description
   * - **by**
     - ``"by"``
     - (none)
     - N-terminal b-ion, C-terminal y-ion
   * - **ax**
     - ``"ax"``
     - H⁻²
     - N-terminal a-ion, C-terminal x-ion
   * - **cz**
     - ``"cz"``
     - (none)
     - N-terminal c-ion, C-terminal z-ion
   * - **ay**
     - ``"ay"``
     - C⁻¹O⁻¹
     - N-terminal a-ion, C-terminal y-ion
   * - **az**
     - ``"az"``
     - C⁻¹H⁻³N⁻¹O⁻¹
     - N-terminal a-ion, C-terminal z-ion
   * - **bx**
     - ``"bx"``
     - CH⁻²O
     - N-terminal b-ion, C-terminal x-ion
   * - **bz**
     - ``"bz"``
     - H⁻³N⁻¹
     - N-terminal b-ion, C-terminal z-ion
   * - **cx**
     - ``"cx"``
     - CHNO
     - N-terminal c-ion, C-terminal x-ion
   * - **cy**
     - ``"cy"``
     - H₃N
     - N-terminal c-ion, C-terminal y-ion
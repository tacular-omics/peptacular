Fragmentation presets
=====================

A preset names an activation method and gives back the keyword arguments
:func:`peptacular.fragment` already takes, so you do not have to write the ion series, the
neutral losses and the charge rule by hand.

.. testcode::

   import peptacular as pt

   print(pt.fragmentation_preset("hcd", charge_state=2))

.. testoutput::

   {'ion_types': ('b', 'y'), 'neutral_deltas': ('H2O', 'NH3'), 'max_ndeltas': 1, 'charges': (1,)}

Splat it into :func:`peptacular.fragment`:

.. testcode::

   fragments = pt.fragment("ELVISLIVESK", **pt.fragmentation_preset("etd", charge_state=2))
   print(len(fragments), sorted({f.ion_type.value for f in fragments}))

.. testoutput::

   55 ['c', 'z.']

A preset is a default, not a measurement
----------------------------------------

Each preset names the series its method usually gives for a tryptic-sized peptide, from the
literature cited below. The series you actually see depend on the instrument, the activation
energy, the precursor charge and the peptide, so treat a preset as a starting point and
override it when your data say otherwise.

The presets
-----------

Every preset sets ``ion_types``, ``neutral_deltas`` and ``max_ndeltas=1`` (at most one neutral
loss per ion). ``charge_state`` adds ``charges``, always ``1`` to ``max(1, charge_state - 1)``,
because the fragments of a ``z+`` precursor cannot each carry ``z`` charges. Leave
``charge_state`` out and :func:`peptacular.fragment` uses the annotation's own charge.

.. list-table::
   :header-rows: 1
   :widths: 12 30 16 42

   * - preset
     - ion_types
     - neutral_deltas
     - primary literature
   * - ``cid`` (alias ``cad``)
     - b, y
     - H2O, NH3
     - Roepstorff and Fohlman 1984, Biol Mass Spectrom 11:601 (10.1002/bms.1200111109), the
       letter that names a, b, c, x, y, z; Paizs and Suhai 2005, Mass Spectrom Rev 24:508-548
       (10.1002/mas.20024) for b and y dominating low-energy CID and for the water and ammonia
       losses
   * - ``hcd`` (aliases ``beam-type cid``, ``pasef``)
     - b, y
     - H2O, NH3
     - Olsen et al. 2007, Nat Methods 4:709-712 (10.1038/nmeth1060)
   * - ``etd``
     - c, z.
     - H2O, NH3
     - Syka et al. 2004, PNAS 101:9528-9533 (10.1073/pnas.0402700101)
   * - ``ecd``
     - c, z.
     - H2O, NH3
     - Zubarev et al. 1998, JACS 120:3265-3266 (10.1021/ja973478k)
   * - ``ethcd`` (aliases ``etcid``, ``etcad``)
     - b, y, c, z.
     - H2O, NH3
     - Frese et al. 2012, Anal Chem 84:9668-9673 (10.1021/ac3025366); Swaney et al. 2007, Anal
       Chem 79:477-485 (10.1021/ac061457f) for the b and y ions that collisional activation of
       the electron-transfer survivor adds
   * - ``ead``
     - a, b, c, x, y, z.
     - H2O, NH3
     - Baba et al. 2021, J Am Soc Mass Spectrom 32:1964-1975 (10.1021/jasms.0c00425)
   * - ``uvpd``
     - a, b, c, x, y, z
     - H2O
     - Shaw et al. 2013, JACS 135:12646-12651 (10.1021/ja4029654); Brodbelt 2020, Chem Rev
       120:3328-3380 (10.1021/acs.chemrev.9b00440)

``etd`` holds no b or y series on purpose. The b and y ions seen in an electron-transfer
experiment come from collisional activation of the charge-reduced survivor (Swaney et al. 2007,
and Han et al. 2007, Rapid Commun Mass Spectrom 21:1567-1573, 10.1002/rcm.2994), which is what
``ethcd`` is for.

:data:`peptacular.FRAGMENTATION_PRESET_CITATIONS` holds the same citations keyed by preset name,
so a script can print the reference beside its numbers.

Overriding one entry
--------------------

The preset is a plain dict and the caller's key wins:

.. testcode::

   hcd = pt.fragmentation_preset("hcd", charge_state=2)
   phospho = pt.fragment("ELVIS[Phospho]LIVESK", **{**hcd, "neutral_deltas": ("H3PO4",)})
   print(len(phospho), sorted({f.ion_type.value for f in phospho}))

.. testoutput::

   34 ['b', 'y']

A method that is not listed is a new key in :data:`peptacular.FRAGMENTATION_PRESETS`:

.. testcode::

   pt.FRAGMENTATION_PRESETS["my-method"] = {"ion_types": ("b", "y", "c"), "neutral_deltas": (), "max_ndeltas": 1}
   print(len(pt.fragment("ELVISLIVESK", **pt.fragmentation_preset("my-method", charge_state=2))))

.. testoutput::

   33

Using a preset with fast_fragment
---------------------------------

:func:`peptacular.fast_fragment` takes only ``ion_types``, ``charges`` and ``monoisotopic``, so a
whole preset cannot be splatted into it. Filter the preset first, and note that ``fast_fragment``
covers the six main series, so a preset carrying a radical series such as ``z.`` is rejected:

.. testcode::

   preset = pt.fragmentation_preset("hcd", charge_state=2)
   fast = {k: v for k, v in preset.items() if k in ("ion_types", "charges", "monoisotopic")}
   print(sorted(key[0].value for key in pt.fast_fragment("ELVISLIVESK", **fast)))

.. testoutput::

   ['b', 'y']

Mapping an instrument's activation type
---------------------------------------

peptacular reads no files, so it does not know what a spectrum was fragmented with. A reader does:
`spxtacular <https://spxtacular.readthedocs.io/>`_ carries the method on
``MsnSpectrum.activation_type``, normalised from the PSI-MS dissociation-method accessions, and its
names line up with the preset names above apart from the aliases in
:data:`peptacular.FRAGMENTATION_PRESET_ALIASES`. Map it there and pass the preset here.

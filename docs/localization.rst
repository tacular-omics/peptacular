Localisation isomers
====================

A modification whose position is not certain can be written three ways in ProForma.
:func:`peptacular.localization_isomers` expands each one into its concrete placements, and
:func:`peptacular.site_determining_ions` says which fragment ions tell the isomers apart.

Candidate positions come from the ProForma string alone. peptacular has no table of which
residues a modification can sit on, and does not look one up in tacular, UNIMOD or anywhere
else. Where a site list matters, say it yourself with :func:`peptacular.candidate_sites`.

Expanding ambiguity
-------------------

.. testcode::

   import peptacular as pt

   # A range: the mod goes on each residue inside it.
   print([a.serialize() for a in pt.localization_isomers("PEP(ST)[Phospho]IDE")])

   # An unknown position with no range: any residue, whatever its letter.
   print([a.serialize() for a in pt.localization_isomers("[Phospho]?PEST")])

   # ^2 places two copies on two different residues.
   print([a.serialize() for a in pt.localization_isomers("[Phospho]^2?STY")])

.. testoutput::

   ['PEPS[Phospho]TIDE', 'PEPST[Phospho]IDE']
   ['P[Phospho]EST', 'PE[Phospho]ST', 'PES[Phospho]T', 'PEST[Phospho]']
   ['S[Phospho]T[Phospho]Y', 'S[Phospho]TY[Phospho]', 'ST[Phospho]Y[Phospho]']

Groups and scores
~~~~~~~~~~~~~~~~~

A ``#label`` group lists the residues the mod may be on, optionally with a score each. The
group's label stays on the placed mod together with the score of the residue it landed on,
so the score lives in the mod string of each isomer:

.. testcode::

   isomers = pt.localization_isomers("PEPS[Phospho#g1(0.8)]T[#g1(0.2)]IDE")
   print([a.serialize() for a in isomers])

.. testoutput::

   ['PEPS[Phospho#g1(0.8)]TIDE', 'PEPST[Phospho#g1(0.2)]IDE']

A residue with no score gets the bare label (``Phospho#g1``). Ranges and unknown-position
mods carry no score. Cross-link (``#XL1``) and branch (``#BRANCH``) labels are not groups and
are left alone; a group label on a terminus, a range or an unknown-position mod raises
:class:`~peptacular.UnsupportedOperationError`.

Order, duplicates and limits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The output order is fixed: groups (by first residue), then ranges, then unknown-position
  mods, each placed from the N- to the C-terminus, the last one varying fastest.
- Isomers that come out identical, such as two copies of the same mod swapped between
  ranges, are returned once.
- ``max_isomers=`` (keyword-only) raises :class:`~peptacular.PeptacularError` as soon as
  the expansion passes that many isomers, before building the rest.
- An annotation with no ambiguity gives a one-item list.

.. testcode::

   try:
       pt.localization_isomers("[Phospho]^3?" + "S" * 40, max_isomers=1000)
   except pt.PeptacularError as error:
       print(type(error).__name__)

.. testoutput::

   PeptacularError

The same expansion is available as a method, ``annotation.localization_isomers()``.

Placing a mod on chosen residues
--------------------------------

:func:`peptacular.candidate_sites` puts one copy of a mod on each unmodified residue whose
letter you list. ``residues`` is required and has no default:

.. testcode::

   for position, isomer in pt.candidate_sites("PEPS[Phospho]TYK", "Phospho", residues="STY"):
       print(position, isomer.serialize())

.. testoutput::

   4 PEPS[Phospho]T[Phospho]YK
   5 PEPS[Phospho]TY[Phospho]K

Site-determining ions
---------------------

For each isomer, :func:`peptacular.site_determining_ions` returns the fragment ions whose
m/z is more than ``tolerance`` away from every ion of every other isomer, so a peak there is
evidence for that isomer alone. The ions are the :class:`~peptacular.Fragment` objects
``fragment()`` returns, computed on its fast prefix-sum path.

.. testcode::

   isomers = pt.localization_isomers("PEP(ST)[Phospho]IDE")
   ions = pt.site_determining_ions(isomers, ion_types=("b", "y"), charges=(1, 2), tolerance=10, unit="ppm")
   for isomer, frags in zip(isomers, ions):
       print(isomer.serialize(), [f"{f.ion_type}{f.position}+{f.charge_state} {f.mz:.3f}" for f in frags])

.. testoutput::

   PEPS[Phospho]TIDE ['b4+1 491.154', 'y4+1 477.219', 'b4+2 246.081', 'y4+2 239.113']
   PEPST[Phospho]IDE ['b4+1 411.187', 'y4+1 557.185', 'b4+2 206.097', 'y4+2 279.096']

``tolerance=None`` (the default) compares m/z values exactly, within 1e-6 Da. ``unit`` is
``"da"`` or ``"ppm"``. With one isomer every ion is returned.

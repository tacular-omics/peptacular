Interoperability
================

Peptacular can exchange peptide annotations with other Python proteomics
packages through optional adapters. Third-party packages are imported only
when their adapter is called, so the normal Peptacular installation remains
small.

Installation
------------

Install one integration or the complete phase-one set:

.. code-block:: bash

   pip install "peptacular[pyteomics]"
   pip install "peptacular[psm-utils]"
   pip install "peptacular[alphabase]"
   pip install "peptacular[interop]"

Pyteomics
---------

Pyteomics and Peptacular both model ProForma annotations, so mutually
supported features can be exchanged through ProForma text:

.. code-block:: python

   import peptacular as pt
   from peptacular.interop.pyteomics import from_pyteomics, to_pyteomics

   annotation = pt.parse("[Acetyl]-PEM[Oxidation]TIDE/2")
   pyteomics_value = to_pyteomics(annotation)
   restored = from_pyteomics(pyteomics_value)

   assert restored == annotation

Elemental compositions can also be exchanged with
``to_pyteomics_composition`` and ``from_pyteomics_composition``. A
``pyteomics.mass.Composition`` does not preserve ``ChargedFormula.charge``.
the outbound adapter therefore rejects charged formula values.

``psm_utils``
-------------

The ``psm_utils`` adapter connects Peptacular annotations to its
``Peptidoform`` and PSM I/O ecosystem:

.. code-block:: python

   import peptacular as pt
   from peptacular.interop.psm_utils import from_psm_utils, to_psm_utils

   annotation = pt.parse("PEM[Oxidation]TIDE/2")
   peptidoform = to_psm_utils(annotation)
   restored = from_psm_utils(peptidoform)

Only the peptidoform is converted. Spectrum identifiers, scores, proteins, and
other PSM metadata remain on the ``psm_utils.PSM`` object.

AlphaBase
---------

AlphaBase's native peptide representation is a pandas DataFrame containing
``sequence``, ``mods``, ``mod_sites``, and ``charge`` columns. The batch
adapter returns that representation directly and refines it with AlphaBase's
public API:

.. code-block:: python

   import peptacular as pt
   from alphabase.spectral_library.base import SpecLibBase
   from peptacular.interop.alphabase import (
       from_alphabase_dataframe,
       to_alphabase_dataframe,
   )

   annotations = [pt.parse("[Acetyl]-PEM[Oxidation]TIDE/2")]
   precursor_df = to_alphabase_dataframe(annotations)

   # The result can be assigned directly to an AlphaBase spectral library.
   library = SpecLibBase()
   library.precursor_df = precursor_df

   restored = from_alphabase_dataframe(precursor_df)

For one annotation, ``to_alphabase_row`` returns a plain row dictionary and
``from_alphabase_row`` accepts a mapping. These helpers do not
invent an AlphaBase peptide class. they expose the columns used by AlphaBase's
actual DataFrame model.

The AlphaBase representation cannot encode every ProForma feature. By default,
``to_alphabase_row`` and ``to_alphabase_dataframe`` raise
``InteropConversionError`` instead of silently losing
information. Callers can explicitly request warning or drop behavior:

.. code-block:: python

   from peptacular.interop import LossPolicy

   row = to_alphabase_row(annotation, loss_policy=LossPolicy.WARN)

Global isotope modifications, labile and unlocalized modifications, ambiguous
intervals or residues, annotation names, charge-adduct identities, and
localized modifications without a usable name are unsupported. Fixed
modifications are expanded onto matching residues before conversion.

Compatibility and errors
------------------------

``InteropConversionError`` identifies values that the target representation
cannot preserve. ``MissingOptionalDependencyError`` includes the exact extra
to install when an integration package is absent. Lossy AlphaBase conversion
uses ``LossyConversionWarning`` when its policy is ``warn``.

Pyteomics and ``psm_utils`` ultimately parse the serialized ProForma value.
Their supported subsets may differ from Peptacular's. When a target parser
rejects a feature, the adapter wraps its parser failure in an
``InteropConversionError`` while preserving the original exception as its
cause.

Release compatibility notes
---------------------------

The text adapters check a ProForma round trip before returning a target object.
If the target parser changes any annotation data, conversion raises
``InteropConversionError``. In particular, Pyteomics 5.0.1 changes negative
integer charges during serialization, so those outbound conversions are
rejected. The same guard applies to ``psm_utils``, which uses Pyteomics.

AlphaBase modification conversion currently supports single, unqualified
modification names registered in AlphaBase. CV accessions, explicit CV name
prefixes, scores, and additional modification tags require an explicit lossy
policy. Exact integral charges such as ``2``, ``2.0``, and ``"2"`` are accepted.
Fractional, boolean, and non-finite charges raise an error.

AlphaBase refines its tables by sorting on peptide length. Returned DataFrame
row order can therefore differ from input order. The inbound adapter preserves
the DataFrame's current order. Convert individual rows when you need to keep
an external association with the original input order.

The default installation includes none of these optional packages. Each extra
may install substantial transitive dependencies. Install only the integration
you need.

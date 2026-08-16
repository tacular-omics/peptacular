Stable JSON representation
==========================

Peptacular provides a lossless, versioned JSON representation for
``ProFormaAnnotation`` and the structured ProForma component model. The format
preserves modification counts, ambiguity scores, intervals, charge carriers,
cross-links, and compound peptidoform ions.

.. code-block:: python

   import peptacular as pt

   annotation = pt.ProFormaAnnotation.parse("PEM[Oxidation]TIDE/2")

   mapping = annotation.to_dict()
   restored = pt.ProFormaAnnotation.from_dict(mapping)

   text = annotation.to_json(indent=2)
   restored = pt.ProFormaAnnotation.from_json(text)

Every root document contains ``schema_version`` and ``$schema`` metadata.
Decoding is strict: unsupported versions, unknown object types, missing fields,
and extra fields are rejected. The decoder uses a closed type registry and never
imports a class selected by JSON input.

The bundled JSON Schema is available without filesystem assumptions:

.. code-block:: python

   schema = pt.get_proforma_json_schema()

The current schema version is ``1.0``. Additive or breaking representation
changes require a new schema version; existing decoders do not silently accept a
document from an unknown version.

Installation
============

Requirements
------------

* Python 3.12 or higher
* `tacular <https://tacular.readthedocs.io/>`_ (installed automatically) for amino acid,
  element and modification data

Install from PyPI
-----------------

.. code-block:: bash

   pip install peptacular

Optional extras
---------------

.. code-block:: bash

   pip install "peptacular[mcp]"         # local MCP server for AI agents, see the MCP page
   pip install "peptacular[pyteomics]"   # Pyteomics adapter
   pip install "peptacular[psm-utils]"   # psm_utils adapter
   pip install "peptacular[alphabase]"   # AlphaBase adapter
   pip install "peptacular[interop]"     # all three interoperability adapters

See :doc:`interoperability` and :doc:`mcp` for details.

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/tacular-omics/peptacular.git
   cd peptacular
   pip install -e .

Using uv (recommended for development)
--------------------------------------

.. code-block:: bash

   git clone https://github.com/tacular-omics/peptacular.git
   cd peptacular
   uv sync

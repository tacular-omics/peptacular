import os
import sys

from peptacular import __version__

sys.path.insert(0, os.path.abspath(".."))

project = "Peptacular"
copyright = "2024, Patrick Tyler Garrett"
author = "Patrick Tyler Garrett"
release = __version__
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # For Google/NumPy style docstrings
    "sphinx_autodoc_typehints",  # Uses your type hints automatically
    "sphinx.ext.viewcode",  # Adds source code links
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",  # Renders LaTeX math in .. math:: directives
    "myst_parser",
]

# Optional: configure source suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
# The ProForma JSON schema, plus the repo-root llms.txt and llms-full.txt
# (llmstxt.org), are served at the site root.
html_extra_path = [
    "../src/peptacular/schemas/proforma-json-v1.schema.json",
    "../llms.txt",
    "../llms-full.txt",
]

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "inherited-members": False,
    "show-inheritance": True,
}

autodoc_typehints = "description"  # Or 'signature' to put types in signature
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Examples that need pandas or polars are skipped when the library is not installed
# (``:skipif: pd is None``). peptacular itself depends on neither.
doctest_global_setup = """
try:
    import pandas as pd
except ImportError:
    pd = None
try:
    import polars as pl
except ImportError:
    pl = None
"""

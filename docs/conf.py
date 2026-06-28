import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "Peptacular"
copyright = "2024, Patrick Tyler Garrett"
author = "Patrick Tyler Garrett"

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
html_static_path = ["_static"]

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "inherited-members": False,
    "show-inheritance": True,
}

autodoc_typehints = "description"  # Or 'signature' to put types in signature
napoleon_google_docstring = True
napoleon_numpy_docstring = True

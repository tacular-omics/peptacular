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
try:
    import numpy as np
except ImportError:
    np = None
"""


def _hide_private_bases(app, name, obj, options, bases):
    """Drop private base classes (``ProFormaAnnotation``'s ``_ModAccessMixin``) from "Bases:"."""
    bases[:] = [base for base in bases if not base.__name__.startswith("_")] or [object]


def _mixin_source(app, modname):
    """Let viewcode link ``ProFormaAnnotation.<method>`` to methods defined on ``_ModAccessMixin``."""
    if modname != "peptacular.annotation._mod_access":
        return None
    from sphinx.pycode import ModuleAnalyzer

    analyzer = ModuleAnalyzer.for_module(modname)
    analyzer.find_tags()
    tags = dict(analyzer.tags)
    prefix = "_ModAccessMixin."
    tags.update({"ProFormaAnnotation." + name[len(prefix) :]: tag for name, tag in analyzer.tags.items() if name.startswith(prefix)})
    return analyzer.code, tags


def setup(app):
    app.connect("autodoc-process-bases", _hide_private_bases)
    app.connect("viewcode-find-source", _mixin_source)

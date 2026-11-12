"""Sphinx configuration for the SCI_py documentation.

Build locally from the repository root:

    pip install -r docs/requirements.txt
    sphinx-build -b html docs docs/_build/html

then open docs/_build/html/index.html.
"""
import os
import sys

# Make the package importable without installing it.
sys.path.insert(0, os.path.abspath("../SCI_py"))

project = "SCI_py"
author = "CompOmics"
copyright = "2026, " + author

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

# Docstrings are numpy style.
napoleon_google_docstring = False
napoleon_numpy_docstring = True

autosummary_generate = True
autosummary_imported_members = False
autodoc_default_options = {"members": True, "member-order": "bysource"}
autodoc_typehints = "description"
add_module_names = False

# Heavy or optional dependencies are stubbed so the docs build without them.
autodoc_mock_imports = ["anndata", "upsetplot", "sklearn", "umap", "statsmodels", "yaml", "toml"]

templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "furo"
html_title = "SCI_py"

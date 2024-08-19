# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/cas'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cas-tools'
copyright = '2024, Huseyin Kir, Ugur Bayindir'
author = 'Huseyin Kir, Ugur Bayindir'
release = '1.0.7'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    "sphinx.ext.autodoc",
    "sphinx.ext.githubpages",
    "sphinx_rtd_theme",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "myst_parser",
]

autodoc_member_order = 'bysource'
autodoc_typehints = 'description'

source_suffix = [".rst", ".md"]
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

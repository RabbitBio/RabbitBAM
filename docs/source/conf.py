# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sphinx_rtd_theme

project = 'RabbitBAM'
copyright = '2025, RabbitBio'
author = 'RabbitBio'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

html_theme_options = {
    "navigation_depth": 3,
}

extensions = [
    'sphinx.ext.autodoc',          # Automatically document APIs
    'sphinx.ext.napoleon',         # Support for NumPy/Google-style docstrings
    'sphinx.ext.viewcode',         # Add links to source code
    'sphinx.ext.todo',             # Support for TODO directives
    'sphinx.ext.autosummary',      # Generate summary tables for modules/classes
    'sphinx.ext.coverage',         # Coverage checks for documentation
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

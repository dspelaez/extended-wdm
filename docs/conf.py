# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import ewdm

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ewdm'
copyright = '2024, Daniel Pelaez-Zapata'
author = 'Daniel Pelaez-Zapata'
version = ewdm.__version__
release = ewdm.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx_gallery.gen_gallery'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# configure napoleon
napoleon_google_docstring = True
napoleon_include_init_with_doc = True

# gallery configuration
sphinx_gallery_conf = {
    'examples_dirs': 'gallery',    # path to your example scripts
    'gallery_dirs': 'auto_gallery' # path to where to save gallery generated output
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
# html_logo = "_static/logo_circled.png"
# html_theme_options = {
    # 'logo_only': True,
    # 'display_version': False,
# }

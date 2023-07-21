# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import subprocess
from datetime import date

# sys.path.insert(0, os.path.abspath('.'))
# sys.path.append(os.path.abspath('exts'))

tags = subprocess.getoutput('git describe --tags --abbrev=0')
print(tags)

# -- Project information -----------------------------------------------------
project = 'Q-POP'
copyright = f"{date.today().year}, DOE COMMS"
author = 'Xiaoxing Cheng'

# The full version, including alpha/beta/rc tags
release = 'v0.0.1'
# release = tags


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # 'breathe', 'exhale',
    'sphinxcontrib.bibtex', 'myst_parser', "sphinx_design", "sphinx_copybutton", "sphinxcontrib.mermaid", "sphinx_multiversion"]

source_suffix = ['.rst', '.md']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
html_sidebars = {
    'bd-sidebar__bottom': [
        'versioning.html',
    ],
}
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'node_modules', 'api']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'
html_theme_options = {
    "home_page_in_toc": True,
    # "github_url": "https://github.com/muprosoftware/muprosdk",
    "repository_url": "https://github.com/DOE-COMMS/Q-POP-Diffraction",
    "use_repository_button": True,
    "repository_branch": "main",
    "path_to_docs": "docs",
    "use_edit_page_button": True,
    "use_sidenotes": True,
    "use_issues_button": True,
    "use_source_button": True,
    "use_download_button": True,
    "icon_links_label": "Quick Links",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/DOE-COMMS",
            "icon": "fa-brands fa-github",
            "type": "fontawesome"
        }, {
            "name": "Email",
            "url": "mailto:lqc3@psu.edu",
            "icon": "fa-regular fa-envelope",
            "type": "fontawesome"
        },  {
            "name": "DOE COMMS",
            "url": "https://sites.psu.edu/doecomms/news/",
            "icon": "fa-regular fa-file",
            "type": "fontawesome"
        }
        # }, {
        #     "name": "PDF",
        #     "url": "https://download.muprosoftware.com/muprosdk-latest.pdf",
        #     "icon": "fa-regular fa-file-pdf",
        #     "type": "fontawesome"
        # }
    ]
}

# html_logo = "_static/mupro-leftright.png"
html_favicon = "_static/logo.png"
html_title = "Q-POP"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = ["custom.css"]

# Tell sphinx what the primary language being documented is.
# primary_domain = 'fortran'

# Tell sphinx what the pygments highlight language should be.
# highlight_language = 'fortran'
bibtex_bibfiles = ['refs.bib']

# ------ myst related settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

# myst_heading_anchors = 4
myst_all_links_external = True
myst_number_code_blocks = ["typescript"]
myst_heading_anchors = 2
myst_footnote_transition = True
myst_dmath_double_inline = True
# ------ myst related settings

# Setup the breathe extension
# breathe_projects = {
#     "MUPRO Phase-Field SDK": "./_doxygen/xml"
# }
# breathe_default_project = "MUPRO Phase-Field SDK"

# # Setup the exhale extension
# exhale_args = {
#     # These arguments are required
#     "containmentFolder":     "./api",
#     "rootFileName":          "library_root.rst",
#     "doxygenStripFromPath":  "..",
#     # Heavily encouraged optional argument (see docs)
#     "rootFileTitle":         "Library API",
#     # Suggested optional arguments
#     "createTreeView":        True,
#     # TIP: if using the sphinx-bootstrap-theme, you need
#     # "treeViewIsBootstrap": True,
#     "exhaleExecutesDoxygen": True,
#     "exhaleUseDoxyfile": True,
#     # "exhaleDoxygenStdin":    "INPUT = ../../dev"
# }

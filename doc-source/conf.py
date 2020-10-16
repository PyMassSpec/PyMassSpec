#!/usr/bin/env python3

# This file is managed by 'repo_helper'. Don't edit it directly.

# stdlib
import os
import re
import sys

sys.path.append(os.path.abspath('.'))
sys.path.append(os.path.abspath('..'))

# 3rd party
import ipynb2rst

# this package
from __pkginfo__ import __version__

nitpicky = True
exclude_patterns = ['../pyms-demo/old/', 'demo_rst/reading_andi.rst', 'chapter09.rst', 'chapter10.rst']

github_username = "domdfcoding"
github_repository = "PyMassSpec"
github_url = f"https://github.com/{github_username}/{github_repository}"

rst_prolog = f""".. |pkgname| replace:: PyMassSpec
.. |pkgname2| replace:: ``PyMassSpec``
.. |browse_github| replace:: `Browse the GitHub Repository <{github_url}>`__
"""

author = "PyMassSpec Authors"
project = "PyMassSpec"
slug = re.sub(r'\W+', '-', project.lower())
release = version = __version__
copyright = "2019-2020 Dominic Davis-Foster"  # pylint: disable=redefined-builtin
language = 'en'
package_root = "pyms"

extensions = [
		'sphinx_toolbox',
		'sphinx_toolbox.more_autodoc',
		'sphinx_toolbox.more_autosummary',
		'sphinx_toolbox.tweaks.param_dash',
		'sphinx.ext.intersphinx',
		'sphinx.ext.mathjax',
		'sphinxcontrib.httpdomain',
		'sphinxcontrib.extras_require',
		'sphinx.ext.todo',
		'sphinxemoji.sphinxemoji',
		'notfound.extension',
		'sphinx_copybutton',
		'sphinxcontrib.default_values',
		'sphinxcontrib.toctree_plus',
		'seed_intersphinx_mapping',
		'autodocsumm',
		'nbsphinx',
		'enum_tools.autoenum',
		'sphinx.ext.autosectionlabel',
		]

sphinxemoji_style = 'twemoji'
todo_include_todos = bool(os.environ.get("SHOW_TODOS", 0))
gitstamp_fmt = "%d %b %Y"

templates_path = ['_templates']
html_static_path = ['_static']
source_suffix = '.rst'
exclude_patterns = []

master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'

intersphinx_mapping = {
		'python': ('https://docs.python.org/3/', None),
		'sphinx': ('https://www.sphinx-doc.org/en/stable/', None),
		}

html_theme = 'domdf_sphinx_theme'
html_theme_options = {
		'includehidden': False,
		'logo_only': False,
		}
html_theme_path = ["../.."]
html_show_sourcelink = True  # True will show link to source

html_context = {
		'display_github': True,
		'github_user': 'domdfcoding',
		'github_repo': 'PyMassSpec',
		'github_version': 'master',
		'conf_py_path': '/doc-source/',
		}

htmlhelp_basename = slug

latex_documents = [('index', f'{slug}.tex', project, author, 'manual')]
man_pages = [('index', slug, project, [author], 1)]
texinfo_documents = [('index', slug, project, author, slug, project, 'Miscellaneous')]

toctree_plus_types = {
		"class",
		"function",
		"method",
		"data",
		"enum",
		"flag",
		"confval",
		"directive",
		"role",
		"confval",
		"protocol",
		"typeddict",
		"namedtuple",
		}

add_module_names = False
hide_none_rtype = True
all_typevars = True
overloads_location = "bottom"


autodoc_exclude_members = [   # Exclude "standard" methods.
		"__dict__",
		"__class__",
		"__dir__",
		"__weakref__",
		"__module__",
		"__annotations__",
		"__orig_bases__",
		"__parameters__",
		"__subclasshook__",
		"__init_subclass__",
		"__attrs_attrs__",
		"__init__",
		"__new__",
		"__getnewargs__",
		"__abstractmethods__",
		"__hash__",
		]
autodoc_default_options = {
		'members': None,  # Include all members (methods).
		'special-members': None,
		"autosummary": None,
		"show-inheritance": None,
		'exclude-members': ','.join(autodoc_exclude_members),
		}

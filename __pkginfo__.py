#  This file is managed by `git_helper`. Don't edit it directly
#  Copyright (C) 2020 Dominic Davis-Foster <dominic@davis-foster.co.uk>
#
#  This file is distributed under the same license terms as the program it came with.
#  There will probably be a file called LICEN[S/C]E in the same directory as this file.
#
#  In any case, this program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This script based on https://github.com/rocky/python-uncompyle6/blob/master/__pkginfo__.py
#

# stdlib
import pathlib

__all__ = [
		"__copyright__",
		"__version__",
		"modname",
		"pypi_name",
		"py_modules",
		"entry_points",
		"__license__",
		"__author__",
		"short_desc",
		"author",
		"author_email",
		"github_username",
		"web",
		"github_url",
		"project_urls",
		"repo_root",
		"long_description",
		"install_requires",
		"extras_require",
		"classifiers",
		"keywords",
		"import_name",
		]

__copyright__ = """
2019-2020 Dominic Davis-Foster <dominic@davis-foster.co.uk>
"""

__version__ = "2.2.21"

modname = "PyMassSpec"
pypi_name = "PyMassSpec"
import_name = "pyms"
py_modules = []
entry_points = {
		"console_scripts": []
		}

__license__ = "GNU General Public License v2 (GPLv2)"

short_desc = "Python Toolkit for Mass Spectrometry"

__author__ = author = "Dominic Davis-Foster"
author_email = "dominic@davis-foster.co.uk"
github_username = "domdfcoding"
web = github_url = f"https://github.com/domdfcoding/PyMassSpec"
project_urls = {
		"Documentation": f"https://PyMassSpec.readthedocs.io",
		"Issue Tracker": f"{github_url}/issues",
		"Source Code": github_url,
		}

repo_root = pathlib.Path(__file__).parent

# Get info from files; set: long_description
long_description = (repo_root / "README.rst").read_text(encoding="utf-8").replace("2.2.21", __version__) + '\n'
conda_description = """Python Toolkit for Mass Spectrometry


Before installing please ensure you have added the following channels: bioconda, conda-forge, domdfcoding"""
__all__.append("conda_description")

install_requires = (repo_root / "requirements.txt").read_text(encoding="utf-8").split('\n')
extras_require = {'all': []}

classifiers = [
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Developers',
		'Intended Audience :: Education',
		'Intended Audience :: End Users/Desktop',
		'Intended Audience :: Science/Research',
		'Operating System :: OS Independent',
		'Topic :: Education',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Chemistry',
		'Topic :: Software Development :: Libraries :: Python Modules',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: Implementation :: CPython',
		'Programming Language :: Python :: 3.7',
		'Programming Language :: Python :: 3.8',
		'Programming Language :: Python',
		'Programming Language :: Python :: 3 :: Only',
		'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',

		]

keywords = ""

long_description = long_description.replace(':ref:`here <pyms-demo>`.', '`here <pyms_demo_>`__.')
long_description = long_description.replace(':ref:`here <pyms-demo/data-files>`.', '`here <datafiles_>`__.')
long_description = long_description.replace(':ref:`Demos and Examples <pyms-demo>`', '`Demos and Examples <pyms_demo_>`__')
long_description = long_description.replace(':ref:`pyms-demo/20a/`', '`pyms-demo/20a/`')
long_description += '''

.. _pyms_demo: https://pymassspec.readthedocs.io/en/master/pyms-demo/introduction.html#pyms-demo
.. _datafiles: https://pymassspec.readthedocs.io/en/master/pyms-demo/data-files.html
'''

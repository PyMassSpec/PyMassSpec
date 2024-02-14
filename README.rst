************
PyMassSpec
************

.. start short_desc

**Python Toolkit for Mass Spectrometry**

.. end short_desc

.. start shields

.. list-table::
	:stub-columns: 1
	:widths: 10 90

	* - Docs
	  - |docs| |docs_check|
	* - Tests
	  - |actions_linux| |actions_windows| |actions_macos| |coveralls|
	* - PyPI
	  - |pypi-version| |supported-versions| |supported-implementations| |wheel|
	* - Anaconda
	  - |conda-version| |conda-platform|
	* - Activity
	  - |commits-latest| |commits-since| |maintained| |pypi-downloads|
	* - QA
	  - |codefactor| |actions_flake8| |actions_mypy|
	* - Other
	  - |license| |language| |requires|

.. |docs| image:: https://img.shields.io/readthedocs/pymassspec/latest?logo=read-the-docs
	:target: https://pymassspec.readthedocs.io/en/latest
	:alt: Documentation Build Status

.. |docs_check| image:: https://github.com/PyMassSpec/PyMassSpec/workflows/Docs%20Check/badge.svg
	:target: https://github.com/PyMassSpec/PyMassSpec/actions?query=workflow%3A%22Docs+Check%22
	:alt: Docs Check Status

.. |actions_linux| image:: https://github.com/PyMassSpec/PyMassSpec/workflows/Linux/badge.svg
	:target: https://github.com/PyMassSpec/PyMassSpec/actions?query=workflow%3A%22Linux%22
	:alt: Linux Test Status

.. |actions_windows| image:: https://github.com/PyMassSpec/PyMassSpec/workflows/Windows/badge.svg
	:target: https://github.com/PyMassSpec/PyMassSpec/actions?query=workflow%3A%22Windows%22
	:alt: Windows Test Status

.. |actions_macos| image:: https://github.com/PyMassSpec/PyMassSpec/workflows/macOS/badge.svg
	:target: https://github.com/PyMassSpec/PyMassSpec/actions?query=workflow%3A%22macOS%22
	:alt: macOS Test Status

.. |actions_flake8| image:: https://github.com/PyMassSpec/PyMassSpec/workflows/Flake8/badge.svg
	:target: https://github.com/PyMassSpec/PyMassSpec/actions?query=workflow%3A%22Flake8%22
	:alt: Flake8 Status

.. |actions_mypy| image:: https://github.com/PyMassSpec/PyMassSpec/workflows/mypy/badge.svg
	:target: https://github.com/PyMassSpec/PyMassSpec/actions?query=workflow%3A%22mypy%22
	:alt: mypy status

.. |requires| image:: https://dependency-dash.repo-helper.uk/github/PyMassSpec/PyMassSpec/badge.svg
	:target: https://dependency-dash.repo-helper.uk/github/PyMassSpec/PyMassSpec/
	:alt: Requirements Status

.. |coveralls| image:: https://img.shields.io/coveralls/github/PyMassSpec/PyMassSpec/master?logo=coveralls
	:target: https://coveralls.io/github/PyMassSpec/PyMassSpec?branch=master
	:alt: Coverage

.. |codefactor| image:: https://img.shields.io/codefactor/grade/github/PyMassSpec/PyMassSpec?logo=codefactor
	:target: https://www.codefactor.io/repository/github/PyMassSpec/PyMassSpec
	:alt: CodeFactor Grade

.. |pypi-version| image:: https://img.shields.io/pypi/v/PyMassSpec
	:target: https://pypi.org/project/PyMassSpec/
	:alt: PyPI - Package Version

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/PyMassSpec?logo=python&logoColor=white
	:target: https://pypi.org/project/PyMassSpec/
	:alt: PyPI - Supported Python Versions

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/PyMassSpec
	:target: https://pypi.org/project/PyMassSpec/
	:alt: PyPI - Supported Implementations

.. |wheel| image:: https://img.shields.io/pypi/wheel/PyMassSpec
	:target: https://pypi.org/project/PyMassSpec/
	:alt: PyPI - Wheel

.. |conda-version| image:: https://img.shields.io/conda/v/domdfcoding/PyMassSpec?logo=anaconda
	:target: https://anaconda.org/domdfcoding/PyMassSpec
	:alt: Conda - Package Version

.. |conda-platform| image:: https://img.shields.io/conda/pn/domdfcoding/PyMassSpec?label=conda%7Cplatform
	:target: https://anaconda.org/domdfcoding/PyMassSpec
	:alt: Conda - Platform

.. |license| image:: https://img.shields.io/github/license/PyMassSpec/PyMassSpec
	:target: https://github.com/PyMassSpec/PyMassSpec/blob/master/LICENSE
	:alt: License

.. |language| image:: https://img.shields.io/github/languages/top/PyMassSpec/PyMassSpec
	:alt: GitHub top language

.. |commits-since| image:: https://img.shields.io/github/commits-since/PyMassSpec/PyMassSpec/v2.4.1
	:target: https://github.com/PyMassSpec/PyMassSpec/pulse
	:alt: GitHub commits since tagged version

.. |commits-latest| image:: https://img.shields.io/github/last-commit/PyMassSpec/PyMassSpec
	:target: https://github.com/PyMassSpec/PyMassSpec/commit/master
	:alt: GitHub last commit

.. |maintained| image:: https://img.shields.io/maintenance/yes/2024
	:alt: Maintenance

.. |pypi-downloads| image:: https://img.shields.io/pypi/dm/PyMassSpec
	:target: https://pypi.org/project/PyMassSpec/
	:alt: PyPI - Downloads

.. end shields

PyMassSpec is a Python_ package for processing gas chromatography-mass spectrometry data.
PyMassSpec provides a framework and a set of components for rapid development and testing of methods for processing of chromatography--mass spectrometry data.
PyMassSpec can be used interactively through the Python shell, in a `Jupyter Notebook <https://jupyter.org/>`_, or the functions can be collected into scripts when it is preferable to perform data processing in the batch mode.

|

Forked from the original PyMS Repository: https://github.com/ma-bio21/pyms.
Originally by Andrew Isaac, Sean O'Callaghan and Vladimir Likić. The original publication can be found here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115

The original project seems to have been abandoned as there has been no activity since 2017.

|

.. contents:: Table of Contents
    :local:



The PyMassSpec project
=========================

The directory structure of PyMassSpec is as follows:

.. code-block:: text

    /
    ├── pyms:      The PyMassSpec code
    │
    ├── pyms-data: Example GC-MS data files
    │
    ├── pyms-demo: Examples of how to use PyMassSpec
    │
    ├── tests: pytest tests
    │
    └── doc-source: Sphinx source for documentation

Features
=========

Installation
==============

PyMassSpec can be installed with the following command:

.. code-block:: bash

    $ pip --user install PyMassSpec

This will also install the following dependencies:

.. code-block:: bash

    numpy >= 1.16.2
    scipy >= 1.2.1
    pymzml >= 2.2.1
    matplotlib >= 3.0.2
    openpyxl >= 2.6.2
    netCDF4 >= 1.5.0
    biopython >= 1.74
    deprecation >= 2.0.6


PyMassSpec can also make use of 'mpi4py' if it is installed. See https://mpi4py.readthedocs.io/en/stable/ for further information.


Usage
=======

A tutorial illustrating various PyMassSpec features in detail is provided
in subsequent chapters of this User Guide. The commands executed
interactively are grouped together by example, and can be found
`here <https://pymassspec.readthedocs.io/en/master/pyms-demo/introduction.html#pyms-demo>`__.

.. If you are viewing this source, the examples can be found in the pyms-demo directory, and the data files in pyms-data

The data used in the PyMassSpec documentation and examples is available
`here <https://pymassspec.readthedocs.io/en/master/pyms-demo/data-files.html>`__.

In the "`Demos and Examples`_" section there
is a page corresponding to each example, coded with the chapter number
(ie. "``pyms-demo/20a/``" corresponds to the Example 20a, from Chapter 2).

Each example has a script named 'proc.py' which contains the commands given in the example.
These scripts can be run with the following command:

.. code-block:: bash

    $ python3 proc.py

Example processing GC-MS data
-------------------------------

Download the file ``gc01_0812_066.jdx`` and save it in the folder ``data``.
This file contains GC-MS data in the the JCAMP-DX format.

First the raw data is loaded:

    >>> from pyms.GCMS.IO.JCAMP import JCAMP_reader
    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
    -> Reading JCAMP file 'Data/gc01_0812_066.jdx'
    >>> data
    <pyms.GCMS.Class.GCMS_data at 0x7f3ec77da0b8>

The intensity matrix object is then built by binning the data:

    >>> from pyms.IntensityMatrix import build_intensity_matrix_i
    >>> im = build_intensity_matrix_i(data)

In this example, we show how to obtain the dimensions of the
newly created intensity matrix, then loop over all ion chromatograms,
and for each ion chromatogram apply Savitzky-Golay noise filter
and tophat baseline correction:

    >>> n_scan, n_mz = im.size
    >>> from pyms.Noise.SavitzkyGolay import savitzky_golay
    >>> from pyms.TopHat import tophat
    >>> for ii in range(n_mz):
    ...     print("working on IC", ii)
    ...     ic = im.get_ic_at_index(ii)
    ...     ic1 = savitzky_golay(ic)
    ...     ic_smooth = savitzky_golay(ic1)
    ...     ic_base = tophat(ic_smooth, struct="1.5m")
    ...     im.set_ic_at_index(ii, ic_base)

The resulting noise and baseline corrected ion chromatogram is saved back into the intensity matrix.

Further examples can be found in the `documentation`_

Contributing
==============

Contributions are very welcome. Tests can be run with `pytest`_.
Please ensure the coverage is at least |coveralls|
before you submit a pull request.

For further information see the section `Contributing to PyMassSpec`_

License
=========
PyMassSpec is Free and Open Source software released under the `GNU General Public License version 2 <GPL_>`__.


Issues
========

If you encounter any problems, please `file an issue`_ along with a
detailed description.


.. _`documentation`: https://pymassspec.readthedocs.io
.. _`Contributing to PyMassSpec`: https://pymassspec.readthedocs.io/en/master/Contributing/Contributing.html
.. _`pytest`: https://pytest.org/
.. _`file an issue`: https://github.com/domdfcoding/pymassspec/issues
.. _Python: https://www.python.org/
.. _GPL: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
.. _Demos and Examples: https://pymassspec.readthedocs.io/en/master/pyms-demo/introduction.html#pyms-demo

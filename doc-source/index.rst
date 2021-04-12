===========
PyMassSpec
===========

.. start short_desc

.. documentation-summary::

.. end short_desc

.. start shields

.. only:: html

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

	.. |docs| rtfd-shield::
		:project: pymassspec
		:alt: Documentation Build Status

	.. |docs_check| actions-shield::
		:workflow: Docs Check
		:alt: Docs Check Status

	.. |actions_linux| actions-shield::
		:workflow: Linux
		:alt: Linux Test Status

	.. |actions_windows| actions-shield::
		:workflow: Windows
		:alt: Windows Test Status

	.. |actions_macos| actions-shield::
		:workflow: macOS
		:alt: macOS Test Status

	.. |actions_flake8| actions-shield::
		:workflow: Flake8
		:alt: Flake8 Status

	.. |actions_mypy| actions-shield::
		:workflow: mypy
		:alt: mypy status

	.. |requires| requires-io-shield::
		:alt: Requirements Status

	.. |coveralls| coveralls-shield::
		:alt: Coverage

	.. |codefactor| codefactor-shield::
		:alt: CodeFactor Grade

	.. |pypi-version| pypi-shield::
		:project: PyMassSpec
		:version:
		:alt: PyPI - Package Version

	.. |supported-versions| pypi-shield::
		:project: PyMassSpec
		:py-versions:
		:alt: PyPI - Supported Python Versions

	.. |supported-implementations| pypi-shield::
		:project: PyMassSpec
		:implementations:
		:alt: PyPI - Supported Implementations

	.. |wheel| pypi-shield::
		:project: PyMassSpec
		:wheel:
		:alt: PyPI - Wheel

	.. |conda-version| image:: https://img.shields.io/conda/v/domdfcoding/PyMassSpec?logo=anaconda
		:target: https://anaconda.org/domdfcoding/PyMassSpec
		:alt: Conda - Package Version

	.. |conda-platform| image:: https://img.shields.io/conda/pn/domdfcoding/PyMassSpec?label=conda%7Cplatform
		:target: https://anaconda.org/domdfcoding/PyMassSpec
		:alt: Conda - Platform

	.. |license| github-shield::
		:license:
		:alt: License

	.. |language| github-shield::
		:top-language:
		:alt: GitHub top language

	.. |commits-since| github-shield::
		:commits-since: v2.3.0
		:alt: GitHub commits since tagged version

	.. |commits-latest| github-shield::
		:last-commit:
		:alt: GitHub last commit

	.. |maintained| maintained-shield:: 2021
		:alt: Maintenance

	.. |pypi-downloads| pypi-shield::
		:project: PyMassSpec
		:downloads: month
		:alt: PyPI - Downloads

.. end shields

PyMassSpec is a Python_ package for processing gas chromatography-mass spectrometry data.
PyMassSpec provides a framework and a set of components for rapid development and testing of methods for processing of chromatography--mass spectrometry data.
PyMassSpec can be used interactively through the Python shell, in a `Jupyter Notebook <https://jupyter.org/>`_, or the functions can be collected into scripts when it is preferable to perform data processing in the batch mode.

|

Forked from the original PyMS Repository: https://github.com/ma-bio21/pyms.
Originally by Andrew Isaac, Sean O'Callaghan and Vladimir Likić. The original publication can be found here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115

The original project seems to have been abandoned as there has been no activity since 2017.


.. contents:: Table of Contents
	:local:



The PyMassSpec project
=========================

The directory structure of PyMassSpec is as follows:

.. code-block:: text

	/
	├── pyms:     The PyMassSpec code
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

.. start installation

.. installation:: PyMassSpec
	:pypi:
	:github:
	:anaconda:
	:conda-channels: bioconda, conda-forge, domdfcoding

.. end installation


.. PyMassSpec can also make use of 'mpi4py' if it is installed. See https://mpi4py.readthedocs.io/en/stable/ for further information.


Usage
=======

A tutorial illustrating various PyMassSpec features in detail is provided
in subsequent chapters of this User Guide. The commands executed
interactively are grouped together by example, and can be found
:ref:`here <pyms-demo>`.

.. If you are viewing this source, the examples can be found in the pyms-demo directory, and the data files in pyms-data

The data used in the PyMassSpec documentation and examples is available
:ref:`here <PyMassSpec test and example data files>`.

In the ":ref:`Demos and Examples <pyms-demo>`" section there
is a page corresponding to each example, coded with the chapter number
(ie. ":ref:`pyms-demo/20a/`" corresponds to the Example 20a, from Chapter 2).

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
	...	 print("working on IC", ii)
	...	 ic = im.get_ic_at_index(ii)
	...	 ic1 = savitzky_golay(ic)
	...	 ic_smooth = savitzky_golay(ic1)
	...	 ic_base = tophat(ic_smooth, struct="1.5m")
	...	 im.set_ic_at_index(ii, ic_base)

The resulting noise and baseline corrected ion chromatogram is saved back into the intensity matrix.

Further examples can be found in the `documentation`_

License
=========
PyMassSpec is Free and Open Source software released under the `GNU General Public License version 2 <GPL_>`__.


Issues
========

If you encounter any problems, please `file an issue`_ along with a
detailed description.


.. _`documentation`: https://pymassspec.readthedocs.io
.. _`pytest`: https://pytest.org/
.. _`file an issue`: https://github.com/domdfcoding/pymassspec/issues
.. _Python: https://www.python.org/
.. _GPL: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html


Installation
-------------

.. start installation

.. installation:: PyMassSpec
	:pypi:
	:github:
	:anaconda:
	:conda-channels: bioconda, conda-forge, domdfcoding

.. end installation


.. toctree::
	:hidden:

	Home<self>

.. toctree::
	:numbered:
	:caption: User Guide
	:maxdepth: 2

	10_gcms_raw_data_model
	20_gcms_data_derived_objects
	30_data_filtering
	40_peak_detection_and_representation
	50_peak_alignment_by_dynamic_programming
	60_display

.. toctree::
	:maxdepth: 2
	:caption: Documentation

	pyms/Base
	pyms/BillerBiemann
	pyms/Display
	pyms/DPA
	pyms/Experiment
	pyms/Gapfill
	pyms/GCMS
	pyms/IntensityMatrix
	pyms/IonChromatogram
	pyms/json
	pyms/Mixins
	pyms/Spectrum
	pyms/Noise
	pyms/Peak
	pyms/Simulator
	pyms/TopHat
	pyms/Utils
	changelog


.. toctree::
	:maxdepth: 2
	:caption: Demos and Examples

	pyms-demo/introduction
	pyms-demo/data-files

.. toctree::
	:hidden:

	pyms-demo/20e
	pyms-demo/32
	pyms-demo/55
	pyms-demo/56
	pyms-demo/64
	pyms-demo/70a
	pyms-demo/70b
	pyms-demo/70c
	pyms-demo/71
	pyms-demo/90
	pyms-demo/91
	pyms-demo/92
	pyms-demo/93
	pyms-demo/94
	pyms-demo/95
	pyms-demo/A1
	pyms-demo/A2
	pyms-demo/x10


.. toctree::
	:maxdepth: 3
	:caption: Contributing

	contributing
	Source
	StyleGuide


.. start links

.. only:: html

	View the :ref:`Function Index <genindex>` or browse the `Source Code <_modules/index.html>`__.

	`Browse the GitHub Repository <https://github.com/PyMassSpec/PyMassSpec>`__

.. end links

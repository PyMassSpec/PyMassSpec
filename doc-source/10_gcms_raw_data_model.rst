*********************
GC-MS Raw Data Model
*********************

.. contents:: Table of Contents
    :local:

Introduction
=============

PyMassSpec can read gas chromatography-mass spectrometry (GC-MS) data stored in
Analytical Data Interchange for Mass Spectrometry (ANDI-MS), [#ANDI-MS]_
and Joint Committee on Atomic and Molecular Physical Data (JCAMP-DX) [#JCAMP-DX]_
formats. The information contained in the data files can vary significantly
depending on the instrument, vendor's software, or conversion utility.
PyMassSpec makes the following assumptions about the information contained in the data file:

* The data contain the m/z and intensity value pairs across a scan.
* Each scan has a retention time.

Internally, PyMassSpec stores the raw data from ANDI files or JCAMP files as a
:class:`~pyms.GCMS.Class.GCMS_data` object.

.. include:: demo_rst/reading_jcamp.rst

.. note:: This example is in :file:`pyms-demo/jupyter/reading_jcamp.ipynb`. There is also an example in that directory for reading ANDI-MS files.

.. include:: demo_rst/comparing_datasets.rst

.. note:: This example is in :file:`pyms-demo/jupyter/comparing_datasets.ipynb`.


.. The file ``gc01_0812_066.jdx`` (located in "data") is a GC-MS experiment
.. converted from Agilent ChemStation format to JCAMP format using File
.. Translator Pro by ChemSW, Inc.


.. The file ``gc01_0812_066.cdf`` is a GC-MS experiment converted to ANDI-MS
.. format from Agilent ChemStation (the same data as in the example 20a above).

.. rubric:: Footnotes

.. [#ANDI-MS] ANDI-MS was developed by the Analytical Instrument Association
.. [#JCAMP-DX] JCAMP-DX is maintained by the International Union of Pure and Applied Chemistry

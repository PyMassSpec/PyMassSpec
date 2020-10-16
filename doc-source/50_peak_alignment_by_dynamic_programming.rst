***************************************
Peak alignment by dynamic programming
***************************************

.. contents:: Table of Contents
    :local:

PyMS provides functions to align GC-MS peaks by dynamic programming [1]_.
The peak alignment by dynamic programming uses both peak apex retention time and mass spectra.
This information is determined from the raw GC-MS data by applying a series of processing steps to produce data that can then be aligned and used for statistical analysis.
The details are described in this chapter.

Preparation of multiple experiments for peak alignment by dynamic programming
=============================================================================


.. include:: demo_rst/Experiment.rst

.. note:: This example is in `pyms-demo/jupyter/Experiment.ipynb`.

.. include:: demo_rst/Multiple_Experiments.rst

.. note:: This example is in `pyms-demo/jupyter/Multiple_Experiments.ipynb`.


Dynamic Programming Alignment
========================================================================

.. include:: demo_rst/DPA.rst

.. note:: These examples are in `pyms-demo/jupyter/DPA.ipynb`.

Common Ion Area Quantitation
================================
.. note:: This example is in :ref:`pyms-demo/64 <demo-64>`

The ``area.csv`` file produced in the preceding section lists the total area of each peak in the alignment.
The total area is the sum of the areas of each of the individual ions in the peak.
While this approach produces broadly accurate results, it can result in errors where neighbouring peaks or unfiltered noise add to the peak in some way.

One alternative to this approach is to pick a single ion which is common to a particular peak (compound), and to report only the area of this ion for each occurrence of that peak in the alignment.
Using the method :meth:`common_ion() <pyms.DPA.Alignment.Alignment.common_ion>` of the class :class:`~pyms.DPA.Alignment.Alignment`, ``PyMassSpec`` can select an ion for each aligned peak which is both abundant and occurs most often for that peak.
We call this the 'Common Ion Algorithm' (CIA).

To use this method it is essential that the individual ion
areas have been set (see section :ref:`individual_ion_areas`).


Using the Common Ion Algorithm
----------------------------------

When using the CIA for area quantitation, a different method of the class :class:`~pyms.DPA.Alignment.Alignment` is used to write the area matrix; :meth:`write_common_ion_csv() <pyms.DPA.Alignment.Alignment.write_common_ion_csv>`.
This requires a list of the common ions for each peak in the alignment.
This list is generated using the Alignment class method :meth:`common_ion() <pyms.DPA.Alignment.Alignment.common_ion>`.

Continuing from the previous example, the following invokes common ion filtering on previously created alignment object 'A9':

    >>> common_ion_list = A9.common_ion()

The variable 'common_ion_list' is a list of the common ion for each
peak in the alignment. This list is the same length as the
alignment. To write peak areas using common ion quantitation:

   >>> A9.write_common_ion_csv('output/area_common_ion.csv',common_ion_list)

References
============

.. [1] Robinson MD, De Souza DP, Keen WW, Saunders EC, McConville MJ, Speed TP, and Likic VA. A dynamic programming approach for the alignment of signal peaks in multiple gas chromatography-mass spectrometry experiments. `BMC Bioinformatics`, 8:419, 2007

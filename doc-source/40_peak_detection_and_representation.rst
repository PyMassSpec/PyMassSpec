***********************************
Peak detection and representation
***********************************

.. contents:: Table of Contents
    :local:

.. include:: demo_rst/Peak.rst

.. note:: This example is in :file:`pyms-demo/jupyter/Peak.ipynb`.


Peak Detection
================

The general use of a :class:`~pyms.Peak.Class.Peak` object is to extract them from the GC-MS data and build a list of peaks. In |pkgname|, the function for peak detection is based on the method of Biller and Biemann (1974) [1]_.
The basic process is to find all maximising ions in a pre-set window of scans, for a given scan.
The ions that maximise at a given scan are taken to belong to the same peak.

The function is :py:meth:`BillerBiemann() <pyms.BillerBiemann.BillerBiemann>`. in :mod:`pyms.BillerBiemann`.
The function has parameters for the window width for detecting the local maxima (``points``), and the number of ``scans`` across which neighbouring, apexing, ions are combined and considered as belonging to the same peak.
The number of neighbouring scans to combine is related to the likelihood of detecting a peak apex at a single scan or several neighbouring scans.
This is more likely when there are many scans across the peak.
It is also possible, however, when there are very few scans across the peak.
The scans are combined by taking all apexing ions to have occurred at the scan that had to greatest TIC prior to combining scans.

.. include:: demo_rst/Peak_Detection.rst

.. note:: This example is in :file:`pyms-demo/jupyter/Peak_Detection.ipynb`.

.. include:: demo_rst/Peak_Filtering_Noise_Analysis.rst

.. note:: This example is in :file:`pyms-demo/jupyter/Peak_Filtering_Noise_Analysis.ipynb`.

.. include:: demo_rst/Peak_Area_Estimation.rst

.. note:: This example is in :file:`pyms-demo/jupyter/Peak_Area_Estimation.ipynb`.


Individual Ion Areas
------------------------

.. note:: This example is in :ref:`pyms-demo/56 <demo-56>`

While the previous approach uses the sum of all areas in the peak to estimate the peak area, the user may also choose to record the area of each individual ion in each peak.

This can be useful when the intention is to later perform quantitation based on the area of a single characteristic ion for a particular compound.
It is also essential if using the Common Ion Algorithm for quantitation, outlined in the section :ref:`Common Ion Area Quantitation`.

To set the area of each ion for each peak, the following code is used:

    >>> from pyms.Peak.Function import peak_top_ion_areas
    >>> for peak in peak_list:
    ...     area_dict = peak_top_ions_areas(intensity_matrix, peak)
    ...     peak.set_ion_areas(area_dict)
    ...

This will set the areas of the 5 most abundant ions in each peak.
If it is desired to record more than the top five ions, the argument ``num_ions=x`` should be supplied, where ``x`` is the number of most abundant ions to be recorded.
For example:

.. code-block:: python

    ...     area_dict = peak_top_ions_areas(intensity_matrix, peak, num_ions=10)

will record the 10 most abundant ions for each peak.

The individual ion areas can be set instead of, or in addition to the total area for each peak.

Reading the area of a single ion in a peak
-------------------------------------------

If the individual ion areas have been set for a peak, it is possible to read the area of an individual ion for the peak.
For example:

>>> peak.get_ion_area(101)

will return the area of the :math:`m/z` value 101 for the peak.
If the area of that ion has not been set (i.e. it was not one of the most abundant ions), the function will return ``None``.

References
============

.. [1] Biller JE and Biemann K. Reconstructed mass spectra, a novel approach for the utilization of gas chromatograph–mass spectrometer data. *Anal. Lett.*, 7:515–528, 1974

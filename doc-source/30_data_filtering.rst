****************
Data Filtering
****************

.. contents:: Table of Contents

Introduction
=============

In this chapter filtering techniques that allow pre-processing of GC-MS data for analysis and comparison to other pre-processed GC-MS data are covered.

Time strings
==============

Before considering the filtering techniques, the mechanism for representing retention times is outlined here.

A time string is the specification of a time interval, that takes the format ``NUMBERs`` or ``NUMBERm`` for time interval in seconds or minutes. For example, these are valid time strings: ``10s`` (10 seconds) and ``0.2m`` (0.2 minutes).

.. include:: demo_rst/IntensityMatrix_Resizing.rst

.. note:: This example is in `pyms-demo/jupyter/IntensityMatrix_Resizing.ipynb`.

|

.. include:: demo_rst/NoiseSmoothing.rst

.. note:: This example is in `pyms-demo/jupyter/NoiseSmoothing.ipynb`.


.. include:: demo_rst/BaselineCorrection.rst

.. note:: This example is in `pyms-demo/jupyter/BaselineCorrection.ipynb`.

.. include:: demo_rst/IntensityMatrix_Preprocessing.rst

The resulting :class:`~pyms.IntensityMatrix.IntensityMatrix` object can be "dumped" to a file for later
retrieval. There are general perpose object file handling methods in
:py:meth:`pyms.Utils.IO <pyms.Utils.IO>`. For example;

    >>> from pyms.Utils.IO import dump_object
    >>> dump_object(im, "output/im-proc.dump")

.. note:: This example is in `pyms-demo/jupyter/IntensityMatrix_Preprocessing.ipynb`.


References
============

.. [1] Serra J. `Image Analysis and Mathematical Morphology`. Academic Press, Inc, Orlando, 1983. ISBN 0126372403

.. [2] Sauve AC and Speed TP. Normalization, baseline correction and alignment of high-throughput mass spectrometry data. `Procedings Gensips`, 2004

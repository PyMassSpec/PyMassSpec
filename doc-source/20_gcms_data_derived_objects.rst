 ****************************
GC-MS data derived objects
****************************

.. contents:: Table of Contents
    :local:

In the raw GC-MS data, consecutive scans do not necessarily contain the same
mass per charge (mass) values. For data processing, it is often necessary to
convert the data to a matrix with a set number of masses and scans.
In ``PyMassSpec`` the resulting object is called intensity matrix.
In this chapter the methods for converting the raw GC-MS data to an
intensity matrix object are illustrated.

IntensityMatrix Object
=======================

The general scheme for converting raw mass values is to bin intensity values
based on the interval the corresponding mass belongs to. The general procedure
is as follows:

* Set the interval between bins, lower and upper bin boundaries.

* Calculate the number of bins to cover the range of all masses.

* Centre the first bin at the minimum mass found for all the raw data.

* Sum intensities whose masses are in a given bin.

A mass, :math:`m`, is considered to belong to a bin when :math:`c - l \le m < c + u`,
where :math:`c` is the centre of the bin, :math:`l` is the lower boundary and :math:`u` is
the upper boundary of the bin. The default bin interval is one with a lower
and upper boundary of :math:`\pm0.5`.

A function to bin masses to the nearest integer is also available. The default
bin interval is one with a lower boundary of :math:`-0.3` and upper boundary of
:math:`+0.7` (as per the NIST library).

Discussion of Binning Boundaries
------------------------------------

For any chemical element :math:`X`, let :math:`w(x)` be the atomic weight of :math:`X`, and

:math:`\delta(X) = \frac{w(X) - \{w(X)\}}{w(X)}`

where :math:`\{a\}` is the integer value of :math:`a` (rounded to the nearest integer).

For example, for hydrogen :math:`\delta(^1\rm{H}) = \frac{1.007825032 - 1}{1.007825032} = 0.0076`.
Similarly :math:`\delta(^{12}\rm{C}) = 0`, :math:`\delta(^{14}\rm{N}) = 0.00022`,
:math:`\delta(^{16}\rm{O}) = -0.00032`, etc.

Let also :math:`\Delta(X) = w(X) - \{w(x)\}`. Then :math:`-0.023 <\Delta(^{31}\rm{P}),
\Delta(^{28}\rm{Si}) < 0`.

Let a compound undergo GC-MS and let Y be one of it's fragments. If Y consists
of :math:`k_{1}`atoms of type :math:`X_{1}`,
:math:`k_{2}` atoms of type :math:`X_{2}`,.....,
:math:`k_{r}` atoms of type :math:`X_{r}`, then
:math:`\Delta(Y) = k_{1}*\Delta(X_{1}) + k_{2}*\Delta(X_{2}) + ....+ k_{r}*\Delta(X_{r})`.

The fragment will usually not contain more than 2 or 3 P or Si atoms and if it's
molecular weight is less than 550 it may not contain more than 35 O atoms, so
:math:`\Delta(Y) \geq -0.023*5 - 0.00051*35 = -0.133`.

On the other hand, of Y contains :math:`k` H atoms and :math:`m` N atoms, then :math:`\Delta(Y)
\leq k*0.00783 + m*0.00051`. Since for each two hydrogen atoms at least one
carbon (or heavier) atom is needed, giving the limit of no more than 80 hydrogen
atoms. Therefore in this case (i.e. H and C atoms only) :math:`\Delta(Y) \leq 80*0.00783 = 0.63`.
If carbon is replaced by any heavier atom, at least 2 hydrogen atoms
will be eliminated and :math:`\Delta(Y)` will become even smaller.

If the molecular weight of :math:`Y` does not exceed 550 (typically the largest mass
scanned for in a GC-MS setup) then :math:`\mathbf{-0.133 \leq \Delta(Y) \leq 0.63}`.
This means that if we set our binning boundaries to :math:`(-0.3, 0.7)` or :math:`(-0.2, 0.8)`
the opportunity for having a fragment whose molecular weight is very close to
the boundary is minimised.

Since the resolution of MS is at least 0.1 dalton, we may assume that it's error
does not exceed 0.05, and MS accuracy will not cause additional problems.

.. include:: demo_rst/IntensityMatrix.rst

.. note:: This example is in `pyms-demo/jupyter/IntensityMatrix.ipynb`.

.. include:: demo_rst/MassSpectrum.rst

.. note:: This example is in `pyms-demo/jupyter/MassSpectrum.ipynb`.

.. include:: demo_rst/IonChromatogram.rst

.. note:: This example is in `pyms-demo/jupyter/IonChromatogram.ipynb`.

Writing IonChromatogram object to a file
--------------------------------------------

.. note:: This example is in :ref:`pyms-demo/31 <demo-31>`

The method :meth:`write() <pyms.IonChromatogram.IonChromatogram.write>`
of an :class:`~pyms.IonChromatogram.IonChromatogram`
object allows the ion chromatogram to be saved to a file:

    >>> tic.write("output/tic.dat", minutes=True)
    >>> im.get_ic_at_mass(73).write("output/ic.dat", minutes=True)

The flag ``minutes=True`` indicates that retention time will be saved in minutes.
The ion chromatogram object saved with with the
:meth:`write() <pyms.IonChromatogram.IonChromatogram.write>`
method is a plain ASCII file which contains a pair of
(retention time, intensity) per line.

.. code-block:: bash

    $ head tic.dat
    5.0930 2.222021e+07
    5.0993 2.212489e+07
    5.1056 2.208650e+07
    5.1118 2.208815e+07
    5.1181 2.200635e+07
    5.1243 2.200326e+07
    5.1306 2.202363e+07
    5.1368 2.198357e+07
    5.1431 2.197408e+07
    5.1493 2.193351e+07

Saving data
=============

.. note:: This example is in :ref:`pyms-demo/32 <demo-32>`

A matrix of intensity values can be saved to a file with the function
:meth:`save_data() <pyms.Utils.IO.save_data>`
from :meth:`pyms.Utils.IO <pyms.Utils.IO>`. A matrix of intensity values can
be returned from an :class:`~pyms.IntensityMatrix.IntensityMatrix`
with the method
:attr:`~pyms.IntensityMatrix.IntensityMatrix.intensity_array`.
For example,

    >>> from pyms.Utils.IO import save_data
    >>> mat = im.intensity_array
    array([[22128.,     0., 10221., ...,     0.,   470.,     0.],
           [22040.,     0., 10335., ...,   408.,     0.,   404.],
           [21320.,     0., 10133., ...,   492.,     0.,   422.],
           ...,
           [    0.,     0.,     0., ...,     0.,     0.,     0.],
           [    0.,     0.,     0., ...,     0.,     0.,     0.],
           [    0.,     0.,     0., ...,     0.,     0.,     0.]])
    >>> save_data("output/im.dat", mat)

It is also possible to save the list of masses (from
:attr:`im.mass_list <pyms.IntensityMatrix.IntensityMatrix.mass_list>`
and the list of retention times (from
:attr:`im.time_list <pyms.IntensityMatrix.IntensityMatrix.time_list>`
using the :meth:`save_data() <pyms.Utils.IO.save_data>` function.
For convenience, the intensity values, mass list and time list,
can be saved with the method
:meth:`export_ascii() <pyms.IntensityMatrix.IntensityMatrix.export_ascii>`.
For example,

    >>> im.export_ascii("output/data")

will create ``data.im.dat``, ``data.rt.dat`` and ``data.mz.dat``, where these
are the intensity matrix, retention time vector, and :math:`m/z` vector. By default
the data is saved as space separated data with a ".dat" extension. It is
also possible to save the data as comma separated data with a ".csv"
extension with the command:

    >>> im.export_ascii("output/data", "csv")

Additionally, the entire :class:`IntensityMatrix <pyms.IntensityMatrix.IntensityMatrix>`
can be exported to LECO CSV format.

    >>> im.export_leco_csv("output/data_leco.csv")

This facility is useful for import into other analytical software packages. The format has a header line specifying the column heading information as:

.. code-block:: text

    scan, retention time, mass1, mass2, ...

and then each row as the intensity data.


Importing ASCII data
====================

.. note:: This example is in :ref:`pyms-demo/32 <demo-32>`

The LECO CSV format data can be imported directly into an `~pyms.IntensityMatrix.IntensityMatrix` object. The data must follow the format outlined above. For example, the file saved above can be read and compared to the original:

.. code-block:: python

    >>> from pyms.IntensityMatrix import IntensityMatrix
    >>> iim = IntensityMatrix([0],[0],[[0]])
    >>> iim.import_leco_csv("output/data_leco.csv")
    >>> im.size
    >>> iim.size

The line :meth:`IntensityMatrix([0],[0],[[0]]) <pyms.IntensityMatrix.IntensityMatrix>` is required to create an empty :class:`~pyms.IntensityMatrix.IntensityMatrix` object.

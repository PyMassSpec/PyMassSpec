===============
Changelog
===============

Changes in v2.2.22-beta2
--------------------------

* :class:`pyms.Spectrum.Scan` and :class:`pyms.Spectrum.MassSpectrum` can now accept any values for ``mass`` and ``intensity`` that that can be converted to a :class:`float` or :class:`int`. This includes strings representing numbers. Previously only :class:`int` and :class:`float` values were permitted.

* If the mass and intensity values supplied to a :class:`pyms.Spectrum.Scan` or a :class:`pyms.Spectrum.MassSpectrum` are :class:`float`, :class:`int`, or a data type derived from :class:`numpy.number`, the data is stored in that type. For other data types, such as strings, :class:`decimal.Decimal` etc., the data is stored as :class:`float`.

  If the data contains values in mixed types then, in most cases, all values will be converted to :class:`float`. If you wish to control this behaviour you should construct a :class:`numpy.ndarray` with the desired type. See https://numpy.org/devdocs/user/basics.types.html for a list of types.

* A :exc:`TypeError` is no longer raised when creating a :class:`pyms.Spectrum.Scan` or a :class:`pyms.Spectrum.MassSpectrum` with a :class:`float`, :class:`int` etc. rather than a sequence. Instead, value is treated as being the sole element in a list.

* Passing a non-numeric string or a list of non-numeric strings to :class:`pyms.Spectrum.Scan` or :class:`pyms.Spectrum.MassSpectrum` now raises a :exc:`ValueError` and not a :exc:`TypeError` as in previous versions.

* :meth:`pyms.Peak.Class.Peak.ion_areas` now accepts dictionary keys as :class:`float` as well as :class:`int`.


Changes in v2.2.22-beta1
--------------------------

* :func:`~pyms.GCMS.IO.ANDI.ANDI_reader` and :class:`pyms.Spectrum.Scan` were modified to allow ANDI-MS files to be read if the data either:

    - had the *m/z* data stored from highest *m/z* to lowest; or
    - contained 0-length scans.

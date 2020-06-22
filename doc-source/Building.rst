*********************
Building from source
*********************

To build the ``PyMassSpec`` package from source using ``setuptools``, run the following command:

.. code-block:: bash

    $ python3 setup.py sdist bdist_wheel

``setuptools`` is configured using the file  :download:`setup.py <../setup.py>`.


Different formats are available for built distributions

.. list-table::
   :header-rows: 1

   * - Format
     - Description
     - Notes
   * - ``gztar``
     - gzipped tar file (``.tar.gz``)
     - default on Unix
   * - ``bztar``
     - bzipped tar file (``.tar.bz2``)
     -
   * - ``xztar``
     - bzipped tar file (``.tar.bz2``)
     -
   * - ``tar``
     - tar file (``.tar``)
     -
   * - ``zip``
     - zip file (``.zip``)
     - default on Windows
   * - ``wininst``
     - self-extracting ZIP file for Windows
     -
   * - ``msi``
     - Microsoft Installer
     -


|

**setup.py**


.. literalinclude:: ../setup.py
    :language: python
    :linenos:

**__pkginfo__.py**


.. literalinclude:: ../__pkginfo__.py
    :language: python
    :linenos:
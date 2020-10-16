*****************************
PyMassSpec coding Style Guide
*****************************

.. contents:: Table of Contents
    :local:


This document provides specific style conventions for PyMassSpec.
It should be read in conjunction with :PEP:`8` "Style Guide for Python Code",
by Guido van Rossum and Barry Warsaw

General
=========

Grouping commands and using newlines
-------------------------------------

Sort functions and class methods alphabetically, with dunder methods at the top.

Return copy.copy or copy.deepcopy only when this will not impact performance or otherwise absolutely necessary. Alternatively, use numpy.array().tolist().



Organise commands into logical groups, and separate if necessary with newlines to
improve readability.


Example:

.. code-block:: python

    # -- snip --
    if not isinstance(file_name, str):
        raise TypeError("'file_name' must be a string")

    try:
        file = CDF(file_name)
        self.__file_name = file_name
        self.__file_handle = file
    except CDFError:
        error("Cannot open file '%s'" % file_name)

    print(" -> Processing netCDF file '%s'" % (self.__file_name))

    self.__set_min_max_mass(file)
    self.__set_intensity_list(file)
    # -- snip --

In block statements (such as for loops and if statements), do not use the blank
line in a single group of statements; use one blank line to separate if the
block contains more than one group of statements.

Examples:

.. code-block:: python

    # -- snip --
    td_list = []
    for ii in range(len(time_list)-1):
        td = time_list[ii+1]-time_list[ii]
        td_list.append(td)
    # -- snip --

.. code-block:: python

    # -- snip ---
    if len(time_list) > len(intensity_matrix):

        self.set_scan_index()
        scan_index_list = self.__scan_index_list

        count = 0
        while len(intensity_matrix) < len(time_list):
            count = count + 1
            scan = numpy.repeat([0], max_mass - min_mass + 1)
            intensity_matrix.insert(0,scan)
    # -- snip ---

File pointers
---------------
Use ``fp`` for file pointer variables. If simultaneous use of two or more file
pointers is required, use ``fp1``, ``fp2``, etc.

Example:

.. code-block:: python

    fp1 = open('some_file.txt','w')
    fp2 = open('another.txt','w')


Short Comments
---------------

If a comment is short, the period at the end is best omitted. Longer comments of
block comments generally consist of one or more paragraphs built out of complete
sentences, and each sentence should end with a period.

Imports
=========

Grouping
----------

Group imports as:

#. Standard library imports
#. External module imports
#. Other PyMassSpec subpackage imports
#. This subpackage imports

Separate each group by a blank line.

Import forms
-------------

For standard library modules, always import the entire module name space. i.e.

.. code-block:: python

      import os
      ...
      os.path()



Naming Styles
===============

Variable names
----------------

Global variable names should be prefixed with an underscore to prevent their
export from the module.

For Specific variable names:

    - Use ``file_name`` instead of ``filename``
    - Use ``fp`` for file pointer, i.e.

        .. code-block:: python

            fp = open(file_name, 'r')

Module names
-------------
Module names should be short, starting with an uppercase letter (i.e. Utils.py).

Class names
------------
Class names use the CapWords convention. Classes for internal use have a leading
underscore in addition.

Exception Names
-----------------
Exceptions should be handled via the function
:py:meth:`pyms.Utils.Error.error() <pyms.Utils.Error.error>`.

Function Names
----------------
Function names should be lowercase, with words separated by underscores where
suitable to improve readability.

Method Names
------------------
Method names should follow the same principles as the function names.

Internal methods and instance variables
-----------------------------------------
Use one leading underscore only for internal methods and instance variables
which are not intended to be part of the class's public interface.

Class-private names
----------------------
Use two leading underscores to denote class-private names, this includes
class-private methods (eg. ``__privfunc()``).

.. note:: Python "mangles" these names with the class name:
    if class Foo has an attribute named ``__a``, it cannot be accessed by ``Foo.__a``.
    (it still could be accessed by calling ``Foo._Foo__a``.)

Private/public class attributes
---------------------------------
Public attributes should have no leading or trailing underscores. Private
attributes should have two leading underscores, no trailing underscores.
Non-public attributes should have a single leading underscore, no trailing
underscores (the difference between private and non-public is that the
former will never be useful for a derived class, while the latter might be).

Reminder: Python names with specific meanings
------------------------------------------------
* ``_single_leading_underscore``: weak "internal use" indicator (e.g. "``from M import *``" does not import objects whose name starts with an underscore).

* ``single_trailing_underscore_``: used by convention to avoid conflicts with Python keyword, "``Tkinter.Toplevel(master, class_='ClassName')``".

* ``__double_leading_underscore``: class-private names as of Python 1.4.

* ``__double_leading_and_trailing_underscore__``: "magic" objects or attributes that live in user-controlled namespaces, e.g. ``__init__``, ``__import__`` or ``__file__``.

Docstrings
===========

General
---------

* All sub-packages, modules, functions, and classes must have proper Sphinx docstrings

* When designating types for :type and :rtype, use the official names from the 'types' package i.e. ``BooleanType``, ``StringType``, ``FileType`` etc.

* All docstrings must start with a single summary sentence concisely describing the function, and this sentence must not be terminated by a period. Additional description may follow in the form of multi-sentenced paragraphs, separated by a blank line from the summary sentence - Leave one blank line above and below the docstring

* Separate ``:summary``, ``:param``/``:type``, ``:return``/``:rtype``, ``:author`` strings with one blank line

Packages
---------
Package doctrings are defined in ``__init__.py``. This example shows top three lines of ``pyms.__input__.py``:

Example:

.. code-block:: python

      """
      The root of the package pyms
      """

Modules
---------
A summary for the module should be written concisely in a single sentence, enclosed above and below with lines containing only ``"""``

Example:

.. code-block:: python

      """
      Provides general I/O functions
      """

Functions
----------

In all functions the following Sphinx tags must be defined:

    * ``:param``
    * ``:type`` (for all input arguments)
    * ``:return``
    * ``:rtype`` (unless the function returns None)
    * ``:author``

Other fields are optional.


Example:

.. code-block:: python

      def open_for_reading(file_name):

          """
          Opens file for reading, returns file pointer

          :param file_name: Name of the file to be opened for reading
          :type file_name: StringType

          :return: Pointer to the opened file
          :rtype: FileType

          :author: Jake Blues
          """

Classes
---------
* The root class docstring must contain ``:summary`` and ``:author`` fields

* The ``__init__`` method must contain ``:param`` and ``:type`` fields. Other fields are optional.

* Methods docstrings adhere to rules for Functions. Except for special methods (i.e. ``__len__()``, ``__del__()``, etc) which should contain only the ``:summary`` field, and possibly the ``:author`` field.

* Class methods. The rules for functions apply, except that the tag ``:author`` does not need to be defined (if authors are given in the class docstring).

    Examples:

    .. code-block:: python

        class ChemStation:

            """
            ANDI-MS reader for Agilent ChemStation NetCDF files

            :author: Jake Blues
            """

            def __init__(self, file_name):
                """
                :param file_name: The name of the ANDI-MS file
                :type file_name: StringType
                """

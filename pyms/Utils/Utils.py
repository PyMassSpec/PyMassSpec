"""
General utility functions
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019 Dominic Davis-Foster                                   #
#                                                                              #
#    This program is free software; you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License version 2 as         #
#    published by the Free Software Foundation.                                #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program; if not, write to the Free Software               #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                              #
################################################################################

# stdlib
import os
import pathlib
from collections import Sequence
from decimal import Decimal

# 3rd party
from typing import Any

import numpy

_list_types = (Sequence, numpy.core.ndarray)
_path_types = (str, os.PathLike, pathlib.Path)


def is_path(obj):
	"""
	Returns whether the object represents a filesystem path

	:param obj:
	:type obj:

	:return:
	:rtype:
	"""

	if isinstance(obj, _path_types):
		return True
	else:
		return hasattr(obj, " __fspath__")


def is_sequence(obj):
	"""
	Returns whether the object is a :class:`~collections.abc.Sequence`, and not a string

	:param obj:
	:type obj:

	:return:
	:rtype: bool
	"""

	return isinstance(obj, _list_types) and not isinstance(obj, str)


def is_sequence_of(obj: Any, of: Any) -> bool:
    """
	Returns whether the object is a :class:`~collections.abc.Sequence`, and not a string, of the given type

	:param obj:
	:type obj: any
	:param of:
	:type of: any

	:return:
	:rtype: bool
	"""

	return (
			isinstance(obj, _list_types)
			and not isinstance(obj, str)
			and all(isinstance(x, of) for x in obj)
	)


def is_positive_int(arg: str) -> bool:
    """
    Determines if the argument is an integer greater than zero

	:param arg: A string to be evaluate as a positive integer
	:type arg: types.str

	:return: A boolean indicator True or False
	:rtype:  bool

	:author: Milica Ng
	"""

	if not isinstance(arg, int):
		return False
	elif not (arg > 0):
		return False
	else:
		return True


def is_list_of_dec_nums(arg: str) -> bool:
    """
    Determines if the argument is a list of decimal numbers

	:param arg: A string to be evaluate as a list of decimal numbers
	:type arg: str

	:return: A boolean indicator True or False
	:rtype:  bool

	:author: Milica Ng
	"""

	return is_sequence_of(arg, (float, Decimal))

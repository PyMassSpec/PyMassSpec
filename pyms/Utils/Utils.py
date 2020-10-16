"""
General utility functions
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                              #
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
from typing import TYPE_CHECKING, Any, Sequence

# 3rd party
import numpy  # type: ignore

__all__ = ["is_path", "is_sequence", "is_sequence_of", "_number_types", "signedinteger", "is_number"]

if TYPE_CHECKING:
	signedinteger = int
else:
	signedinteger = numpy.signedinteger

_list_types = (Sequence, numpy.core.ndarray)
_path_types = (str, os.PathLike, pathlib.Path)
_number_types = (int, float, signedinteger)


def is_path(obj: Any) -> bool:
	"""
	Returns whether the object represents a filesystem path

	:param obj:
	"""

	if isinstance(obj, _path_types):
		return True
	else:
		return hasattr(obj, " __fspath__")


def is_sequence(obj) -> bool:
	"""
	Returns whether the object is a :class:`~collections.abc.Sequence`,
	and not a string.

	:param obj:
	"""

	return isinstance(obj, _list_types) and not isinstance(obj, str)


def is_sequence_of(obj: Any, of: Any) -> bool:
	"""
	Returns whether the object is a :class:`~collections.abc.Sequence`,
	and not a string, of the given type.

	:param obj:
	:param of:
	"""

	return isinstance(obj, _list_types) and not isinstance(obj, str) and all(isinstance(x, of) for x in obj)


def is_number(obj: Any) -> bool:

	return isinstance(obj, _number_types)

"""
Mixins for PyMassSpec Classes.
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
import math
from typing import List, Optional
from warnings import warn

# 3rd party
import numpy

# this package
from pyms.Utils.Utils import is_number

__all__ = [
		"MaxMinMassMixin",
		"MassListMixin",
		"TimeListMixin",
		"IntensityArrayMixin",
		"GetIndexTimeMixin",
		]


class MaxMinMassMixin:
	"""
	Mixin class to add the ``min_mass`` and ``max_mass`` properties,
	which provide read-only access to the internal
	``_min_mass`` and ``_max_mass`` attributes.
	"""  # noqa: D400

	_min_mass: Optional[float]
	_max_mass: Optional[float]

	@property
	def min_mass(self) -> Optional[float]:
		"""
		Returns the minimum *m/z* value in the spectrum.

		:author: Andrew Isaac
		"""

		return self._min_mass

	@property
	def max_mass(self) -> Optional[float]:
		"""
		Returns the maximum *m/z* value in the spectrum.

		:author: Andrew Isaac
		"""

		return self._max_mass


class MassListMixin(MaxMinMassMixin):
	"""
	Mixin class to add the ``mass_list`` property,
	which returns a copy of the internal ``_mass_list`` attribute.
	"""  # noqa: D400

	_mass_list: List[float]

	@property
	def mass_list(self) -> List[float]:
		"""
		Returns a list of the masses.

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
		"""

		return self._mass_list[:]


class TimeListMixin:
	"""
	Mixin class to add the ``time_list`` property,
	which returns a copy of the internal ``_time_list`` attribute.
	"""  # noqa: D400

	_time_list: List[float]

	@property
	def time_list(self) -> List[float]:
		"""
		Returns a copy of the time list.

		:return: List of retention times

		:authors: Andrew Isaac, Lewis Lee, Vladimir Likic
		"""

		return self._time_list[:]


class IntensityArrayMixin:
	_intensity_array: numpy.ndarray

	@property
	def intensity_array(self) -> numpy.ndarray:
		"""
		Returns a copy of the intensity array.

		:return: Matrix of intensity values.

		:authors: Andrew Isaac, Lewis Lee
		"""

		return numpy.copy(self._intensity_array)

	@property
	def intensity_matrix(self) -> numpy.ndarray:
		"""
		Returns a copy of the intensity matrix.

		:return: Matrix of intensity values.

		:author: Andrew Isaac
		"""

		warn(f"Use 'intensity_array' attribute instead", DeprecationWarning)

		return numpy.copy(self._intensity_array)

	@property
	def intensity_array_list(self) -> List[List[float]]:
		"""
		Returns a copy of the intensity array as a list of lists of floats.

		:return: Matrix of intensity values.

		:author: Andrew Isaac
		"""

		return self._intensity_array.tolist()

	@property
	def matrix_list(self) -> numpy.ndarray:
		"""
		Returns the intensity matrix as a list of lists of floats.

		:return: Matrix of intensity values

		:author: Andrew Isaac
		"""
		warn(f"Use 'intensity_array_list' attribute instead", DeprecationWarning)
		return self.intensity_array


class GetIndexTimeMixin:
	_min_rt: float
	_max_rt: float
	_time_list: List[float]

	def get_index_at_time(self, time: float) -> int:
		"""
		Returns the nearest index corresponding to the given time.

		:param time: Time in seconds

		:return: Nearest index corresponding to given time

		:authors: Lewis Lee, Tim Erwin, Vladimir Likic

		.. versionchanged:: 2.3.0

			Now returns ``-1`` if no index is found.
		"""

		if not is_number(time):
			raise TypeError("'time' must be a number")

		if (time < self._min_rt) or (time > self._max_rt):
			raise IndexError(
					f"time {time:.2f} is out of bounds (min: {self._min_rt:.2f}, max: {self._max_rt:.2f})"
					)

		time_list = self._time_list
		time_diff_min = self._max_rt
		ix_match = -1

		for ix in range(len(time_list)):

			time_diff = math.fabs(time - time_list[ix])

			if time_diff < time_diff_min:
				ix_match = ix
				time_diff_min = time_diff

		return ix_match

	def get_time_at_index(self, ix: int) -> float:
		"""
		Returns time at given index.

		:param ix:

		:authors: Lewis Lee, Vladimir Likic
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an integer")

		if ix < 0 or ix > len(self._time_list) - 1:
			raise IndexError("index out of bounds")

		return self._time_list[ix]

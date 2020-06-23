"""
Mixins for pyms Classes
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
from numbers import Number
from typing import List
from warnings import warn

# 3rd party
import deprecation  # type: ignore
import numpy  # type: ignore

# this package
from pyms import __version__


class MaxMinMassMixin:
	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'max_mass' attribute instead")
	def get_max_mass(self) -> float:
		"""
		Get the max mass value

		:return: The maximum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.max_mass

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'min_mass' attribute instead")
	def get_min_mass(self) -> float:
		"""
		Get the min mass value over all scans

		:return: The minimum mass of all the data
		:rtype: float

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.min_mass

	@property
	def min_mass(self) -> float:
		"""
		Returns the minimum m/z value in the spectrum

		:return: Minimum m/z
		:rtype: float

		:author: Andrew Isaac
		"""

		return self._min_mass

	@property
	def max_mass(self) -> float:
		"""
		Returns the maximum m/z value in the spectrum

		:return: Maximum m/z
		:rtype: float

		:author: Andrew Isaac
		"""

		return self._max_mass


class MassListMixin(MaxMinMassMixin):

	@property
	def mass_list(self) -> List:
		"""
		Returns a list of the masses

		:return: mass list
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self._mass_list[:]

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'mass_list' attribute instead")
	def get_mass_list(self) -> List:
		"""
		Returns a list of the masses

		:return: mass list
		:rtype: list

		:author: Qiao Wang
		:author: Andrew Isaac
		:author: Vladimir Likic
		"""

		return self.mass_list


class TimeListMixin:

	@property
	def time_list(self) -> List:
		"""
		Returns a copy of the time list

		:return: List of retention times
		:rtype: list

		:author: Andrew Isaac
		:author: Lewis Lee
		:author: Vladimir Likic
		"""

		return self._time_list[:]

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'time_list' attribute instead")
	def get_time_list(self) -> List:
		"""
		Returns a copy of the time list

		:return: List of retention times
		:rtype: list

		:author: Andrew Isaac
		"""

		return self.time_list


class IntensityArrayMixin:

	@property
	def intensity_array(self) -> List:
		"""
		Returns a copy of the intensity array

		:return: Matrix of intensity values
		:rtype: list

		:author: Andrew Isaac
		:author: Lewis Lee
		"""

		return numpy.copy(self._intensity_array)

	@property
	def intensity_matrix(self) -> List:
		"""
		Returns a copy of the intensity matrix

		:return: Matrix of intensity values
		:rtype: list

		:author: Andrew Isaac
		"""

		warn(f"Use 'intensity_array' attribute instead", DeprecationWarning)

		return numpy.copy(self._intensity_array)

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use 'intensity_array' attribute instead")
	def get_intensity_array(self) -> numpy.array:
		"""
		Returns the entire intensity array

		:return: Intensity array
		:rtype: numpy.ndarray

		:author: Lewis Lee
		:author: Vladimir Likic
		"""

		return self.intensity_array

	@property
	def intensity_array_list(self) -> List:
		"""
		Returns a copy of the intensity array as a list of lists of floats

		:return: Matrix of intensity values
		:rtype: list

		:author: Andrew Isaac
		"""

		return self._intensity_array.tolist()

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details=f"Use 'matrix_list' attribute instead")
	def get_matrix_list(self) -> List:
		"""
		Returns a copy of the intensity matrix as a list of lists of floats

		:return: Matrix of intensity values
		:rtype: list

		:author: Andrew Isaac
		"""

		return self.intensity_array

	@property
	def matrix_list(self) -> List:
		"""
		Returns a the intensity matrix as a list of lists of floats

		:return: Matrix of intensity values
		:rtype: list

		:author: Andrew Isaac
		"""
		warn(f"Use 'intensity_array' attribute instead", DeprecationWarning)
		return self.intensity_array


class GetIndexTimeMixin:
	def get_index_at_time(self, time: float) -> int:
		"""
		Returns the nearest index corresponding to the given time

		:param time: Time in seconds
		:type time: float

		:return: Nearest index corresponding to given time
		:rtype: int

		:author: Lewis Lee
		:author: Tim Erwin
		:author: Vladimir Likic
		"""

		if not isinstance(time, Number):
			raise TypeError("'time' must be a number")

		if (time < self._min_rt) or (time > self._max_rt):
			raise IndexError(f"time {time:.2f} is out of bounds (min: {self._min_rt:.2f}, max: {self._max_rt:.2f})")

		time_list = self._time_list
		time_diff_min = self._max_rt
		ix_match = None

		for ix in range(len(time_list)):

			time_diff = math.fabs(time - time_list[ix])

			if time_diff < time_diff_min:
				ix_match = ix
				time_diff_min = time_diff

		return ix_match

	def get_time_at_index(self, ix: int) -> float:
		"""
		Returns time at given index

		:param ix: An index
		:type ix: int

		:return: Time value
		:rtype: float

		:author: Lewis Lee
		:author: Vladimir Likic
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an integer")

		if ix < 0 or ix > len(self._time_list) - 1:
			raise IndexError("index out of bounds")

		return self._time_list[ix]

"""
Classes to model a GC-MS Ion Chromatogram
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
import copy
import pathlib
import warnings
from typing import Any, List, Optional, Sequence, Union

# 3rd party
import numpy  # type: ignore

# this package
from pyms.Base import pymsBaseClass
from pyms.Mixins import GetIndexTimeMixin, IntensityArrayMixin, TimeListMixin
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Utils import is_number, is_path, is_sequence

__all__ = ["IonChromatogram"]


class IonChromatogram(pymsBaseClass, TimeListMixin, IntensityArrayMixin, GetIndexTimeMixin):
	r"""
	Models an ion chromatogram

	An ion chromatogram is a set of intensities as a function of retention time.
	This can can be either *m/z* channel intensities (for example, ion
	chromatograms at ``m/z = 65``\), or cumulative intensities over all measured *m/z*.
	In the latter case the ion chromatogram is total ion chromatogram (TIC).

	The nature of an IonChromatogram object can be revealed by inspecting
	the value of the attribute 'mass'. This is set to the *m/z* value of the
	ion chromatogram, or to :py:obj:`None` for TIC.

	:param ia: Ion chromatogram intensity values
	:param time_list: A list of ion chromatogram retention times
	:param mass: Mass of ion chromatogram (:py:obj:`None` if TIC)

	:authors: Lewis Lee, Vladimir Likic, Dominic Davis-Foster (type assertions and properties)
	"""

	def __init__(self, ia: numpy.ndarray, time_list: List[float], mass: Optional[float] = None):
		if not isinstance(ia, numpy.ndarray):
			raise TypeError("'ia' must be a numpy array")

		if not is_sequence(time_list) or not all(is_number(time) for time in time_list):
			raise TypeError("'time_list' must be a list of numbers")

		if len(ia) != len(time_list):
			raise ValueError("Intensity array and time list differ in length")

		if mass is not None and not is_number(mass):
			raise TypeError("'mass' must be a number")

		self._intensity_array = ia
		self._time_list = time_list
		self._mass: Optional[float] = mass
		self._time_step = self.__calc_time_step()
		self._min_rt = min(time_list)
		self._max_rt = max(time_list)

	def __len__(self) -> int:
		"""
		Returns the length of the IonChromatogram object

		:return: Length of ion chromatogram
		:rtype: int

		:authors: Lewis Lee, Vladimir Likic
		"""

		return self._intensity_array.size

	def __sub__(self, other: "IonChromatogram") -> "IonChromatogram":
		"""
		Subtracts another IC from the current one

		:param other: Another IC
		"""

		ia_for_sub = other.intensity_array

		for i in range(self._intensity_array.size):
			self._intensity_array[i] = self._intensity_array[i] - ia_for_sub[i]

		return self

	def __eq__(self, other: Any) -> bool:
		"""
		Return whether this IonChromatogram object is equal to another object

		:param other: The other object to test equality with
		"""

		if isinstance(other, self.__class__):
			return (
					self.time_list == other.time_list
					and all(numpy.equal(self.intensity_array, other.intensity_array)) and self.mass == other.mass
					)

		return NotImplemented

	def __copy__(self) -> "IonChromatogram":
		"""
		Returns a new IonChromatogram containing a copy of the data in this object
		"""

		return IonChromatogram(
				ia=numpy.copy(self._intensity_array),
				time_list=self._time_list[:],
				mass=copy.copy(self._mass),
				)

	def __deepcopy__(self, memodict={}):
		return self.__copy__()

	def get_intensity_at_index(self, ix: int) -> float:
		"""
		Returns intensity at given index

		:param ix: An index

		:return: Intensity value

		:authors: Lewis Lee, Vladimir Likic
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an integer")

		if ix < 0 or ix > self._intensity_array.size - 1:
			raise IndexError("index out of bounds")

		return self._intensity_array[ix]

	@IntensityArrayMixin.intensity_array.setter  # type: ignore
	def intensity_array(self, ia: Union[Sequence, numpy.ndarray]):
		"""
		Sets the value for the intensity array

		:param ia: An array of new intensity values

		:author: Vladimir Likic
		"""

		if not is_sequence(ia):
			raise TypeError("'intensity_array' must be a Sequence")

		if not isinstance(ia, numpy.ndarray):
			ia = numpy.array(ia)

		self._intensity_array = ia

	def is_tic(self) -> bool:
		"""
		Returns whether the ion chromatogram is a total ion chromatogram (TIC)

		:authors: Lewis Lee, Vladimir Likic
		"""

		return self._mass is None

	@property
	def mass(self) -> Optional[float]:
		"""
		Returns the m/z channel of the IC

		:return: m/z channel of the IC

		:author: Sean O'Callaghan
		"""

		if self._mass is None:
			warnings.warn("TIC has no m/z label", Warning)

		return self._mass

	@property
	def time_step(self) -> float:
		"""
		Returns the time step

		:return: Time step

		:authors: Lewis Lee, Vladimir Likic
		"""

		return self._time_step

	def __calc_time_step(self) -> float:
		"""
		Calculates the time step

		:return: Time step value

		:authors: Lewis Lee, Vladimir Likic
		"""

		td_list = []
		for ii in range(len(self._time_list) - 1):
			td = self._time_list[ii + 1] - self._time_list[ii]
			td_list.append(td)

		td_array = numpy.array(td_list)
		time_step = td_array.mean()

		return time_step

	def write(self, file_name: Union[str, pathlib.Path], minutes: bool = False, formatting: bool = True):
		"""
		Writes the ion chromatogram to the specified file

		:param file_name: The name of the output file
		:param minutes: A boolean value indicating whether to write time in minutes
		:param formatting: Whether to format the numbers in the output.

		:authors: Lewis Lee, Vladimir Likic, Dominic Davis-Foster (pathlib support)
		"""

		if not is_path(file_name):
			raise TypeError("'file_name' must be a string or a PathLike object")

		file_name = prepare_filepath(file_name)

		with file_name.open("w") as fp:

			time_list = copy.deepcopy(self._time_list)

			if minutes:
				for ii in range(len(time_list)):
					time_list[ii] = time_list[ii] / 60.0

			for ii in range(len(time_list)):
				if formatting:
					fp.write(f"{time_list[ii]:8.4f} {self._intensity_array[ii]:#.6e}\n")
				else:
					fp.write(f"{time_list[ii]} {self._intensity_array[ii]}\n")

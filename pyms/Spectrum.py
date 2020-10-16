"""
Classes to model Mass Spectra and Scans
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
import pathlib
import re
import warnings
from typing import Any, Iterator, List, Optional, Sequence, Tuple, Union

# 3rd party
import numpy  # type: ignore

# this package
from pyms.Base import pymsBaseClass
from pyms.Mixins import MassListMixin
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.jcamp import xydata_tags
from pyms.Utils.Utils import is_path, is_sequence

__all__ = ["array_as_numeric", "Scan", "MassSpectrum", "normalize_mass_spec"]


def array_as_numeric(array: Union[Sequence, numpy.ndarray]) -> numpy.ndarray:
	"""
	Convert the given numpy array to a numeric data type.

	If the data in the array is already in a numeric data type no changes will be made.

	If ``array`` is a python :class:`~python:collections.abc.Sequence` then it will first be
	converted to a numpy array.

	:param array:
	"""

	if not isinstance(array, numpy.ndarray):
		array = numpy.array(array)

	if not numpy.issubdtype(array.dtype, numpy.number):
		# Array doesn't contain numerical data
		array = array.astype(numpy.float64)

	return array


class Scan(pymsBaseClass, MassListMixin):
	"""
	Generic object for a single Scan's raw data

	:param mass_list: A sequence of mass values
	:param intensity_list: A sequence intensity values

	:authors: Andrew Isaac, Qiao Wang, Vladimir Likic, Dominic Davis-Foster
	"""

	def __init__(
			self,
			mass_list: Union[Sequence[float], numpy.ndarray],
			intensity_list: Union[Sequence[float], numpy.ndarray],
			):
		mass_list = list(array_as_numeric(mass_list))
		intensity_list = list(array_as_numeric(intensity_list))

		if len(mass_list) != len(intensity_list):
			raise ValueError("'mass_list' is not the same size as 'intensity_list'")

		sorted_mass_list = sorted(mass_list)

		if sorted_mass_list != mass_list:
			# Mass list isn't in ascending order
			if sorted_mass_list[::-1] == mass_list:
				# Mass list is in descending order
				mass_list = mass_list[::-1]
				intensity_list = intensity_list[::-1]
			else:
				warnings.warn(
						"""Unknown sort order for mass list; it doesn't appear to be in either ascending or descending order.
Please report this at https://github.com/domdfcoding/pymassspec/issues and upload an example data file if possible.
"""
						)

		self._mass_list = mass_list
		self._intensity_list = intensity_list

		if self:
			self._min_mass = min(mass_list)
			self._max_mass = max(mass_list)
		else:
			self._min_mass = None
			self._max_mass = None

	def __len__(self) -> int:
		"""
		Returns the length of the object

		:authors: Andrew Isaac, Qiao Wang, Vladimir Likic
		"""

		return len(self._mass_list)

	def __bool__(self) -> bool:
		return bool(self._mass_list)

	def __eq__(self, other: Any) -> bool:
		"""
		Return whether this object is equal to another object

		:param other: The other object to test equality with
		"""

		if isinstance(other, self.__class__):
			return (
					numpy.array_equal(self._intensity_list, other.intensity_list)
					and numpy.array_equal(self._mass_list, other.mass_list)
					)

		return NotImplemented

	def __copy__(self) -> "Scan":
		"""
		Returns a copy of the object
		"""

		return self.__class__(self._mass_list[:], self._intensity_list[:])

	def __deepcopy__(self, memodict={}) -> "Scan":
		return self.__copy__()

	@property
	def __dict__(self):
		return {
				"intensity_list": self.intensity_list,
				"mass_list": self.mass_list,
				}

	def __iter__(self):
		for key, value in self.__dict__.items():
			yield key, value

	def __getstate__(self):
		return self.__dict__

	def __setstate__(self, state):
		self.__init__(**state)  # type: ignore

	def iter_peaks(self) -> Iterator[Tuple[float, float]]:
		"""
		Iterate over the peaks in the mass spectrum
		"""

		for mz, intensity in zip(self.mass_list, self.intensity_list):
			yield mz, intensity

	@property
	def intensity_list(self) -> List:
		"""
		Returns a copy of the intensity list

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
		"""

		return self._intensity_list[:]

	@property
	def mass_spec(self) -> List:
		"""
		Returns the intensity list

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
		"""

		return self._intensity_list

	@classmethod
	def from_dict(cls, dictionary):
		return cls(**dictionary)


class MassSpectrum(Scan):
	"""
	Models a binned mass spectrum

	:param mass_list: mass values
	:param intensity_list: intensity values

	:authors: Andrew Isaac, Qiao Wang, Vladimir Likic, Dominic Davis-Foster
	"""

	def __init__(
			self,
			mass_list: Union[Sequence[float], numpy.ndarray],
			intensity_list: Union[Sequence[float], numpy.ndarray],
			):
		Scan.__init__(self, mass_list, intensity_list)

	@Scan.intensity_list.setter  # type: ignore
	def intensity_list(self, value: List[float]):
		"""
		Set the intensity values for the spectrum

		:param value: list of intensity value for each mass in ``mass_list``
		"""

		value = array_as_numeric(value)

		# if not isinstance(value, _list_types) or not isinstance(value[0], Number):
		# 	raise TypeError("'intensity_list' must be a list of numbers")

		# if not len(self.mass_list) == len(value):
		# 	raise ValueError("'mass_list' and 'intensity_list' are not the same size")

		self._intensity_list = list(value)

	@Scan.mass_spec.setter  # type: ignore
	def mass_spec(self, value: List[float]):
		"""
		Set the intensity values for the spectrum

		:param value: list of intensity value for each mass in `mass_list`
		"""

		value = array_as_numeric(value)

		# if not isinstance(value, _list_types) or not isinstance(value[0], Number):
		# 	raise TypeError("'intensity_list' must be a list of numbers")

		# if not len(self.mass_list) == len(value):
		# 	raise ValueError("'mass_list' and 'intensity_list' are not the same size")

		self._intensity_list = list(value)

	@MassListMixin.mass_list.setter  # type: ignore
	def mass_list(self, value: List[float]):
		"""
		Set the mass values for the spectrum

		:param value: list of mass values for the spectrum
		"""

		value = array_as_numeric(value)

		# if not isinstance(value, _list_types) or not isinstance(value[0], Number):
		# 	raise TypeError("'mass_list' must be a list of numbers")

		# if not len(self.mass_list) == len(value):
		# 	raise ValueError("'mass_list' and 'intensity_list' are not the same size")

		self._mass_list = list(value)

		if self:
			self._min_mass = min(value)
			self._max_mass = max(value)
		else:
			self._min_mass = None
			self._max_mass = None

	def crop(
			self,
			min_mz: Optional[float] = None,
			max_mz: Optional[float] = None,
			inplace: bool = False,
			) -> "MassSpectrum":
		"""
		Crop the Mass Spectrum between the given mz values

		:param min_mz: The minimum mz for the new mass spectrum
		:param max_mz: The maximum mz for the new mass spectrum
		:param inplace: Whether the cropping should be applied this instance or to a copy (default behaviour).

		:return: The cropped Mass Spectrum
		"""

		if min_mz is None:
			min_mz = self.min_mass

		if max_mz is None:
			max_mz = self.max_mass

		min_mz_idx = self.intensity_list.index(min_mz)
		max_mz_idx = self.intensity_list.index(max_mz) + 1

		return self.icrop(min_mz_idx, max_mz_idx, inplace)

	def icrop(
			self,
			min_index: int = 0,
			max_index: int = -1,
			inplace: bool = False,
			) -> "MassSpectrum":
		"""
		Crop the Mass Spectrum between the given indices

		:param min_index: The minimum index for the new mass spectrum
		:param max_index: The maximum index for the new mass spectrum
		:param inplace: Whether the cropping should be applied this instance or to a copy (default behaviour).

		:return: The cropped Mass Spectrum
		"""

		cropped_intensity_list = self.intensity_list[min_index:max_index]
		cropped_mass_list = self.mass_list[min_index:max_index]

		if inplace:
			self.intensity_list = cropped_intensity_list
			self.mass_list = cropped_mass_list
			return self
		else:
			return self.__class__(
					intensity_list=cropped_intensity_list,
					mass_list=cropped_intensity_list,
					)

	def n_largest_peaks(self, n: int):
		"""
		Returns the indices of the n largest peaks in the Mass Spectrum

		:param n: The number of peaks to return the indices for
		"""

		# Make copies of the intensity_list
		intensity_list = self.intensity_list

		largest_indices = []

		for i in range(0, n):
			max_int_index = max(range(len(intensity_list)), key=intensity_list.__getitem__)

			del intensity_list[max_int_index]

			largest_indices.append(max_int_index)

		return largest_indices

	def get_intensity_for_mass(self, mass: float) -> float:
		"""
		Returns the intensity for the given mass.

		:param mass:
		"""

		mass_idx = self._mass_list.index(mass)
		return self._intensity_list[mass_idx]

	def get_mass_for_intensity(self, intensity: float) -> float:
		"""
		Returns the mass for the given intensity.
		If more than one mass has the given intensity, the first mass is returned.

		:param intensity:
		"""

		intensity_idx = self._intensity_list.index(intensity)
		return self._mass_list[intensity_idx]

	@classmethod
	def from_jcamp(cls, file_name: Union[str, pathlib.Path]) -> "MassSpectrum":
		"""
		Create a MassSpectrum from a JCAMP-DX file

		:param file_name: Path of the file to read

		:return: MassSpectrum

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic, David Kainer, Dominic Davis-Foster
		"""

		if not is_path(file_name):
			raise TypeError("'file_name' must be a string or a PathLike object")

		file_name = prepare_filepath(file_name, mkdirs=False)

		print(f" -> Reading JCAMP file '{file_name}'")
		lines_list = file_name.open('r')
		xydata = []
		last_tag = None

		for line in lines_list:

			if line.strip():
				if line.startswith("##"):
					# key word or information
					fields = line.split('=', 1)
					current_tag = fields[0] = fields[0].lstrip("##").upper()
					last_tag = fields[0]

					if current_tag.upper().startswith("END"):
						break

				else:
					if last_tag in xydata_tags:
						line_sub = re.split(r",| ", line.strip())
						for item in line_sub:
							if not len(item.strip()) == 0:
								xydata.append(float(item.strip()))

		# By this point we should have all of the xydata
		if len(xydata) % 2 == 1:
			# TODO: This means the data is not in x, y pairs
			#  Make a better error message
			raise ValueError("data not in pair !")

		mass_list = []
		intensity_list = []
		for i in range(len(xydata) // 2):
			mass_list.append(xydata[i * 2])
			intensity_list.append(xydata[i * 2 + 1])

		return cls(mass_list, intensity_list)

	@classmethod
	def from_mz_int_pairs(cls, mz_int_pairs: Sequence[Tuple[float, float]]):
		"""
		Construct a MassSpectrum from a list of (m/z, intensity) tuples.

		:param mz_int_pairs:
		"""

		err_msg = "`mz_int_pairs` must be a list of (m/z, intensity) tuples."

		if (
				not is_sequence(mz_int_pairs) or not is_sequence(mz_int_pairs[0])
				# or not isinstance(mz_int_pairs[0][0], Number)
				):
			raise TypeError(err_msg)

		if not len(mz_int_pairs[0]) == 2:
			raise ValueError(err_msg)

		mass_list = []
		intensity_list = []
		for mass, intensity in mz_int_pairs:
			mass_list.append(float(mass))
			intensity_list.append(float(intensity))

		return cls(mass_list, intensity_list)

	def __bool__(self) -> bool:
		if self.mass_list or self.mass_spec:
			return True
		return False


def normalize_mass_spec(
		mass_spec: MassSpectrum,
		relative_to: Optional[float] = None,
		inplace: bool = False,
		max_intensity: float = 100,
		) -> MassSpectrum:
	"""
	Normalize the intensities in the given Mass Spectrum to values between 0 and ``max_intensity``,
	which by default is 100.0.

	:param mass_spec: The Mass Spectrum to normalize
	:param relative_to: The largest intensity in the original data set.
		If not None the intensities are computed relative to this value.
		If None the value is calculated from the mass spectrum.
		This can be useful when normalizing several mass spectra to each other.
	:param inplace: Whether the normalization should be applied to the
		:class:`~pyms.Spectrum.MassSpectrum` object given, or to a copy (default behaviour).
	:param max_intensity: The maximum intensity in the normalized spectrum.
		If omitted the range 0-100.0 is used.
		If an integer the normalized intensities will be integers.

	:return: The normalized mass spectrum
	"""

	if relative_to is None:
		relative_to = max(mass_spec.intensity_list)

	normalized_intensity_list = [(x / float(relative_to)) * max_intensity for x in mass_spec.intensity_list]

	if isinstance(max_intensity, int):
		normalized_intensity_list = [round(x) for x in normalized_intensity_list]

	if inplace:
		mass_spec.intensity_list = normalized_intensity_list
		return mass_spec
	else:
		normalized_mass_spec = MassSpectrum(mass_spec.mass_list, normalized_intensity_list)

		return normalized_mass_spec

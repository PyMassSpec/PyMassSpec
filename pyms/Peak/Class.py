"""
Provides a class to model signal peak.
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
import warnings
from typing import Dict, List, Optional, Sequence, Tuple, Union, cast, overload
from warnings import warn

# this package
from pyms.Base import pymsBaseClass
from pyms.IntensityMatrix import BaseIntensityMatrix
from pyms.Spectrum import MassSpectrum
from pyms.Utils.Utils import is_number, is_sequence

__all__ = ["AbstractPeak", "Peak", "ICPeak"]


class AbstractPeak(pymsBaseClass):
	"""
	Models a signal peak.

	:param rt: Retention time.
	:param minutes: Retention time units flag. If :py:obj:`True`, retention time
		is in minutes; if :py:obj:`False` retention time is in seconds.
	:param outlier: Whether the peak is an outlier.

	:authors: Vladimir Likic, Andrew Isaac,
		Dominic Davis-Foster (type assertions and properties), David Kainer (outlier flag)

	.. versionadded:: 2.3.0
	"""

	def __init__(
			self,
			rt: Union[int, float] = 0.0,
			minutes: bool = False,
			outlier: bool = False,
			):

		if not is_number(rt):
			raise TypeError("'rt' must be a number")

		if minutes:
			rt = rt * 60.0

		# basic peak attributes
		self.is_outlier = outlier
		self._rt = float(rt)
		self._pt_bounds: Optional[Tuple[int, int, int]] = None
		self._area: Optional[float] = None
		self._ion_areas: Dict[float, float] = {}

		self.make_UID()

	def __eq__(self, other) -> bool:  # noqa: MAN001
		"""
		Return whether this Peak object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return (
					self.UID == other.UID and self.bounds == other.bounds and self.rt == other.rt
					and self.area == other.area
					)

		return NotImplemented

	# def __copy__(self):
	# 	#return pickle.loads(pickle.dumps(self))
	#
	# 	if self._mass_spectrum is None:
	# 		peak = Peak(rt=copy.copy(self._rt),
	# 					ms=copy.copy(self._ic_mass),
	# 					minutes=self._minutes,
	# 					outlier=self.is_outlier)
	# 	else:
	# 		peak = Peak(rt=copy.copy(self._rt),
	# 					ms=copy.copy(self._mass_spectrum),
	# 					minutes=self._minutes,
	# 					outlier=self.is_outlier)
	# 	if self._area is not None:
	# 		peak.area = self.area
	# 	if self._pt_bounds is not None:
	# 		peak.bounds = copy.copy(self.bounds)
	# 	if self._ic_mass is not None:
	# 		peak.ic_mass = 0+self.ic_mass
	#
	# 	return peak
	#
	# def __deepcopy__(self, memodict={}):
	# 	return self.__copy__()

	@property
	def area(self) -> Optional[float]:
		"""
		The area under the peak.

		:author: Andrew Isaac
		"""

		return self._area

	@area.setter
	def area(self, value: float) -> None:
		"""
		Sets the area under the peak.

		:param value: The peak area

		:author: Andrew Isaac
		"""

		if not is_number(value):
			raise TypeError("'Peak.area' must be a positive number")
		elif value < 0:
			raise ValueError("'Peak.area' must be a positive number")
		elif value <= 0:
			warnings.warn(f"Peak {self!r} @ {self.rt} seconds has zero area.")

		self._area = value

	@property
	def bounds(self) -> Optional[Tuple[int, int, int]]:
		"""
		The peak boundaries in points.

		:return: A 3-element tuple containing the left, apex, and right
			peak boundaries in points. Left and right are offsets.

		:author: Andrew Isaac
		"""

		return self._pt_bounds

	@bounds.setter
	def bounds(self, value: Sequence[int]) -> None:
		"""
		Sets peak boundaries in points.

		:param value: A 3-element tuple containing the left, apex, and right
			peak boundaries in points. Left and right are offsets.
		"""

		if not is_sequence(value):
			raise TypeError("'Peak.bounds' must be a Sequence")

		if len(value) != 3:
			raise ValueError("'Peak.bounds' must have exactly 3 elements")

		for index, item in enumerate(value):
			if not isinstance(item, int):
				raise TypeError(f"'Peak.bounds' element #{index} must be an integer")

		self._pt_bounds = cast(Tuple[int, int, int], tuple(value[:3]))

	def get_ion_area(self, ion: float) -> Optional[float]:
		"""
		Returns the area of a single ion chromatogram under the peak.

		:param ion: The ion to calculate the area for.

		:return: The area of the ion under this peak.
		"""

		try:
			return self._ion_areas[ion]
		except KeyError:
			return None

	@property
	def ion_areas(self) -> Dict:
		"""
		Returns a copy of the ion areas dict.

		:return: The dictionary of ``ion: ion area`` pairs
		"""  # noqa: D400

		if len(self._ion_areas) == 0:
			raise ValueError("no ion areas set")

		return copy.deepcopy(self._ion_areas)

	@ion_areas.setter
	def ion_areas(self, value: Dict) -> None:
		"""
		Sets the ``ion: ion area`` pairs dictionary.

		:param value: The dictionary of ion:ion_area pairs
		"""

		if not isinstance(value, dict) or not is_number(list(value.keys())[0]):
			raise TypeError("'Peak.ion_areas' must be a dictionary of ion:ion_area pairs")

		self._ion_areas = value

	def make_UID(self) -> None:
		"""
		Create a unique peak ID (UID).

		The UID comprises the retention time of the peak to two decimal places.
		Subclasses may define a more unique ID.

		:author: Andrew Isaac
		"""

		self._UID = f"{self._rt:.2f}"

	@property
	def rt(self) -> float:
		"""
		The retention time of the peak, in seconds.
		"""

		return self._rt

	def set_bounds(self, left: int, apex: int, right: int) -> None:
		"""
		Sets peak boundaries in points.

		:param left: Left peak boundary, in points offset from apex
		:param apex: Apex of the peak, in points
		:param right: Right peak boundary, in points offset from apex
		"""

		self.bounds = (left, apex, right)

	def set_ion_area(self, ion: int, area: float) -> None:
		"""
		Sets the area for a single ion.

		:param ion: the ion whose area is being entered.
		:param area: the area under the IC of ion.

		:author: Sean O'Callaghan
		"""

		if not isinstance(ion, int):
			raise TypeError("'ion' must be an integer")

		if not is_number(area):
			raise TypeError("'area' must be a number")

		self._ion_areas[ion] = area

	@property
	def UID(self) -> str:
		"""
		Return the unique peak ID (UID), either:

		- Integer masses of top two intensities and their ratio (as ``Mass1-Mass2-Ratio*100``); or
		- the single mass as an integer and the retention time.

		:return: UID string

		:author: Andrew Isaac
		"""  # noqa: D400

		return self._UID

	# def __dict__(self):
	#
	# 	return {
	# 			"UID": self.UID,
	# 			"bounds": self.bounds,
	# 			"area": self.area,
	# 			"is_outlier": self.is_outlier,
	# 			"ion_areas": self.ion_areas,
	# 			"mass_spectrum": self.mass_spectrum,
	# 			"rt": self.rt,
	# 			"ic_mass": self.ic_mass,
	#
	#
	#
	# 			}
	#
	# def __iter__(self):
	# 	for key, value in self.__dict__().items():
	# 		yield key, value


class Peak(AbstractPeak):
	"""
	Subclass of :class:`~.Peak` representing a peak in a mass spectrum.

	:param rt: Retention time.
	:param ms: The mass spectrum at the apex of the peak.
	:param minutes: Retention time units flag. If :py:obj:`True`, retention time
		is in minutes; if :py:obj:`False` retention time is in seconds.
	:param outlier: Whether the peak is an outlier.

	:authors: Vladimir Likic, Andrew Isaac,
		Dominic Davis-Foster (type assertions and properties), David Kainer (outlier flag)

	.. versionchanged:: 2.3.0

		Functionality related to single ion peaks has moved to the :class:`~.ICPeak` class.
		The two classes share a common base class, :class:`AbstractPeak`, which can be used
		in type checks for functions that accept either type of peak.

	.. versionchanged:: 2.3.0

		If the ``ms`` argument is unset an empty mass spectrum is used,
		rather than :py:obj:`None` in previous versions.

	.. TODO:: Change type hint of ``ms`` to Optional[MassSpectrum] once __new__ removed.
	"""

	_mass_spectrum: MassSpectrum

	def __init__(
			self,
			rt: Union[int, float] = 0.0,
			ms: Union[float, MassSpectrum, None] = None,
			minutes: bool = False,
			outlier: bool = False,
			):

		if ms is None:
			ms = MassSpectrum([], [])
		elif not isinstance(ms, MassSpectrum):
			raise TypeError("'ms' must be a MassSpectrum object")

		self._mass_spectrum = ms

		super().__init__(rt, minutes, outlier)

	@overload
	def __new__(
			cls,
			rt: float,
			ms: Optional[MassSpectrum],
			minutes: bool = ...,
			outlier: bool = ...,
			) -> "Peak": ...

	@overload
	def __new__(  # type: ignore[misc]
		cls,
		rt: float,
		ms: float,
		minutes: bool = ...,
		outlier: bool = ...,
		) -> "ICPeak": ...

	def __new__(
			cls,
			rt: float = 0.0,
			ms: Union[float, MassSpectrum, None] = None,
			minutes: bool = False,
			outlier: bool = False,
			):

		if is_number(ms):
			warnings.warn(
					"Creating a Peak object for a single ion chromatogram is "
					"deprecated in 2.3.0 and will be removed in 3.0.0. Use the ICPeak class instead."
					)
			return ICPeak(rt, cast(float, ms), minutes, outlier)
		else:
			obj = super().__new__(cls)
			obj.__init__(rt, ms, minutes, outlier)  # type: ignore[misc]  # unhappy about calling __init__
			return obj

	def __eq__(self, other) -> bool:  # noqa: MAN001
		"""
		Return whether this Peak object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return (
					self.UID == other.UID and self.bounds == other.bounds and self.rt == other.rt
					and self.mass_spectrum == other.mass_spectrum and self.area == other.area
					)

		return NotImplemented

	def crop_mass(self, mass_min: float, mass_max: float) -> None:
		"""
		Crops mass spectrum.

		:param mass_min: Minimum mass value.
		:param mass_max: Maximum mass value.

		:author: Andrew Isaac
		"""

		if not self._mass_spectrum:
			raise ValueError("Mass spectrum is unset.")

		if not is_number(mass_min) or not is_number(mass_max):
			raise TypeError("'mass_min' and 'mass_max' must be numbers")

		if mass_min >= mass_max:
			raise ValueError("'mass_min' must be less than 'mass_max'")

		mass_list = self._mass_spectrum.mass_list

		if mass_min < min(mass_list):
			raise ValueError(f"'mass_min' is less than the smallest mass: {min(mass_list)}")

		if mass_max > max(mass_list):
			raise ValueError(f"'mass_max' is greater than the largest mass: {max(mass_list)}")

		# pre build mass_list and list of indices
		new_mass_list = []
		new_mass_spec = []
		mass_spec = self._mass_spectrum.mass_spec
		for ii in range(len(mass_list)):
			mass = mass_list[ii]
			if mass_min <= mass <= mass_max:
				new_mass_list.append(mass)
				new_mass_spec.append(mass_spec[ii])

		self._mass_spectrum.mass_list = new_mass_list
		self._mass_spectrum.mass_spec = new_mass_spec

		if len(new_mass_list) == 0:
			raise ValueError("mass spectrum is now empty")
		elif len(new_mass_list) < 10:
			warn("peak mass spectrum contains < 10 points", Warning)

		# update UID
		self.make_UID()

	def get_int_of_ion(self, ion: int) -> float:
		"""
		Returns the intensity of a given ion in this peak.

		:param ion: The m/z value of the ion of interest
		"""

		try:
			index = self._mass_spectrum.mass_list.index(ion)
			return self._mass_spectrum.mass_spec[index]
		except (ValueError, IndexError):
			raise IndexError(
					f"'ion' out of range of mass spectrum (range "
					f"{min(self._mass_spectrum.mass_list)} to "
					f"{max(self._mass_spectrum.mass_list)})"
					)

	def get_third_highest_mz(self) -> int:
		"""
		Returns the *m/z* value with the third highest intensity.
		"""

		if not self._mass_spectrum:
			raise ValueError("Mass spectrum is unset.")

		mass_list = self._mass_spectrum.mass_list
		mass_spec = self._mass_spectrum.mass_spec
		# find top two masses
		best = 0
		best_ii = 0
		best2_ii = 0
		best3_ii = 0
		for ii, intensity in enumerate(mass_spec):
			if intensity > best:
				best = intensity
				best3_ii = best2_ii
				best2_ii = best_ii
				best_ii = ii

		return int(mass_list[best3_ii])

	def make_UID(self) -> None:
		"""
		Create a unique peak ID (UID):

		- Integer masses of top two intensities and their ratio (as ``Mass1-Mass2-Ratio*100``); or

		:author: Andrew Isaac
		"""  # noqa: D400

		if self._mass_spectrum:
			mass_list = self._mass_spectrum.mass_list
			mass_spec = self._mass_spectrum.mass_spec
			# find top two masses
			best = 0
			best_ii = 0
			best2_ii = 0
			for ii in range(len(mass_spec)):
				if mass_spec[ii] > best:
					best = mass_spec[ii]
					best2_ii = best_ii
					best_ii = ii
			ratio = int(100 * mass_spec[best2_ii] / best)
			self._UID = f"{int(mass_list[best_ii]):d}-{int(mass_list[best2_ii]):d}-{ratio:d}-{self._rt:.2f}"
		else:
			super().make_UID()

	@property
	def mass_spectrum(self) -> MassSpectrum:
		"""
		The mass spectrum at the apex of the peak.
		"""

		return copy.copy(self._mass_spectrum)

	@mass_spectrum.setter
	def mass_spectrum(self, value: MassSpectrum) -> None:
		"""
		Sets the mass spectrum for the apex of the peak.
		"""

		if not isinstance(value, MassSpectrum):
			raise TypeError("'Peak.mass_spectrum' must be a MassSpectrum object")

		self._mass_spectrum = value
		self.make_UID()

	def null_mass(self, mass: float) -> None:
		"""
		Ignore given mass in spectra.

		:param mass: Mass value to remove

		:author: Andrew Isaac
		"""

		if not self._mass_spectrum:
			raise ValueError("Mass spectrum is unset.")

		if not is_number(mass):
			raise TypeError("'mass' must be a number")

		mass_list = self._mass_spectrum.mass_list

		if mass < min(mass_list) or mass > max(mass_list):
			raise IndexError("'mass' not in mass range:", min(mass_list), "to", max(mass_list))

		best = max(mass_list)
		ix = 0
		for ii in range(len(mass_list)):
			tmp = abs(mass_list[ii] - mass)
			if tmp < best:
				best = tmp
				ix = ii

		self._mass_spectrum.mass_spec[ix] = 0

		# update UID
		self.make_UID()

	def find_mass_spectrum(self, data: BaseIntensityMatrix, from_bounds: float = False) -> None:
		"""
		.. TODO:: What does this function do?

		Sets the peak's mass spectrum from the data.

		Clears the single ion chromatogram mass.

		:param data:
		:param from_bounds: Whether to use the attribute :attr:`pyms.Peak.Class.Peak.pt_bounds`
			or to find the peak apex from the peak retention time.
		"""  # noqa: D400

		if not isinstance(data, BaseIntensityMatrix):
			raise TypeError("'data' must be an IntensityMatrix")

		if from_bounds:
			if self._pt_bounds is None:
				raise NameError("pt_bounds not set for this peak")
			else:
				pt_apex = self._pt_bounds[1]
		else:
			# get the index of peak apex from peak retention time
			pt_apex = data.get_index_at_time(self._rt)

		# set the mass spectrum
		self._mass_spectrum = data.get_ms_at_index(pt_apex)

		# TODO: something about this function for ICPeak
		# clear single ion chromatogram mass
		# self._ic_mass = None
		self.make_UID()

	def top_ions(self, num_ions: int = 5) -> List[float]:
		"""
		Computes the highest #num_ions intensity ions.

		:param num_ions: The number of ions to be recorded.

		:return: A list of the ions with the highest intensity.

		:authors: Sean O'Callaghan, Dominic Davis-Foster (type assertions)
		"""

		if not self._mass_spectrum:
			raise ValueError("Mass spectrum is unset.")

		if not isinstance(num_ions, int):
			raise TypeError("'n_top_ions' must be an integer")

		intensity_list = self.mass_spectrum.mass_spec
		mass_list = self.mass_spectrum.mass_list

		ic_tuple = zip(intensity_list, mass_list)

		sorted_ic = sorted(ic_tuple)
		top_ic = sorted_ic[-num_ions:]

		top_ions = []

		for entry in top_ic:
			top_ions.append(entry[1])

		return top_ions


class ICPeak(AbstractPeak):
	"""
	Subclass of :class:`~.Peak` representing a peak in an ion chromatogram for a single mass.

	:param rt: Retention time.
	:param mass: The mass of the ion.
	:param minutes: Retention time units flag. If :py:obj:`True`, retention time
		is in minutes; if :py:obj:`False` retention time is in seconds.
	:param outlier: Whether the peak is an outlier.

	:authors: Vladimir Likic, Andrew Isaac,
		Dominic Davis-Foster (type assertions and properties), David Kainer (outlier flag)

	.. versionadded:: 2.3.0
	"""

	_ic_mass: Optional[float]

	def __init__(
			self,
			rt: Union[int, float] = 0.0,
			mass: Optional[float] = None,
			minutes: bool = False,
			outlier: bool = False,
			):

		if mass and not is_number(mass):
			raise TypeError("'ms' must be a number")

		self._ic_mass = mass

		super().__init__(rt, minutes, outlier)

	def __eq__(self, other) -> bool:  # noqa: MAN001
		"""
		Return whether this Peak object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return (
					self.UID == other.UID and self.bounds == other.bounds and self.rt == other.rt
					and self.area == other.area and self.ic_mass == other.ic_mass
					)

		return NotImplemented

	@property
	def ic_mass(self) -> Optional[float]:
		"""
		The mass for a single ion chromatogram peak.

		:return: The mass of the single ion chromatogram that the peak is from
		"""

		return self._ic_mass

	@ic_mass.setter
	def ic_mass(self, value: float) -> None:
		"""
		Sets the mass for a single ion chromatogram peak and clears the mass spectrum.

		:param value: The mass of the ion chromatogram that the peak is from
		"""

		if not is_number(value):
			raise TypeError("'Peak.ic_mass' must be a number")

		self._ic_mass = value
		self.make_UID()

	def make_UID(self) -> None:
		"""
		Create a unique peak ID (UID):

		- the single mass as an integer and the retention time.

		:author: Andrew Isaac
		"""  # noqa: D400

		if self._ic_mass is not None:
			self._UID = f"{int(self._ic_mass):d}-{self._rt:.2f}"
		else:
			super().make_UID()

"""
Provides a class to model signal peak
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
from numbers import Number
from warnings import warn

# 3rd party
import deprecation

# this package
from pyms import __version__
from pyms.Base import pymsBaseClass
from pyms.IntensityMatrix import IntensityMatrix
from pyms.Spectrum import MassSpectrum
from pyms.Utils.Utils import is_sequence


class Peak(pymsBaseClass):
	"""
	Models a signal peak

	:param rt: Retention time
	:type rt: int or float
	:param ms: Either an ion mass, a mass spectrum or None
	:type ms: int or float or class`pyms.Spectrum.MassSpectrum`, optional
	:param minutes: Retention time units flag. If True, retention time
		is in minutes; if False retention time is in seconds
	:type minutes: bool, optional
	:param outlier: Whether the peak is an outlier
	:type outlier: bool, optional

	:author: Vladimir Likic
	:author: Andrew Isaac
	:author: Dominic Davis-Foster (type assertions and properties)
	:author: David Kainer (outlier flag)
	"""

	def __init__(self, rt=0.0, ms=None, minutes=False, outlier=False):
		"""
		Models a signal peak
		"""

		if not isinstance(rt, Number):
			raise TypeError("'rt' must be a number")

		if ms and not isinstance(ms, MassSpectrum) and not isinstance(ms, Number):
			raise TypeError("'ms' must be a number or a MassSpectrum object")

		if minutes:
			rt = rt * 60.0

		# basic peak attributes
		self.is_outlier = outlier
		self._rt = float(rt)
		self._pt_bounds = None
		self._area = None
		self._ion_areas = {}

		# these two attributes are required for
		# setting the peak mass spectrum
		if isinstance(ms, MassSpectrum):
			# mass spectrum
			self._mass_spectrum = ms
			self._ic_mass = None
		else:
			# single ion chromatogram properties
			self._mass_spectrum = None
			self._ic_mass = ms

		self.make_UID()

	def __eq__(self, other):
		"""
		Return whether this Peak object is equal to another object

		:param other: The other object to test equality with
		:type other: object

		:rtype: bool
		"""

		if isinstance(other, self.__class__):
			return self.UID == other.UID \
					and self.bounds == other.bounds \
					and self.rt == other.rt \
					and self.mass_spectrum == other.mass_spectrum \
					and self.area == other.area

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
	def area(self):
		"""
		Gets the area under the peak

		:return: The peak area
		:rtype: float

		:author: Andrew Isaac
		"""

		return self._area

	@area.setter
	def area(self, value):
		"""
		Sets the area under the peak

		:param value: The peak area
		:type value: float

		:author: Andrew Isaac
		"""

		if not isinstance(value, Number):
			raise TypeError("'Peak.area' must be a positive number")
		elif value <= 0:
			raise ValueError("'Peak.area' must be a positive number")

		self._area = value

	@property
	def bounds(self):
		"""
		Gets peak boundaries in points

		:return: A list containing left, apex, and right
			peak boundaries in points, left and right are offsets
		:rtype: list

		:author: Andrew Isaac
		"""

		return self._pt_bounds

	@bounds.setter
	def bounds(self, value):
		"""
		Sets peak boundaries in points

		:param value: A list containing left, apex, and right
			peak boundaries in points, left and right are offsets
		:type value: ~collections.abc.Sequence
		"""

		if not is_sequence(value):
			raise TypeError("'Peak.bounds' must be a Sequence")

		if len(value) != 3:
			raise ValueError("'Peak.bounds' must have exactly 3 elements")

		for index, item in enumerate(value):
			if not isinstance(item, int):
				raise TypeError(f"'Peak.bounds' element #{index} must be an integer")

		self._pt_bounds = value

	def crop_mass(self, mass_min, mass_max):
		"""
		Crops mass spectrum

		:param mass_min: Minimum mass value
		:type mass_min: Number
		:param mass_max: Maximum mass value
		:type mass_max: Number

		:author: Andrew Isaac
		"""

		if not isinstance(mass_min, Number) or not isinstance(mass_max, Number):
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

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.area` instead")
	def get_area(self):
		"""
		Gets the area under the peak

		:return: The peak area
		:rtype: float

		:author: Andrew Isaac
		"""

		return self.area

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.ic_mass` instead")
	def get_ic_mass(self):
		"""
		Returns the mass for a single ion chromatogram peak

		:rtype: float or int
		"""

		return self.ic_mass

	def get_int_of_ion(self, ion):
		"""
		Returns the intensity of a given ion

		:param ion: The m/z value of the ion of interest
		:type ion: int

		:return: The intensity of the given ion in this peak
		:rtype: int
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

	def get_ion_area(self, ion):
		"""
		Returns the area of a single ion chromatogram under the peak

		:param ion: The ion to calculate the area for
		:type ion: int

		:return: The area of the ion under this peak
		:rtype: float
		"""

		try:
			return self._ion_areas[ion]
		except KeyError:
			return None

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.ion_areas` instead")
	def get_ion_areas(self):
		"""
		Returns a copy of the ion areas dict containing ion:ion area pairs

		:rtype: dict
		"""

		return self.ion_areas

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.mass_spectrum` instead")
	def get_mass_spectrum(self):
		"""
		Returns the mass spectrum at the apex of the peak

		:rtype: pyms.Spectrum.MassSpectrum
		"""

		return self.mass_spectrum

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.bounds` instead")
	def get_pt_bounds(self):
		"""
		Returns the peak boundaries in points

		:return: A list containing left, apex, and right peak boundaries in
			points, left and right are offsets
		:rtype: list

		:author: Andrew Isaac
		"""

		return self.bounds

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.rt` instead")
	def get_rt(self):
		"""
		Returns the retention time

		:rtype: float
		"""

		return self.rt

	def get_third_highest_mz(self):
		"""
		Returns the mz value with the third highest intensity

		:rtype: int
		"""

		if self._mass_spectrum is not None:
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

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.UID` instead")
	def get_UID(self):
		"""
		Return the unique peak ID (UID), either:

		- Integer masses of top two intensities and their ratio (as Mass1-Mass2-Ratio*100); or
		- the single mass as an integer and the retention time.

		:return: UID string
		:rtype: str

		:author: Andrew Isaac
		"""

		return self.UID

	@property
	def ic_mass(self):
		"""
		Gets the mass for a single ion chromatogram peak

		:return: The mass of the single ion chromatogram that the peak is from
		:rtype: float or int
		"""

		return self._ic_mass

	@ic_mass.setter
	def ic_mass(self, value):
		"""
		Sets the mass for a single ion chromatogram peak
			Clears the mass spectrum

		:param value: The mass of the ion chromatogram that the peak is from
		:type value: float
		"""

		if not isinstance(value, Number):
			raise TypeError("'Peak.ic_mass' must be a number")
		self._ic_mass = value
		# clear mass spectrum
		self._mass_spectrum = None
		self.make_UID()

	@property
	def ion_areas(self):
		"""
		Returns a copy of the ion areas dict

		:return: The dictionary of ion:ion area pairs
		:rtype: dict
		"""
		if len(self._ion_areas) == 0:
			raise ValueError("no ion areas set")

		return copy.deepcopy(self._ion_areas)

	@ion_areas.setter
	def ion_areas(self, value):
		"""
		Sets the ion:ion area pair dictionary

		:param value: The dictionary of ion:ion_area pairs
		:type value: dict
		"""

		if not isinstance(value, dict) or not isinstance(list(value.keys())[0], Number):
			raise TypeError("'Peak.ion_areas' must be a dictionary of ion:ion_area pairs")
		self._ion_areas = value

	def make_UID(self):
		"""
		Create a unique peak ID (UID), either:

		- Integer masses of top two intensities and their ratio (as Mass1-Mass2-Ratio*100); or
		- the single mass as an integer and the retention time.

		:author: Andrew Isaac
		"""

		if self._mass_spectrum is not None:
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
		elif self._ic_mass is not None:
			self._UID = f"{int(self._ic_mass):d}-{self._rt:.2f}"
		else:
			self._UID = f"{self._rt:.2f}"

	@property
	def mass_spectrum(self):
		"""
		Gets the mass spectrum at the apex of the peak

		:return: The mass spectrum at the apex of the peak
		:rtype: pyms.Spectrum.MassSpectrum
		"""

		return copy.copy(self._mass_spectrum)

	@mass_spectrum.setter
	def mass_spectrum(self, value):
		"""
		Sets the mass spectrum
			Clears the mass for a single ion chromatogram peak

		:param value: The mass spectrum at the apex of the peak
		:rtype: pyms.Spectrum.MassSpectrum
		"""

		if not isinstance(value, MassSpectrum):
			raise TypeError("'Peak.mass_spectrum' must be a MassSpectrum object")

		self._mass_spectrum = value
		# clear ion mass
		self._ic_mass = None
		self.make_UID()

	def null_mass(self, mass):
		"""
		Ignore given mass in spectra

		:param mass: Mass value to remove
		:type mass: int or float

		:author: Andrew Isaac
		"""

		if self._mass_spectrum is None:
			raise NameError("mass spectrum not set for this peak")

		if not isinstance(mass, Number):
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

	@property
	def rt(self):
		"""
		Return the retention time

		:return: Retention time
		:rtype: float
		"""

		return self._rt

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.area` instead")
	def set_area(self, area):
		"""
		Sets the area under the peak

		:param area: The peak area
		:type area: float

		:author: Andrew Isaac
		"""

		if not isinstance(area, Number):
			raise TypeError("'area' must be a positive number")
		elif area <= 0:
			raise ValueError("'area' must be a positive number")

		self._area = area

	def set_bounds(self, left, apex, right):
		"""
		Sets peak boundaries in points

		:param left: Left peak boundary, in points offset from apex
		:param apex: Apex of the peak, in points
		:param right: Right peak boundary, in points offset from apex
		"""

		self.bounds = (left, apex, right)

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.ic_mass` instead")
	def set_ic_mass(self, mz):
		"""
		Sets the mass for a single ion chromatogram peak.
		Clears the mass spectrum.

		:param mz: The mass of the ion chromatogram that the peak is from
		:type mz: float
		"""

		if not isinstance(mz, (float, int)):
			raise TypeError("'mz' must be a number")
		self._ic_mass = mz
		# clear mass spectrum
		self._mass_spectrum = None
		self.make_UID()

	def set_ion_area(self, ion, area):
		"""
		sets the area for a single ion

		:param ion: the ion whose area is being entered
		:type ion: int
		:param area: the area under the IC of ion
		:type area: int or float

		:author: Sean O'Callaghan
		"""

		if not isinstance(ion, int):
			raise TypeError("'ion' must be an integer")

		if not isinstance(area, Number):
			raise TypeError("'area' must be a number")

		self._ion_areas[ion] = area

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.ion_areas` instead")
	def set_ion_areas(self, ion_areas):
		"""
		Set the ion:ion area pair dictionary

		:param ion_areas: The dictionary of ion:ion_area pairs
		:type ion_areas: dict
		"""

		self.ion_areas = ion_areas

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.mass_spectrum` instead")
	def set_mass_spectrum(self, ms):
		"""
		Sets the mass spectrum.
		Clears the mass for a single ion chromatogram peak.

		:param ms: The mass spectrum at the apex of the peak
		:type ms: pyms.Spectrum.MassSpectrum
		"""

		if not isinstance(ms, MassSpectrum):
			raise TypeError("'ms' must be a MassSpectrum object")

		self._mass_spectrum = ms
		# clear ion mass
		self._ic_mass = None
		self.make_UID()

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Peak.Class.Peak.bounds` instead")
	def set_pt_bounds(self, pt_bounds):
		"""
		Sets peak boundaries in points

		:param pt_bounds: A list containing left, apex, and right
			peak boundaries in points, left and right are offsets
		:type pt_bounds: list
		"""

		self.bounds = pt_bounds

	@property
	def UID(self):
		"""
		Return the unique peak ID (UID), either:

		- Integer masses of top two intensities and their ratio (as Mass1-Mass2-Ratio*100); or
		- the single mass as an integer and the retention time.

		:return: UID string
		:rtype: str

		:author: Andrew Isaac
		"""

		return self._UID

	def find_mass_spectrum(self, data, from_bounds=False):
		"""
		TODO: What does this function do?

		Sets peak mass spectrum from the data.
		Clears the single ion chromatogram mass.

		:param data: An IntensityMatrix object
		:type data: pyms.IntensityMatrix.IntensityMatrix
		:param from_bounds: Indicator whether to use the attribute 'pt_bounds'
			or to find the peak apex from the peak retention time
		:type from_bounds: bool
		"""

		if not isinstance(data, IntensityMatrix):
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

		# clear single ion chromatogram mass
		self._ic_mass = None
		self.make_UID()

	def top_ions(self, num_ions=5):
		"""
		Computes the highest #num_ions intensity ions

		:param num_ions: The number of ions to be recorded
		:type num_ions: int

		:return: A list of the #num_ions highest intensity ions, default 5
		:rtype: list, optional

		:author: Sean O'Callaghan
		:author: Dominic Davis-Foster (type assertions)
		"""

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
	#
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

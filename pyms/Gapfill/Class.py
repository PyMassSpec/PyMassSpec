"""
Provides a class for handling Missing Peaks in an output file (i.e. area.csv)
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

# 3rd party
import deprecation

# this package
from pyms import __version__


class MissingPeak:
	"""
	Class to encapsulate a peak object identified as missing in the output area
	matrix fom PyMassSpec.

	:param common_ion: Common ion for the peak across samples in an experiment
	:type common_ion: int
	:param qual_ion_1:
	:type qual_ion_1:
	:param qual_ion_2:
	:type qual_ion_2:
	:param rt: Retention time of the peak. Default ``0.0``
	:type rt: float, optional

	:author: Jairus Bowne
	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (properties)
	"""

	def __init__(self, common_ion: int, qual_ion_1, qual_ion_2, rt=0.0) -> float:
		"""
		Initialise MissingPeak Class
		"""

		self.__common_ion = common_ion
		self.__qual_1 = qual_ion_1
		self.__qual_2 = qual_ion_2
		self.__rt = rt
		self.__exact_rt = 'na'
		self.__common_ion_area = 'na'

	@property
	def common_ion(self) -> int:
		"""
		Returns the common ion for the peak object across an experiment

		:return: Common ion for the peak
		:rtype: int

		:author: Jairus Bowne
		"""

		return self.__common_ion

	@property
	def common_ion_area(self) -> int:
		"""
		Returns the common ion area

		:return: The area of the common ion
		:rtype: int
		"""

		return self.__common_ion_area

	@common_ion_area.setter
	def common_ion_area(self, common_ion_area: int):
		"""
		sets the common ion area calculated by the gap fill algorithm

		:param common_ion_area: The area of the common ion
		:type common_ion_area: int
		"""

		self.__common_ion_area = common_ion_area

	@property
	def exact_rt(self) -> float:
		"""
		Returns the retention time of the peak

		:return: the retention time of the peak
		:rtype: float
		"""

		return self.__exact_rt

	@exact_rt.setter
	def exact_rt(self, rt) -> float:
		"""
		Sets the retention time of a peak

		:param rt: The retention time of the apex of the peak
		:type rt: float
		"""

		self.__exact_rt = rt

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.common_ion` instead")
	def get_common_ion(self) -> int:
		"""
		Returns the common ion for the peak object across an experiment

		:return: Common ion for the peak
		:rtype: int

		:author: Jairus Bowne
		"""

		return self.common_ion

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.get_common_ion_area` instead")
	def get_common_ion_area(self) -> int:
		"""
		returns the common ion area

		:return common_ion_area: The area of the common ion
		:rtype: int
		"""
		return self.common_ion_area

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.exact_rt` instead")
	def get_exact_rt(self) -> float:
		"""
		returns the retention time of the peak

		:return: the retention time of the peak
		:rtype: float
		"""

		return self.exact_rt

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.qual_ion1` instead")
	def get_qual_ion1(self) Int:
		"""
		Returns the top (most abundant) ion for the peak object

		:return: Most abundant ion
		:rtype: int

		:author: Jairus Bowne
		"""

		return self.qual_ion1

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.qual_ion2` instead")
	def get_qual_ion2(self) -> int
		"""
		Returns the second most abundant ion for the peak object

		:return: Second most abundant ion
		:rtype: int

		:author: Jairus Bowne
		"""

		return self.qual_ion2

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.rt` instead")
	def get_rt(self) -> float:
		"""
		returns the retention time of the peak

		:return: the retention time of the peak
		:rtype: float
		"""
		return self.rt

	@property
	def qual_ion1(self) -> int:
		"""
		Returns the top (most abundant) ion for the peak object

		:return: Most abundant ion
		:rtype: int

		:author: Jairus Bowne
		"""

		# TODO: Consider the abundance of ions when some (i.e. 73, 147) have
		#  been im.null_mass()'d. Is there a way to determine whether that
		#  has been done to generate the original peak list?

		return self.__qual_1

		# return int(string.split(self.__UID, '-')[0])

	@property
	def qual_ion2(self) -> int:
		"""
		Returns the second most abundant ion for the peak object

		:return: Second most abundant ion
		:rtype: int

		:author: Jairus Bowne
		"""

		# TODO: Consider the abundance of ions when some (i.e. 73, 147) have
		# 	been im.null_mass()'d. Is there a way to determine whether that
		# 	has been done to generate the original peak list?

		return self.__qual_1

		# return int(string.split(self.__UID, '-')[0])

	@property
	def rt(self) -> float:
		"""
		returns the retention time of the peak

		:return: the retention time of the peak
		:rtype: float
		"""

		return self.__rt

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.common_ion_area` instead")
	def set_common_ion_area(self, common_ion_area: int):
		"""
		sets the common ion area calculated by the gap fill algorithm.

		:param common_ion_area: The area of the common ion
		:type common_ion_area: int
		"""

		self.common_ion_area = common_ion_area

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.MissingPeak.exact_rt` instead")
	def set_exact_rt(self, rt) -> float:
		"""
		sets the retention time of a peak

		:param rt: The retention time of the apex of the peak
		:type rt: float
		"""

		self.__exact_rt = rt


class Sample:
	"""
	A collection of MissingPeak objects

	:param sample_name: the experiment code/name
	:type sample_name: str

	:param matrix_position: position along x-axis where sample is located
	:type matrix_position: int

	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (properties)
	"""

	def __init__(self, sample_name, matrix_position: int):
		"""
		A collection of MissingPeak objects
		"""

		self.__sample_name = sample_name
		self.__matrix_position = matrix_position
		self.__missing_peak_list = []

	def add_missing_peak(self, missing_peak):
		"""
		Add a new MissingPeak object to the Sample

		:param missing_peak: The missing peak object to be added
		:type missing_peak: pyms.Gapfill.Class.MissingPeak
		"""

		# TODO: Do some checking here!!!

		self.__missing_peak_list.append(missing_peak)

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.Sample.missing_peaks` instead")
	def get_missing_peaks(self) -> MissingPeak:
		"""
		Returns a list of the MissingPeak objects in the Sample object

		:return: list of the missing peaks
		:rtype: :class:`list` of :class:`pyms.Gapfill.Class.MissingPeak` objects
		"""
		return self.__missing_peak_list

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.Sample.rt_areas` instead")
	def get_mp_rt_area_dict(self) -> dict:
		"""
		Returns a dictionary containing rt:area pairs

		:return: a dict containing rt:area pairs
		:rtype: dict
		"""

		rt_area_dict = {}
		for peak in self.__missing_peak_list:
			rt = peak.rt
			area = peak.get_common_ion_area()

			rt_area_dict[rt] = area

		return rt_area_dict

	def get_mp_rt_exact_rt_dict(self) -> dict:
		"""
		Returns a dictionary containing average_rt:exact_rt pairs

		:return: a dict of average_rt:exact_rt pairs
		:rtype: dict
		"""

		rt_exact_rt_dict = {}
		for peak in self.__missing_peak_list:
			rt = peak.rt
			exact_rt = peak.get_exact_rt()

			rt_exact_rt_dict[rt] = exact_rt

		return rt_exact_rt_dict

	@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
							current_version=__version__,
							details="Use :attr:`pyms.Gapfill.Class.Sample.name` instead")
	def get_name(self) -> str:
		"""
		Returns the sample name

		:return: The name of the sample
		:rtype: str
		"""

		return self.__sample_name

	@property
	def missing_peaks(self) -> MissingPeak:
		"""
		Returns a list of the MissingPeak objects in the Sample object

		:return: list of the missing peaks
		:rtype: :class:`list` of :class:`pyms.Gapfill.Class.MissingPeak` objects
		"""
		return self.__missing_peak_list

	@property
	def name(self) -> str:
		"""
		Returns the sample name

		:return: The name of the sample
		:rtype: str
		"""

		return self.__sample_name

	@property
	def rt_areas(self) -> dict:
		"""
		returns a dictionary containing rt:area pairs

		:return: a dict containing rt:area pairs
		:rtype: dict
		"""

		rt_area_dict = {}
		for peak in self.__missing_peak_list:
			rt = peak.rt
			area = peak.get_common_ion_area()

			rt_area_dict[rt] = area

		return rt_area_dict

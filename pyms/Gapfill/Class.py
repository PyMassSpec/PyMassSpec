"""
Provides a class for handling Missing Peaks in an output file (i.e. area.csv)
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
from typing import Dict, List, Optional

__all__ = ["MissingPeak", "Sample"]


class MissingPeak:
	"""
	Class to encapsulate a peak object identified as missing in the output area
	matrix fom PyMassSpec.

	:param common_ion: Common ion for the peak across samples in an experiment.
	:param qual_ion_1: The top (most abundant) ion for the peak object
	:param qual_ion_2: The second most abundant ion for the peak object
	:param rt: Retention time of the peak.

	:authors: Jairus Bowne, Sean O'Callaghan, Dominic Davis-Foster
	"""

	#: The area of the common ion
	common_ion_area: Optional[float]

	#: The retention time of the apex of the peak
	exact_rt: Optional[float]

	def __init__(self, common_ion: int, qual_ion_1: int, qual_ion_2: int, rt: float = 0.0):
		self.__common_ion = common_ion
		self.__qual_1 = qual_ion_1
		self.__qual_2 = qual_ion_2
		self.__rt = rt
		self.exact_rt = None
		self.common_ion_area = None

	@property
	def common_ion(self) -> int:
		"""
		Returns the common ion for the peak object across an experiment

		:return: Common ion for the peak

		:author: Jairus Bowne
		"""

		return self.__common_ion

	@property
	def qual_ion1(self) -> int:
		"""
		Returns the top (most abundant) ion for the peak object

		:return: Most abundant ion

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
		"""

		return self.__rt


class Sample:
	"""
	A collection of MissingPeak objects.

	:param sample_name: the experiment code/name.
	:param matrix_position: position along x-axis where sample is located.

	:authors: Sean O'Callaghan, Dominic Davis-Foster (properties)
	"""

	def __init__(self, sample_name: str, matrix_position: int):
		self._sample_name = sample_name
		self._matrix_position = matrix_position
		self._missing_peak_list: List[MissingPeak] = []

	def add_missing_peak(self, missing_peak: MissingPeak) -> None:
		"""
		Add a new MissingPeak object to the Sample

		:param missing_peak: The missing peak object to be added
		"""

		# TODO: Do some checking here

		self._missing_peak_list.append(missing_peak)

	def get_mp_rt_exact_rt_dict(self) -> Dict[float, Optional[float]]:
		"""
		Returns a dictionary containing ``average_rt : exact_rt`` pairs.
		"""

		rt_exact_rt_dict = {}
		for peak in self._missing_peak_list:
			rt = peak.rt
			exact_rt = peak.exact_rt

			rt_exact_rt_dict[rt] = exact_rt

		return rt_exact_rt_dict

	@property
	def missing_peaks(self) -> List[MissingPeak]:
		"""
		Returns a list of the MissingPeak objects in the Sample object
		"""

		return self._missing_peak_list

	@property
	def name(self) -> str:
		"""
		Returns name of the sample.
		"""

		return self._sample_name

	@property
	def rt_areas(self) -> Dict[float, Optional[float]]:
		"""
		returns a dictionary containing ``rt : area`` pairs.
		"""

		rt_area_dict = {}
		for peak in self._missing_peak_list:
			rt = peak.rt
			area = peak.common_ion_area

			rt_area_dict[rt] = area

		return rt_area_dict

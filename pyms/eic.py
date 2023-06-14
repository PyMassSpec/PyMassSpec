"""
Class to model a subset of data from an Intensity Matrix.
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2020 Dominic Davis-Foster                                   #
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
from typing import Iterable, List, Optional, Sequence, Union, cast

# 3rd party
import numpy

# this package
from pyms.IntensityMatrix import BaseIntensityMatrix, IntensityMatrix
from pyms.IonChromatogram import BasePeakChromatogram, ExtractedIonChromatogram, IonChromatogram
from pyms.Utils.Utils import is_number

__all__ = ["ExtractedIntensityMatrix", "build_extracted_intensity_matrix"]


class ExtractedIntensityMatrix(BaseIntensityMatrix):
	"""
	Represents an extracted subset of the chromatographic data.

	:param time_list: Retention time values
	:param mass_list: Binned mass values
	:param intensity_array: List of lists of binned intensity values per scan

	:authors: Dominic Davis-Foster

	.. versionadded:: 2.3.0
	"""

	def get_ic_at_mass(self, mass: Optional[float] = None) -> IonChromatogram:
		"""
		Returns the ion chromatogram for the nearest binned mass to the specified mass.

		If no mass value is given, the function returns the extracted ion chromatogram.

		:param mass: Mass value of an ion chromatogram

		:return: Ion chromatogram for given mass

		:authors: Andrew Isaac, Vladimir Likic
		"""

		if mass is None:
			return self.eic
		elif not is_number(mass):
			raise TypeError("'mass' must be a number")

		if mass < self._min_mass or mass > self._max_mass:
			print("min mass: ", self._min_mass, "max mass:", self._max_mass)
			raise IndexError("mass is out of range")

		ix = self.get_index_of_mass(mass)

		return self.get_ic_at_index(ix)

	@property
	def eic(self) -> ExtractedIonChromatogram:
		"""
		Returns an :class:`~.IonChromatogram` object representing this EIC.
		"""

		intensity_list = []
		for row in self._intensity_array:
			intensity_list.append(sum(row))

		return ExtractedIonChromatogram(numpy.array(intensity_list), self._time_list[:], self.mass_list)

	@property
	def bpc(self) -> IonChromatogram:
		"""
		Constructs a Base Peak Chromatogram from the data.

		This represents the most intense ion
		-- out of those used to create the :class:`~.ExtractedIntensityMatrix`
		-- for each scan.

		:authors: Dominic Davis-Foster

		.. versionadded:: 2.3.0
		"""

		return BasePeakChromatogram(
				[max(intensities) for intensities in self._intensity_array],
				self._time_list[:],
				)


def build_extracted_intensity_matrix(
		im: IntensityMatrix,
		masses: Sequence[Union[float, Iterable[float]]],
		left_bound: float = 0.5,
		right_bound: float = 0.5,
		) -> ExtractedIntensityMatrix:
	"""
	Given an intensity matrix and a list of masses, construct a
	:class:`~.ExtractedIntensityMatrix` for those masses.

	The masses can either be:

	* single masses (of type :class:`float`),
	* an iterable of masses.

	``left_bound`` and ``right_bound`` specify a range in which to include values for around each mass.
	For example, a mass of ``169`` with bounds of ``0.3`` and ``0.7`` would include every mass
	between ``168.7`` and ``169.7`` (inclusive on both sides).

	Set the bounds to ``0`` to include only the given masses.

	:param im:
	:param masses:
	:param left_bound:
	:param right_bound:
	"""  # noqa: D400

	flat_target_masses: List[float] = []

	for mass in masses:
		if is_number(mass):
			mass = cast(float, mass)
			flat_target_masses.append(mass)
		elif isinstance(masses, (Sequence, Iterable, range)):
			mass = cast(Union[Sequence, Iterable, range], mass)
			flat_target_masses.extend(mass)
		else:
			raise NotImplementedError(f"Unsupported type '{type(mass)}'")

	# get indices of all those masses, taking bounds into account.
	indices = set()

	for idx, mass in enumerate(im.mass_list):
		for target_mass in flat_target_masses:
			if (target_mass - left_bound) <= mass <= (target_mass + right_bound):
				indices.add(idx)

	# construct array of rt vs (intensity for each mass)

	intensity_array = []
	target_indices = sorted(indices)

	for intensities in im._intensity_array:
		intensities_for_scan = []

		for idx in target_indices:
			intensities_for_scan.append(intensities[idx])

		intensity_array.append(intensities_for_scan)

	# Construct the extracted intensity matrix
	return ExtractedIntensityMatrix(
			time_list=im.time_list,
			mass_list=[im._mass_list[idx] for idx in target_indices],
			intensity_array=intensity_array,
			)

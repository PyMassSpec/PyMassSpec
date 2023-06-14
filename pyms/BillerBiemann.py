"""
Functions to perform Biller and Biemann deconvolution.
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
from typing import List, Sequence, Tuple, Union

# 3rd party
import numpy

# this package
from pyms.IntensityMatrix import BaseIntensityMatrix
from pyms.IonChromatogram import IonChromatogram
from pyms.Peak.Class import Peak
from pyms.Peak.List.Function import is_peak_list
from pyms.Spectrum import MassSpectrum
from pyms.Utils.Utils import _number_types, is_number, is_sequence_of

__all__ = [
		"BillerBiemann",
		"get_maxima_indices",
		"get_maxima_list",
		"get_maxima_list_reduced",
		"get_maxima_matrix",
		"num_ions_threshold",
		"rel_threshold",
		"sum_maxima",
		]

#######################
# structure
# 1) find local maxima per ion, store intensity and scan index
# 2) sum across N scans to compensate for scan type
# 3) sum ions belonging to each maxima scan
#######################


def BillerBiemann(im: BaseIntensityMatrix, points: int = 3, scans: int = 1) -> List[Peak]:
	"""
	Deconvolution based on the algorithm of Biller and Biemann (1974).

	:param im:
	:param points: Number of scans over which to consider a maxima to be a peak.
	:param scans: Number of scans to combine peaks from to compensate for spectra skewing.

	:return: List of detected peaks

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(im, BaseIntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	if not isinstance(points, int):
		raise TypeError("'points' must be an integer")

	if not isinstance(scans, int):
		raise TypeError("'scans' must be an integer")

	rt_list = im.time_list
	mass_list = im.mass_list
	peak_list = []
	maxima_im = get_maxima_matrix(im, points, scans)

	for row_idx, row in enumerate(maxima_im):
		if sum(row) > 0:
			rt = rt_list[row_idx]
			ms = MassSpectrum(mass_list, row)
			peak = Peak(rt, ms)
			peak.bounds = (0, row_idx, 0)  # store IM index for convenience
			# TODO: can the bounds be determined from the intensity matrix?
			peak_list.append(peak)

	return peak_list


def get_maxima_indices(ion_intensities: Union[Sequence, numpy.ndarray], points: int = 3) -> List[int]:
	"""
	Returns the scan indices for the apexes of the ion.

	:param ion_intensities: A list of intensities for a single ion.
	:param points: Number of scans over which to consider a maxima to be a peak.

	:author: Andrew Isaac, Dominic Davis-Foster (type assertions)

	**Example:**

	.. code-block:: python

		>>> # A trivial set of data with two clear peaks
		>>> data = [1, 2, 3, 4, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1]
		>>> get_maxima_indices(data)
		[4, 13]
		>>> # Wider window (more points)
		>>> get_maxima_indices(data, points=10)
		[13]

	"""

	if not is_sequence_of(ion_intensities, _number_types):
		raise TypeError("'ion_intensities' must be a sequence of numbers")

	if not isinstance(points, int):
		raise TypeError("'points' must be an integer")

	# find peak inflection points
	# use a 'points' point window
	# for a plateau after a rise, need to check if it is the left edge of a peak
	peak_point = []
	edge = -1
	points = int(points)
	half = int(points / 2)
	points = 2 * half + 1  # ensure odd number of points

	for index in range(len(ion_intensities) - points + 1):

		left = ion_intensities[index:index + half]
		mid = ion_intensities[index + half]
		right = ion_intensities[index + half + 1:index + points]
		# print(left, mid, right)

		if mid > max(left) and mid > max(right):
			# the max value is in the middle
			peak_point.append(index + half)
			edge = -1  # ignore previous rising edge

		elif mid > max(left) and mid == max(right):
			# start of plateau following rise (left of peak?)
			edge = index + half  # ignore previous rising edge, update latest

		elif mid == max(left) and mid > max(right):
			# start of fall from plateau
			if edge > -1:
				centre = int((edge + index + half) / 2)  # mid point
				peak_point.append(centre)
			edge = -1

	return peak_point


def get_maxima_list(ic: IonChromatogram, points: int = 3) -> List[List[float]]:
	"""
	List of retention time and intensity of local maxima for ion.

	:param ic:
	:param points: Number of scans over which to consider a maxima to be a peak.

	:return: A list of retention time and intensity of local maxima for ion.

	:author: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(ic, IonChromatogram):
		raise TypeError("'ic' must be an IonChromatogram object")

	if not isinstance(points, int):
		raise TypeError("'points' must be an integer")

	peak_point = get_maxima_indices(ic.intensity_array, points)
	mlist = []

	for index in range(len(peak_point)):
		rt = ic.get_time_at_index(peak_point[index])
		intensity = ic.get_intensity_at_index(peak_point[index])
		mlist.append([rt, intensity])

	return mlist


def get_maxima_list_reduced(
		ic: IonChromatogram,
		mp_rt: float,
		points: int = 13,
		window: int = 3,
		) -> List[Tuple[float, float]]:
	"""
	List of retention time and intensity of local maxima for ion.

	| Only peaks around a specific retention time are recorded.
	| Created for use with gap filling algorithm.

	:param ic:
	:param mp_rt: The retention time of the missing peak
	:param points: Number of scans over which to consider a maxima to be a peak.
	:param window: The window around ``mp_rt`` where peaks should be recorded.

	:return: A list of 2-element tuple containing the retention time and
		intensity of local maxima for each ion.

	:author: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(ic, IonChromatogram):
		raise TypeError("'ic' must be an IonChromatogram object")

	if not is_number(mp_rt):
		raise TypeError("'mp_rt' must be a number")

	peak_point = get_maxima_indices(ic.intensity_array, points)
	maxima_list = []

	for index in range(len(peak_point)):
		rt = ic.get_time_at_index(peak_point[index])

		if (rt > float(mp_rt) - window) and (rt < float(mp_rt) + window):
			intensity = ic.get_intensity_at_index(peak_point[index])
			maxima_list.append((rt, intensity))

	return maxima_list


def get_maxima_matrix(im: BaseIntensityMatrix, points: int = 3, scans: int = 1) -> numpy.ndarray:
	"""
	Constructs a matrix containing only data for scans in which particular ions apexed.

	The data can be optionally consolidated into the scan within a range
	with the highest total intensity by adjusting the ``scans`` parameter.
	By default this is ``1``, which does not consolidate the data.

	The columns are ion masses and the rows are scans.
	Get matrix of local maxima for each ion.

	:param im:
	:param points: Number of scans over which to consider a maxima to be a peak.
	:param scans: Number of scans to combine peaks from to compensate for spectra skewing.

	:return: A matrix of giving the intensities of ion masses (columns) and for each scan (rows).

	:author: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(im, BaseIntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	if not isinstance(points, int):
		raise TypeError("'points' must be an integer")

	if not isinstance(scans, int):
		raise TypeError("'scans' must be an integer")

	numrows, numcols = im.size  # scans, masses
	# zeroed matrix, size numrows*numcols
	maxima_im = numpy.zeros((numrows, numcols))
	raw_im = im.intensity_array

	# Construct a 2d array which is all zeros apart from the apexing ions
	for col in range(numcols):  # assume all rows have same width
		# 1st, find maxima
		maxima = get_maxima_indices(raw_im[:, col], points)

		# 2nd, fill intensities
		for row in maxima:
			maxima_im[row, col] = raw_im[row, col]

	# combine spectra within 'scans' scans.
	half = int(scans / 2)

	for row_idx in range(numrows):
		# print(f"{row_idx=}")
		# tic = 0
		best = 0
		loc = 0

		# find best in scans
		for ii in range(scans):
			# print(f"{ii=}")

			# Check the scan window around row_idx is not
			# out of range on left (0) or right (numrows)
			if 0 <= row_idx - half + ii < numrows:
				# print(row_idx - half + ii)
				tic = maxima_im[row_idx - half + ii].sum()

				# find the index of the scan in the window with the highest TIC intensity
				if tic > best:
					best = tic
					loc = ii

		for ii in range(scans):
			# Consolidate data in scan with highest TIC
			source_idx = row_idx - half + ii  # the scan to move data from
			dest_idx = row_idx - half + loc  # the scan to move data into
			if 0 <= source_idx < numrows and ii != loc:
				for col in range(numcols):
					maxima_im[dest_idx, col] += maxima_im[source_idx, col]
					maxima_im[source_idx, col] = 0

	return maxima_im


def num_ions_threshold(
		pl: Sequence[Peak],
		n: int,
		cutoff: float,
		copy_peaks: bool = True,
		) -> List[Peak]:
	"""
	Remove Peaks where there are fewer than ``n`` ions with intensities above the given threshold.

	:param pl:
	:param n: Minimum number of ions that must have intensities above the cutoff.
	:param cutoff: The minimum intensity threshold.
	:param copy_peaks: Whether the returned peak list should contain copies of the peaks.

	:return: A new list of Peak objects.

	:author: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not is_peak_list(pl):
		raise TypeError("'pl' must be a list of Peak objects")

	if not isinstance(n, int):
		raise TypeError("'n' must be an integer")

	if not is_number(cutoff):
		raise TypeError("'cutoff' must be a number")

	if copy_peaks:
		pl = copy.deepcopy(pl)

	new_pl = []
	for p in pl:
		ms = p.mass_spectrum
		ia = ms.mass_spec
		ions = 0
		for i in range(len(ia)):
			if ia[i] >= cutoff:
				ions += 1
		if ions >= n:
			new_pl.append(p)

	return new_pl


def rel_threshold(pl: Sequence[Peak], percent: float = 2, copy_peaks: bool = True) -> List[Peak]:
	"""
	Remove ions with relative intensities less than the given relative
	percentage of the maximum intensity.

	:param pl:
	:param percent: Threshold for relative percentage of intensity.
	:default percent: ``2%``
	:param copy_peaks: Whether the returned peak list should contain copies of the peaks.

	:return: A new list of Peak objects with threshold ions.

	:author: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""  # noqa: D400

	peak_list = pl

	if not is_peak_list(peak_list):
		raise TypeError("'pl' must be a list of Peak objects")
	if not is_number(percent):
		raise TypeError("'percent' must be a number > 0")

	if percent <= 0:
		raise ValueError("'percent' must be a number > 0")

	if copy_peaks:
		peak_list = copy.deepcopy(peak_list)

	new_peak_list = []
	for peak in peak_list:
		ms = peak.mass_spectrum
		# if ms is None:
		# 	raise ValueError("The peak has no mass spectrum.")

		ia = ms.mass_spec
		# assume max(ia) big so /100 1st
		cutoff = (max(ia) / 100.0) * float(percent)
		for i in range(len(ia)):
			if ia[i] < cutoff:
				ia[i] = 0
		ms.mass_spec = ia
		peak.mass_spectrum = ms
		new_peak_list.append(peak)

	return new_peak_list


def sum_maxima(im: BaseIntensityMatrix, points: int = 3, scans: int = 1) -> IonChromatogram:
	"""
	Reconstruct the TIC as sum of maxima.

	:param im:
	:param points: Peak if maxima over 'points' number of scans.
	:param scans: Number of scans to combine peaks from to compensate for spectra skewing.

	:return: The reconstructed TIC.

	:author: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(im, BaseIntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	if not isinstance(points, int):
		raise TypeError("'points' must be an integer")

	if not isinstance(scans, int):
		raise TypeError("'scans' must be an integer")

	maxima_im = get_maxima_matrix(im, points)
	sums = []
	numrows = len(maxima_im)
	half = int(scans / 2)

	for row in range(numrows):
		val = 0
		for ii in range(scans):
			if 0 <= row - half + ii < numrows:
				val += maxima_im[row - half + ii].sum()
		sums.append(val)
	tic = IonChromatogram(numpy.array(sums), im.time_list)

	return tic

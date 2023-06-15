"""
Functions related to Peak modification.
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
from typing import Any, List, Optional, Sequence, Union

# 3rd party
import numpy

# this package
from pyms.IntensityMatrix import BaseIntensityMatrix
from pyms.Peak import Peak
from pyms.Spectrum import MassSpectrum
from pyms.Utils.Math import median_outliers
from pyms.Utils.Time import time_str_secs
from pyms.Utils.Utils import is_sequence, is_sequence_of

__all__ = ["composite_peak", "fill_peaks", "is_peak_list", "sele_peaks_by_rt"]


def composite_peak(peak_list: List[Peak], ignore_outliers: bool = False) -> Optional[Peak]:
	"""
	Create a peak that consists of a composite spectrum from all spectra in the list of peaks.

	:param peak_list: A list of peak objects
	:param ignore_outliers:

	:return: The composite peak

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not is_peak_list(peak_list):
		raise TypeError("'peak_list' must be a list of Peak objects")

	first = True
	count = 0
	avg_rt = 0.0
	# new_ms = None

	# DK: first mark peaks in the list that are outliers by RT, but only if there are more than 3 peaks in the list
	if ignore_outliers:
		rts = []
		if len(peak_list) > 3:
			for peak in peak_list:
				rts.append(peak.rt)

			is_outlier = median_outliers(rts)

			for i, val in enumerate(is_outlier):
				if val:
					peak_list[i].is_outlier = True

	# DK: the average RT and average mass spec for the compound peak is now calculated from peaks that are NOT outliers.
	# This should improve the ability to order peaks and figure out badly aligned entries
	for peak in peak_list:
		if peak is not None and ((ignore_outliers and not peak.is_outlier) or not ignore_outliers):
			ms = peak.mass_spectrum

			if ms is None:
				raise ValueError("The peak has no mass spectrum.")

			spec = numpy.array(ms.mass_spec, dtype='d')
			if first:
				avg_spec = numpy.zeros(len(ms.mass_spec), dtype='d')
				mass_list = ms.mass_list
				first = False
			# scale all intensities to [0,100]
			max_spec = max(spec) / 100.0
			if max_spec > 0:
				spec = spec / max_spec
			else:
				spec = spec * 0
			avg_rt += peak.rt
			avg_spec += spec
			count += 1
	if count > 0:
		avg_rt = avg_rt / count
		avg_spec = avg_spec / count
		new_ms = MassSpectrum(mass_list, avg_spec)

		return Peak(avg_rt, new_ms)
	else:
		return None


def fill_peaks(
		data: BaseIntensityMatrix,
		peak_list: List[Peak],
		D: float,
		minutes: bool = False,
		) -> List[Peak]:
	"""
	Gets the best matching Retention Time and spectra from 'data' for each peak
	in the peak list.

	:param data: A data IntensityMatrix that has the same mass range as the
		peaks in the peak list
	:param peak_list: A list of peak objects
	:param D: Peak width standard deviation in seconds.
		Determines search window width.
	:param minutes: Return retention time as minutes

	:return: List of Peak Objects

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""  # noqa: D400

	if not is_peak_list(peak_list):
		raise TypeError("'peak_list' must be a list of Peak objects")

	if not isinstance(D, float):
		raise TypeError("'D' must be a float")

	# Test for best match in range where RT weight is greater than _TOL
	_TOL = 0.001
	cutoff = D * math.sqrt(-2.0 * math.log(_TOL))

	# Penalise for neighboring peaks
	# reweight so RT weight at nearest peak is _PEN
	_PEN = 0.5

	datamat = data.intensity_array
	mass_list = data.mass_list
	datatimes = data.time_list
	minrt = min(datatimes)
	maxrt = max(datatimes)
	rtl = 0.0
	rtr = 0.0
	new_peak_list = []
	for ii in range(len(peak_list)):
		mass_spec = peak_list[ii].mass_spectrum

		if mass_spec is None:
			raise ValueError("The peak has no mass spectrum.")

		spec = mass_spec.mass_spec
		spec = numpy.array(spec, dtype='d')
		rt = peak_list[ii].rt
		sum_spec_squared = numpy.sum(spec**2, axis=0)

		# get neighbour RT's
		if ii > 0:
			rtl = peak_list[ii - 1].rt
		if ii < len(peak_list) - 1:
			rtr = peak_list[ii + 1].rt
		# adjust weighting for neighbours
		rtclose = min(abs(rt - rtl), abs(rt - rtr))
		Dclose = rtclose / math.sqrt(-2.0 * math.log(_PEN))

		if Dclose > 0:
			Dclose = min(D, Dclose)
		else:
			Dclose = D

		# Get bounds
		rtlow = rt - cutoff
		if rtlow < minrt:
			rtlow = minrt
		lowii = data.get_index_at_time(rtlow)

		rtup = rt + cutoff
		if rtup > maxrt:
			rtup = maxrt
		upii = data.get_index_at_time(rtup)

		# Get sub matrix of scans in bounds
		submat = numpy.array(datamat[lowii:upii + 1], dtype='d')
		subrts = numpy.array(datatimes[lowii:upii + 1], dtype='d')

		sum_summat_squared = numpy.sum(submat**2, axis=1)

		# transpose spec (as matrix) for dot product
		spec = numpy.transpose([spec])
		# dot product on rows

		toparr = numpy.dot(submat, spec)
		botarr = numpy.sqrt(sum_spec_squared * sum_summat_squared)

		# convert back to 1-D array
		toparr = toparr.ravel()

		# scaled dot product of each scan
		cosarr = toparr / botarr

		# RT weight of each scan
		rtimearr = numpy.exp(-((subrts - rt) / float(Dclose))**2 / 2.0)

		# weighted scores
		scorearr = cosarr * rtimearr

		# index of best score
		best_ii = scorearr.argmax()

		# Add new peak
		bestrt = subrts[best_ii]
		bestspec = submat[best_ii].tolist()
		ms = MassSpectrum(mass_list, bestspec)
		new_peak_list.append(Peak(bestrt, ms, minutes))

	return new_peak_list


def is_peak_list(peaks: Any) -> bool:
	"""
	Returns whether ``peaks`` is a valid peak list.

	:author: Dominic Davis-Foster
	"""

	return is_sequence_of(peaks, Peak)


def sele_peaks_by_rt(
		peaks: Union[Sequence[Peak], numpy.ndarray],
		rt_range: Sequence[str],
		) -> List[Peak]:
	"""
	Selects peaks from a retention time range.

	:param peaks: A list of peak objects
	:param rt_range: A list of two time strings, specifying lower and upper retention times.

	:return: A list of peak objects
	"""

	if not is_peak_list(peaks):
		raise TypeError("'peaks' must be a Sequence of Peak objects")

	if not is_sequence(rt_range):
		raise TypeError("'rt_range' must be a Sequence")
	else:
		if len(rt_range) != 2:
			raise ValueError("'rt_range' must have exactly two elements")

		if not isinstance(rt_range[0], str) or not isinstance(rt_range[1], str):
			raise TypeError("lower/upper retention time limits must be strings")

	rt_lo = time_str_secs(rt_range[0])
	rt_hi = time_str_secs(rt_range[1])

	if rt_lo >= rt_hi:
		raise ValueError("lower retention time limit must be less than upper")

	peaks_sele = []

	for peak in peaks:
		rt = peak.rt
		if rt_lo < rt < rt_hi:
			peaks_sele.append(peak)

	# print("%d peaks selected" % (len(peaks_sele)))

	return peaks_sele

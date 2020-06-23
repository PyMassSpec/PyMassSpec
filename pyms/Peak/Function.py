"""
Functions related to Peak modification
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

# stdlib
import copy
from math import ceil
from numbers import Number
from statistics import median

# 3rd party
import deprecation  # type: ignore
from numpy import percentile  # type: ignore

# this package
from pyms import __version__
from pyms.IntensityMatrix import IntensityMatrix
from pyms.Peak import Peak
from pyms.Utils.Utils import is_sequence


def peak_sum_area(im, peak, single_ion=False, max_bound=0):
	"""
	Calculate the sum of the raw ion areas based on detected boundaries.

	:param im: The originating IntensityMatrix object
	:type im: pyms.IntensityMatrix.IntensityMatrix
	:param peak: The Peak object
	:type peak: pyms.Peak.Class.Peak
	:param single_ion: whether single ion areas should be returned, default False
	:type single_ion: bool, optional
	:param max_bound: Optional value to limit size of detected bound, default 0
	:type max_bound: int, optional

	:return: Sum of peak apex ions in detected bounds
	:rtype: float

	:author: Andrew Isaac
	:author: Dominic Davis-Foster (type assertions)
	"""

	# TODO: what's the point of single_ion?

	if not isinstance(im, IntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	if not isinstance(peak, Peak):
		raise TypeError("'peak' must be a Peak object")

	if not isinstance(max_bound, int):
		raise TypeError("'max_bound' must be an integer")

	sum_area = 0
	# Use internal values (not copy)
	# mat = im.matrix_list
	mat = im.intensity_array
	ms = peak.mass_spectrum
	rt = peak.rt
	apex = im.get_index_at_time(rt)

	# get peak masses with non-zero intensity
	mass_ii = [ii for ii in range(len(ms.mass_list)) if ms.mass_spec[ii] > 0]

	area_dict = {}
	# get stats on boundaries
	for ii in mass_ii:
		# get ion chromatogram as list
		ia = [mat[scan][ii] for scan in range(len(mat))]
		area, left, right, l_share, r_share = ion_area(ia, apex, max_bound)
		# need actual mass for single ion areas
		actual_mass = ms.mass_list[ii]
		area_dict[actual_mass] = area
		sum_area += area

	if single_ion:
		return sum_area, area_dict
	else:
		return sum_area


def peak_pt_bounds(im, peak):
	"""
	Approximate the peak bounds (left and right offsets from apex).

	:param im: The originating IntensityMatrix object
	:type im: pyms.IntensityMatrix.IntensityMatrix
	:param peak: The Peak object
	:type peak: pyms.Peak.Class.Peak

	:return: Sum of peak apex ions in detected bounds
	:rtype: tuple of ints

	:author: Andrew Isaac
	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster
	"""

	if not isinstance(im, IntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	if not isinstance(peak, Peak):
		raise TypeError("'peak' must be a Peak object")

	# Use internal values (not copy)
	# mat = im.matrix_list
	mat = im.intensity_array
	ms = peak.mass_spectrum
	rt = peak.rt
	apex = im.get_index_at_time(rt)

	# get peak masses with non-zero intensity
	mass_ii = [ii for ii in range(len(ms.mass_list)) if ms.mass_spec[ii] > 0]

	left_list = []
	right_list = []

	# get stats on boundaries
	for ii in mass_ii:
		# get ion chromatogram as list
		ia = [mat[scan][ii] for scan in range(len(mat))]
		area, left, right, l_share, r_share = ion_area(ia, apex, 0)
		left_list.append(left)
		right_list.append(right)

	left_list.sort()
	right_list.sort()

	return int(ceil(percentile(left_list, 95))), int(ceil(percentile(right_list, 95)))


def peak_top_ion_areas(im, peak, n_top_ions=5, max_bound=0):
	"""
	Calculate and return the ion areas of the five most
		abundant ions in the peak

	:param im: The originating IntensityMatrix object
	:type im: pyms.IntensityMatrix.IntensityMatrix
	:param peak: The Peak object
	:type peak: pyms.Peak.Class.Peak
	:param n_top_ions: Number of top ions to return areas for, default 5
	:type n_top_ions: int, optional
	:param max_bound: Optional value to limit size of detected bound, default 0
	:type max_bound: int, optional

	:return: Dictionary of ion:ion_area pairs
	:rtype: dict

	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(im, IntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	if not isinstance(peak, Peak):
		raise TypeError("'peak' must be a Peak object")

	if not isinstance(n_top_ions, int):
		raise TypeError("'n_top_ions' must be an integer")

	if not isinstance(max_bound, int):
		raise TypeError("'max_bound' must be an integer")

	# ms = peak.mass_spectrum
	rt = peak.rt
	apex = im.get_index_at_time(rt)

	ion_areas = {}  # Dictionary to store ion:ion_area pairs

	top_ions = peak.top_ions(n_top_ions)
	# print(top_ions)

	for ion in top_ions:
		ion_chrom = im.get_ic_at_mass(ion)
		# need ia as a list not numpy array so use .tolist()
		ia = ion_chrom.intensity_array.tolist()
		area, left, right, l_share, r_share = ion_area(ia, apex, max_bound)
		# need actual mass for single ion areas
		ion_areas[ion] = area

	return ion_areas


@deprecation.deprecated(deprecated_in="2.0.0", removed_in="2.2.0",
						current_version=__version__,
						details="Use :func:`pyms.Peak.Function.top_ions_v2` instead")
def top_ions_v1(peak, num_ions=5):
	"""
	Computes the highest 5 intensity ions

	:param peak: the peak to be processed
	:type peak: pyms.Peak.Class.Peak
	:param num_ions: The number of ions to be recorded, default 5
	:type num_ions: int, optional

	:return: A list of the top 5 highest intensity ions
	:rtype: list

	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(peak, Peak):
		raise TypeError("'peak' must be a Peak object")

	if not isinstance(num_ions, int):
		raise TypeError("'n_top_ions' must be an integer")

	intensity_list = peak.mass_spectrum.mass_spec
	mass_list = peak.mass_spectrum.mass_list

	intensity_list_sorted = copy.deepcopy(intensity_list)
	intensity_list_sorted.sort()

	top_ions = []
	top_intensities = intensity_list_sorted[-num_ions:]

	for i in range(len(intensity_list)):
		if intensity_list[i] in top_intensities:
			top_ions.append(mass_list[i])

	return top_ions


@deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
						current_version=__version__,
						details="Use :meth:`pyms.Peak.Class.Peak.top_ions` instead")
def top_ions_v2(peak, num_ions=5):
	"""
	Computes the highest #num_ions intensity ions.

	:param peak: The peak to be processed
	:type peak: pyms.Peak.Class.Peak
	:param num_ions: The number of ions to be recorded, default 5
	:type num_ions: int, optional

	:return: A list of the num_ions highest intensity ions
	:rtype: list

	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(peak, Peak):
		raise TypeError("'peak' must be a Peak object")
	if not isinstance(num_ions, int):
		raise TypeError("'n_top_ions' must be an integer")

	intensity_list = peak.mass_spectrum.mass_spec
	mass_list = peak.mass_spectrum.mass_list

	ic_tuple = zip(intensity_list, mass_list)

	sorted_ic = sorted(ic_tuple)
	top_ic = sorted_ic[-num_ions:]

	top_ions = []

	for entry in top_ic:
		top_ions.append(entry[1])

	return top_ions


def ion_area(ia, apex, max_bound=0, tol=0.5):
	"""
	Find bounds of peak by summing intensities until change in sum is less than 'tol' percent of the current area.

	:param ia: List of intensities for a given mass
	:type ia: list
	:param apex: Index of the peak apex.
	:type apex: int
	:param max_bound: Optional value to limit size of detected bound, default 0
	:type max_bound: int, optional
	:param tol: Percentage tolerance of added area to current area, default 0.5
	:type tol: float, optional

	:return: Area, left and right boundary offset, shared left, shared right
	:rtype: tuple

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(ia, list) or not isinstance(ia[0], Number):
		raise TypeError("'ia' must be a list of numbers")
	if not isinstance(apex, int):
		raise TypeError("'apex' must be an integer")
	if not isinstance(max_bound, int):
		raise TypeError("'max_bound' must be an integer")
	if not isinstance(tol, float):
		raise TypeError("'tol' must be a float")

	# Left area
	lhs = ia[:apex + 1]
	lhs.reverse()  # reverse, as search to right is bounds safe
	l_area, left, l_share = half_area(lhs, max_bound, tol)

	# Right area
	rhs = ia[apex:]
	r_area, right, r_share = half_area(rhs, max_bound, tol)
	r_area -= ia[apex]  # counted apex twice for tollerence, now ignore

	# Put it all together
	return l_area + r_area, left, right, l_share, r_share


def half_area(ia, max_bound=0, tol=0.5):
	"""
	Find bound of peak by summing intensities until change in sum is less than
	'tol' percent of the current area.

	:param ia: List of intensities from Peak apex for a given mass
	:type ia: list
	:param max_bound: Optional value to limit size of detected bound, default 0
	:type max_bound: int, optional
	:param tol: Percentage tolerance of added area to current area, default 0.5
	:type tol: float, optional

	:return: Half peak area, boundary offset, shared (True if shared ion)
	:rtype: tuple

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(ia, list) or not isinstance(ia[0], Number):
		raise TypeError("'ia' must be a list of numbers")
	if not isinstance(max_bound, int):
		raise TypeError("'max_bound' must be an integer")
	if not isinstance(tol, float):
		raise TypeError("'tol' must be a float")

	tol = tol / 200.0  # halve and convert from percent

	# Default number of points to sum new area across, for smoothing
	wide = 3

	# start at 0, compare average value of 'wide' points to the right,
	# centre 'wide' points on edge point,
	# and keep moving right until:
	# i) tollerence reached
	# ii) edge area starts increasing
	# iii) bound reached

	#
	# initialise areas and bounds
	shared = False
	area = ia[0]
	edge = float(sum(ia[0:wide])) / wide
	old_edge = 2 * edge  # bigger than expected edge
	index = 1
	if max_bound < 1:
		limit = len(ia)
	else:
		limit = min(max_bound + 1, len(ia))
	# while edge > area * tol and edge < old_edge and index < limit:
	while area * tol < edge < old_edge and index < limit:
		old_edge = edge
		area += ia[index]
		edge = float(sum(ia[index:index + wide])) / wide  # bounds safe
		index += 1
	if edge >= old_edge:
		shared = True
	index -= 1

	return area, index, shared


def median_bounds(im, peak, shared=True):
	"""
	Calculates the median of the left and right bounds found for each apexing peak mass

	:param im: The originating IntensityMatrix object
	:type im: pyms.IntensityMatrix.IntensityMatrix
	:param peak: The Peak object
	:type peak: pyms.Peak.Class.Peak
	:param shared: Include shared ions shared with neighbouring peak, default True
	:type shared: bool, optional

	:return: Median left and right boundary offset in points
	:rtype: tuple

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(im, IntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")
	if not isinstance(peak, Peak):
		raise TypeError("'peak' must be a Peak object")
	if not isinstance(shared, bool):
		raise TypeError("'shared' must be a boolean")

	mat = im.intensity_array
	ms = peak.mass_spectrum
	rt = peak.rt
	apex = im.get_index_at_time(rt)
	# check if RT based index is similar to stored index
	tmp = peak.bounds
	if is_sequence(tmp) and apex - 1 < tmp[1] < apex + 1:
		apex = tmp[1]

	# get peak masses with non-zero intensity
	mass_ii = [ii for ii in range(len(ms.mass_list)) if ms.mass_spec[ii] > 0]

	# get stats on boundaries
	left_list = []
	right_list = []
	for ii in mass_ii:
		# get ion chromatogram as list
		ia = [mat[scan][ii] for scan in range(len(mat))]
		area, left, right, l_share, r_share = ion_area(ia, apex)
		if shared or not l_share:
			left_list.append(left)
		if shared or not r_share:
			right_list.append(right)

	# return medians
	# NB if shared=True, lists maybe empty
	l_med = 0
	r_med = 0
	if len(left_list) > 0:
		l_med = median(left_list)
	if len(right_list) > 0:
		r_med = median(right_list)

	return l_med, r_med

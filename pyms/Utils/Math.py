"""
Provides mathematical functions.
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                              #
#                                                                              #
#    is_float from 'jcamp' by Nathan Hagen								       #
# 	 https://github.com/nzhagen/jcamp										   #
# 	 Licensed under the X11 License											   #
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
from statistics import mean, median
from statistics import stdev as std
from typing import Iterable, List, Sequence, Union, overload

# 3rd party
import numpy

# this package
from pyms.Utils.Utils import is_number, is_sequence

__all__ = [
		"vector_by_step",
		"MAD",
		"rmsd",
		"mad_based_outlier",
		"percentile_based_outlier",
		"median_outliers",
		"is_float",
		"mean",
		"median",
		"std",
		]


def vector_by_step(start: float, stop: float, step: float) -> List[float]:
	"""
	Generates a list by using start, stop, and step values.

	:param start: Initial value
	:param stop: Max value
	:param step: Step

	:author: Vladimir Likic
	"""  # noqa: D400

	if not is_number(start) or not is_number(stop) or not is_number(step):
		raise TypeError("parameters 'start', 'stop', and 'step' must be numbers")

	v = []

	p = start
	while p < stop:
		v.append(p)
		p = p + step

	return v


def MAD(v: Union[Sequence, numpy.ndarray]) -> float:
	"""
	Median absolute deviation.

	:param v: List of values to calculate the median absolute deviation of.

	:return: median absolute deviation

	:author: Vladimir Likic
	"""

	if not is_sequence(v):
		raise TypeError("'v' must be a Sequence")

	m = median(v)
	m_list = []

	for xi in v:
		d = math.fabs(xi - m)
		m_list.append(d)

	mad = median(m_list) / 0.6745

	return mad


def rmsd(list1: Union[Sequence, numpy.ndarray], list2: Union[Sequence, numpy.ndarray]) -> float:
	"""
	Calculates RMSD for the 2 lists.

	:param list1: First data set
	:param list2: Second data set

	:return: RMSD value

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
	"""

	if not is_sequence(list1):
		raise TypeError("'list1' must be a Sequence")

	if not is_sequence(list2):
		raise TypeError("'list2' must be a Sequence")

	total = 0.0
	for i in range(len(list1)):
		total = total + (list1[i] - list2[i])**2
	_rmsd = math.sqrt(total / len(list1))
	return _rmsd


def mad_based_outlier(data: Iterable[float], thresh: float = 3.5):
	"""
	Identify outliers using the median absolute deviation (MAD).

	:param data:
	:param thresh:

	:author: David Kainer
	:url: http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
	"""

	data = numpy.array(data)
	if len(data.shape) == 1:
		data = data[:, None]
	_median = numpy.nanmedian(data)
	diff = numpy.nansum((data - _median)**2, dtype=float, axis=-1)
	diff = numpy.sqrt(diff)
	med_abs_deviation = numpy.nanmedian(diff)

	modified_z_score = 0.6745 * diff / med_abs_deviation

	return modified_z_score > thresh


def percentile_based_outlier(data: Iterable[float], threshold: int = 95):
	"""
	Identify outliers using a percentile.

	:param data:
	:param threshold:

	:author: David Kainer
	:url: http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
	"""

	data = numpy.array(data)
	diff = (100 - threshold) / 2.0
	# nanpercentile only works in numpy 1.9 and up
	# minval, maxval = numpy.nanpercentile(data, [diff, 100 - diff])
	data = numpy.array(data)
	minval, maxval = numpy.percentile(numpy.compress(numpy.isnan(data) is False, data), (diff, 100 - diff))
	return (data < minval) | (data > maxval)


def median_outliers(data, m: float = 2.5):
	"""
	Identify outliers using the median value.

	:param data:
	:param m:

	:author: David Kainer
	:author: eumiro (https://stackoverflow.com/users/449449/eumiro)
	:author: Benjamin Bannier (https://stackoverflow.com/users/176922/benjamin-bannier)
	:url: http://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
	"""

	data = numpy.array(data)
	d = numpy.abs(data - numpy.nanmedian(data))
	mdev = numpy.nanmedian(d)
	s = d / mdev if mdev else 0.
	return s > m


@overload
def is_float(s: str) -> bool: ...


@overload
def is_float(s: List[str]) -> List[bool]: ...


def is_float(s: Union[str, List[str]]) -> Union[bool, List[bool]]:
	"""
	Test if a string, or list of strings, contains a numeric value(s).

	:param s: The string or list of strings to test.

	:return: A single boolean or list of boolean values indicating whether each input can be converted into a float.
	"""

	if isinstance(s, (tuple, list)):
		if not all(isinstance(i, str) for i in s):
			raise TypeError(f"Input {s} is not a list of strings")

		if len(s) == 0:
			raise ValueError(f'Input {s} is empty')
		else:
			return_list = [True] * len(s)
			for i in range(0, len(s)):
				try:
					float(s[i])
				except ValueError:
					return_list[i] = False
		return return_list
	else:
		if not isinstance(s, str):
			raise TypeError(f"Input '{s}' is not a string")

		try:
			float(s)
			return True
		except ValueError:
			return False

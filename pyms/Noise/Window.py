"""
Moving window noise filter.
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
from statistics import median
from typing import Union

# 3rd party
import numpy

# this package
from pyms.GCMS.Function import ic_window_points
from pyms.IntensityMatrix import BaseIntensityMatrix
from pyms.IonChromatogram import IonChromatogram
from pyms.Noise.SavitzkyGolay import _IM

__all__ = ["window_smooth", "window_smooth_im"]

_DEFAULT_WINDOW = 3


def window_smooth(
		ic: IonChromatogram,
		window: Union[int, str] = _DEFAULT_WINDOW,
		use_median: bool = False,
		) -> IonChromatogram:
	"""
	Applies window smoothing on ion chromatogram.

	:param ic:
	:param window: The window selection parameter. This can be an integer
		or time string. If an integer, taken as the number of points. If a
		string, must be in the form ``'<NUMBER>s'`` or ``'<NUMBER>m'``, specifying
		a time in seconds or minutes, respectively
	:param use_median: Whether to use the the mean or median window smoothing.

	:return: Smoothed ion chromatogram

	:authors: Vladimir Likic, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(ic, IonChromatogram):
		raise TypeError("'ic' must be an IonChromatogram object")

	if not isinstance(window, (int, str)):
		raise TypeError("'window' must be a int or string")

	if not isinstance(use_median, bool):
		raise TypeError("'median' must be a Boolean")

	ia = ic.intensity_array

	wing_length = ic_window_points(ic, window, half_window=True)

	if use_median:
		ia_denoise = _median_window(ia, wing_length)
	else:
		ia_denoise = _mean_window(ia, wing_length)

	ic_denoise = copy.deepcopy(ic)
	ic_denoise.intensity_array = ia_denoise

	return ic_denoise


def window_smooth_im(
		im: _IM,
		window: Union[int, str] = _DEFAULT_WINDOW,
		use_median: bool = False,
		) -> _IM:
	"""
	Applies window smoothing on Intensity Matrix.

	Simply wraps around the window smooth function above.

	:param im:
	:type im: :class:`~.BaseIntensityMatrix`
	:param window: The window selection parameter.
	:param use_median: If :py:obj:`True` median window smoothing will be used.
		If :py:obj:`False` mean window smoothing will be used.

	:return: Smoothed Intensity Matrix
	:rtype: :class:`~.BaseIntensityMatrix`

	:authors: Sean O'Callaghan, Vladimir Likic
	"""

	if not isinstance(im, BaseIntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	n_scan, n_mz = im.size

	im_smooth = copy.deepcopy(im)

	for ii in range(n_mz):
		ic = im_smooth.get_ic_at_index(ii)
		ic_smooth = window_smooth(ic, window, use_median)
		im_smooth.set_ic_at_index(ii, ic_smooth)

	return im_smooth


def _mean_window(ia: numpy.ndarray, wing_length: int) -> numpy.ndarray:
	"""
	Applies mean-window averaging on the array of intensities.

	:param ia: Intensity array
	:param wing_length: The number of points on either side of a point
		in the ion chromatogram.

	:return: Smoothed intensity array

	:author: Vladimir Likic
	"""

	ia_denoise = numpy.repeat([0], ia.size)

	index = 0
	end = ia.size - 1

	while index <= end:
		left = index - wing_length
		right = index + wing_length + 1
		if left < 0:
			left = 0
		ia_denoise[index] = ia[left:right].mean()
		index = index + 1

	return ia_denoise


def _median_window(ia: numpy.ndarray, wing_length: int) -> numpy.ndarray:
	"""
	Applies median-window averaging on the array of intensities.

	:param ia: Intensity array
	:param wing_length: An integer value representing the number of
		points on either side of a point in the ion chromatogram

	:return: Smoothed intensity array

	:author: Vladimir Likic
	"""

	ia_denoise = numpy.repeat([0], ia.size)

	index = 0
	end = ia.size - 1

	while index <= end:
		left = index - wing_length
		right = index + wing_length + 1
		if left < 0:
			left = 0
		ia_denoise[index] = median(ia[left:right])
		index = index + 1

	return ia_denoise

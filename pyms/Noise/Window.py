"""
Moving window noise filter
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
from statistics import median

# 3rd party
from typing import Union

import numpy  # type: ignore

# this package
from pyms.GCMS.Function import ic_window_points
from pyms.IntensityMatrix import IntensityMatrix
from pyms.IonChromatogram import IonChromatogram


__DEFAULT_WINDOW = 3


def window_smooth(ic: IonChromatogram, window: Union[int,str] = __DEFAULT_WINDOW, use_median: bool = False) -> IonChromatogram:
    """
    Applies window smoothing on ion chromatogram

    :param ic: The input ion chromatogram
    :type ic: pyms.IonChromatogram.IonChromatogram
    :param window: The window selection parameter. This can be an integer
        or time string. If integer, taken as the number of points. If a
        string, must of the form "<NUMBER>s" or "<NUMBER>m", specifying
        a time in seconds or minutes, respectively
    :type window: int or str, optional
    :param use_median: An indicator whether the mean or median window smoothing
        to be used
    :type use_median: bool, optional

    :return: Smoothed ion chromatogram
    :rtype: pyms.IonChromatogram.IonChromatogram

    :author: Vladimir Likic
    :author: Dominic Davis-Foster (type assertions)
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
        ia_denoise = __median_window(ia, wing_length)
    else:
        ia_denoise = __mean_window(ia, wing_length)

    ic_denoise = copy.deepcopy(ic)
    ic_denoise.intensity_array = ia_denoise

    return ic_denoise


def window_smooth_im(im: IntensityMatrix, window: Union[int,str] = __DEFAULT_WINDOW, use_median: bool = False) -> IntensityMatrix:
    """
    Applies window smoothing on Intensity Matrix

              Simply wraps around the window smooth function above

    :param im: The input Intensity Matrix
    :type im: pyms.IntensityMatrix.IntensityMatrix
    :param window: The window selection parameter.
    :type window: int or str, optional
    :param use_median: An indicator whether the mean or median window smoothing
        to be used
    :type use_median: bool, optional

    :return: Smoothed Intensity Matrix
    :rtype: pyms.IntensityMatrix.IntensityMatrix

    :author: Sean O'Callaghan
    :author: Vladimir Likic
    """

    if not isinstance(im, IntensityMatrix):
        raise TypeError("'im' must be an IntensityMatrix object")

    n_scan, n_mz = im.size

    im_smooth = copy.deepcopy(im)

    for ii in range(n_mz):
        ic = im_smooth.get_ic_at_index(ii)
        ic_smooth = window_smooth(ic, window, use_median)
        im_smooth.set_ic_at_index(ii, ic_smooth)

    return im_smooth


def __mean_window(ia: numpy.core.ndarray, wing_length: int) -> numpy.core.ndarray:
    """
    Applies mean-window averaging on the array of intensities.

    :param ia: Intensity array
    :type ia: numpy.core.ndarray
    :param wing_length: An integer value representing the number of
        points on either side of a point in the ion chromatogram
    :type wing_length: int

    :return: Smoothed intensity array
    :rtype: numpy.core.ndarray

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


def __median_window(ia: IntensityMatrix, wing_length: int) -> numpy.ndarray:
    """
    Applies median-window averaging on the array of intensities.

    :param ia: Intensity array
    :type ia: numpy.core.ndarray
    :param wing_length: An integer value representing the number of
        points on either side of a point in the ion chromatogram
    :type wing_length: int

    :return: Smoothed intensity array
    :rtype: numpy.core.ndarray

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

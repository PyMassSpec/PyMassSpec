"""
Savitzky-Golay noise filter.
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
from typing import TypeVar, Union

# 3rd party
import numpy

# this package
from pyms.GCMS.Function import ic_window_points
from pyms.IntensityMatrix import BaseIntensityMatrix
from pyms.IonChromatogram import IonChromatogram

__all__ = ["savitzky_golay", "savitzky_golay_im"]

_DEFAULT_WINDOW = 7
_DEFAULT_POLYNOMIAL_DEGREE = 2
_IM = TypeVar("_IM", bound=BaseIntensityMatrix)


def savitzky_golay(
		ic: IonChromatogram,
		window: Union[int, str] = _DEFAULT_WINDOW,
		degree: int = _DEFAULT_POLYNOMIAL_DEGREE,
		) -> IonChromatogram:
	"""
	Applies Savitzky-Golay filter on an ion chromatogram.

	:param ic: The input ion chromatogram.
	:param window: The window selection parameter. This can be an integer
		or time string. If an integer, taken as the number of points. If a
		string, must be the form ``'<NUMBER>s'`` or ``'<NUMBER>m'``, specifying
		a time in seconds or minutes, respectively.
	:param degree: degree of the fitting polynomial for the Savitzky-Golay filter.

	:return: Smoothed ion chromatogram.

	:authors: Uwe Schmitt, Vladimir Likic, Dominic Davis-Foster
	"""

	if not isinstance(ic, IonChromatogram):
		raise TypeError("'ic' must be an IonChromatogram object")

	if not isinstance(window, (int, str)):
		raise TypeError("'window' must be either an int or a string")

	if not isinstance(degree, int):
		raise TypeError("'degree' must be an integer")

	ia = ic.intensity_array

	wing_length = ic_window_points(ic, window, half_window=True)

	# print(" -> Applying Savitzky-Golay filter")
	# print("      Window width (points): %d" % ( 2*wing_length+1 ))
	# print("      Polynomial degree: %d" % ( degree ))

	coeff = _calc_coeff(wing_length, degree)
	ia_denoise = _smooth(ia, coeff)

	ic_denoise = copy.deepcopy(ic)
	ic_denoise.intensity_array = ia_denoise

	return ic_denoise


def savitzky_golay_im(
		im: _IM,
		window: Union[int, str] = _DEFAULT_WINDOW,
		degree: int = _DEFAULT_POLYNOMIAL_DEGREE,
		) -> _IM:
	"""
	Applies Savitzky-Golay filter on Intensity Matrix.

	Simply wraps around the Savitzky Golay function above.

	:param im:
	:type im: :class:`~.BaseIntensityMatrix`
	:param window: The window selection parameter.
	:param degree: degree of the fitting polynomial for the Savitzky-Golay filter.

	:return: Smoothed IntensityMatrix.
	:rtype: :class:`~.BaseIntensityMatrix`

	:authors: Sean O'Callaghan, Vladimir Likic, Dominic Davis-Foster
	"""

	if not isinstance(im, BaseIntensityMatrix):
		raise TypeError("'im' must be an IntensityMatrix object")

	if not isinstance(window, (int, str)):
		raise TypeError("'window' must be either an int or a string")

	if not isinstance(degree, int):
		raise TypeError("'degree' must be an integer")

	n_scan, n_mz = im.size

	im_smooth = copy.deepcopy(im)

	for ii in range(n_mz):
		ic = im_smooth.get_ic_at_index(ii)
		ic_smooth = savitzky_golay(ic, window, degree)
		im_smooth.set_ic_at_index(ii, ic_smooth)

	return im_smooth


def _calc_coeff(num_points: int, pol_degree: int, diff_order: int = 0) -> numpy.ndarray:
	"""
	Calculates filter coefficients for symmetric savitzky-golay filter.

	.. seealso::

		Section 14.8: Savitzky-Golay Smoothing Filters in
		Numerical Recipes in C, Second Edition (1992)
		by Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P.

		Published by Cambridge University Press

	:param num_points: Means that 2*num_points+1 values contribute to the smoother
	:param pol_degree: The degree of fitting polynomial
	:param diff_order: The degree of implicit differentiation.  0 means
		that filter results in smoothing of function, 1 means that filter
		results in smoothing the first derivative of function, and so on.
		Always use ``0``

	:return: Filter coefficients

	:author: Uwe Schmitt
	"""

	# setup normal matrix
	A = numpy.zeros((2 * num_points + 1, pol_degree + 1), float)
	for i in range(2 * num_points + 1):
		for j in range(pol_degree + 1):
			A[i, j] = pow(i - num_points, j)

	# calculate diff_order-th row of inv(A^T A)
	ATA = numpy.dot(A.transpose(), A)
	rhs = numpy.zeros((pol_degree + 1), float)
	rhs[diff_order] = 1
	D = numpy.linalg.cholesky(ATA)
	wvec = _resub(D, rhs)

	# calculate filter-coefficients
	coeff = numpy.zeros((2 * num_points + 1), float)
	for n in range(-num_points, num_points + 1):
		x = 0.0
		for m in range(pol_degree + 1):
			x += wvec[m] * pow(n, m)
		coeff[n + num_points] = x

	return coeff


def _resub(D, rhs):
	"""
	Solves ``D D^T = rhs`` by resubstitution.

	D is lower triangle-matrix from cholesky-decomposition.

	:param D:
	:param rhs:

	:author: Uwe Schmitt
	"""

	M = D.shape[0]
	x1 = numpy.zeros((M, ), float)
	x2 = numpy.zeros((M, ), float)

	# resub step 1
	for l in range(M):
		total = rhs[l]
		for n in range(l):
			total -= D[l, n] * x1[n]
		x1[l] = total / D[l, l]

	# resub step 2
	for l in range(M - 1, -1, -1):
		total = x1[l]
		for n in range(l + 1, M):
			total -= D[n, l] * x2[n]
		x2[l] = total / D[l, l]

	return x2


def _smooth(signal: numpy.ndarray, coeff: numpy.ndarray) -> numpy.ndarray:
	"""
	Applies coefficients calculated by :func:`~._calc_coeff()` to signal.

	:param signal:
	:param coeff:

	:author: Uwe Schmitt
	"""

	size = numpy.size(coeff - 1) // 2
	res = numpy.convolve(signal, coeff)
	return res[size:-size]

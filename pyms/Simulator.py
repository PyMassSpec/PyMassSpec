"""
Provides functions for simulation of GCMS data.
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
from typing import List

# 3rd party
import numpy

# this package
from pyms import Peak
from pyms.IntensityMatrix import BaseIntensityMatrix, IntensityMatrix
from pyms.IonChromatogram import IonChromatogram

__all__ = [
		"add_gaussc_noise",
		"add_gaussc_noise_ic",
		"add_gaussv_noise",
		"add_gaussv_noise_ic",
		"chromatogram",
		"gaussian",
		"gcms_sim",
		]


def add_gaussc_noise(im: BaseIntensityMatrix, scale: float) -> None:
	"""
	Adds noise to an IntensityMatrix object.

	:param im: the intensity matrix object
	:param scale: the scale of the normal distribution from which the noise is drawn

	:author: Sean O'Callaghan
	"""

	n_scan, n_mz = im.size

	for i in range(n_mz):
		ic = im.get_ic_at_index(i)
		add_gaussc_noise_ic(ic, scale)
		im.set_ic_at_index(i, ic)


def add_gaussc_noise_ic(ic: IonChromatogram, scale: float) -> None:
	"""
	Adds noise drawn from a normal distribution with constant scale to an ion chromatogram.

	:param ic: The ion Chromatogram.
	:param scale: The scale of the normal distribution.

	:author: Sean O'Callaghan
	"""

	noise = numpy.random.normal(0.0, scale, (len(ic)))

	i_array_with_noise = ic.intensity_array + noise
	ic.intensity_array = i_array_with_noise


def add_gaussv_noise(
		im: BaseIntensityMatrix,
		scale: int,
		cutoff: int,
		prop: float,
		) -> None:
	"""
	Adds noise to an IntensityMatrix object.

	:param im: the intensity matrix object
	:param scale: the scale of the normal distribution from which the noise is drawn
	:param cutoff: The level below which the intensity of the ic at that point
		has no effect on the scale of the noise distribution
	:param scale: The scale of the normal distribution for ic values
	:param prop: For intensity values above the cutoff, the scale is
		multiplied by the ic value multiplied by ``prop``.

	:author: Sean O'Callaghan
	"""

	n_scan, n_mz = im.size

	for i in range(n_mz):
		ic = im.get_ic_at_index(i)
		add_gaussv_noise_ic(ic, scale, cutoff, prop)
		im.set_ic_at_index(i, ic)


def add_gaussv_noise_ic(
		ic: IonChromatogram,
		scale: int,
		cutoff: int,
		prop: float,
		) -> None:
	"""
	Adds noise to an ic. The noise value is drawn from a normal
	distribution, the scale of this distribution depends on the
	value of the ic at the point where the noise is being added

	:param ic: The IonChromatogram
	:param cutoff: The level below which the intensity of the ic at that point
		has no effect on the scale of the noise distribution
	:param scale: The scale of the normal distribution for ic values below the cutoff
		is modified for values above the cutoff
	:param prop: For ic values above the cutoff, the scale is multiplied by the ic
		value multiplied by ``prop``.

	:author: Sean O'Callaghan
	"""  # noqa: D400

	noise = numpy.zeros(len(ic))

	i_array = ic.intensity_array
	# time_list = ic.time_list

	for i in range(len(ic)):
		if i_array[i] < cutoff:
			noise[i] = numpy.random.normal(0.0, scale, 1)
		else:
			noise[i] = numpy.random.normal(0.0, scale * i_array[i] * prop, 1)

	i_array_with_noise = noise + i_array
	ic.intensity_array = i_array_with_noise


def chromatogram(
		n_scan: int,
		x_zero: int,
		sigma: float,
		peak_scale: float,
		) -> numpy.ndarray:
	"""
	Returns a simulated ion chromatogram of a pure component.

	The ion chromatogram contains a single gaussian peak.

	:param n_scan: the number of scans
	:param x_zero: The apex of the peak
	:param sigma: The standard deviation of the distribution
	:param peak_scale: the intensity of the peak at the apex

	:author: Sean O'Callaghan
	"""

	ic = numpy.zeros(n_scan, 'd')

	for i in range(n_scan):
		x = float(i)

		ic[i] = gaussian(x, x_zero, sigma, peak_scale)

	return ic


def gaussian(point: float, mean: int, sigma: float, scale: float) -> float:
	"""
	Calculates a point on a gaussian density function.

	``f = s*exp(-((x-x0)^2)/(2*w^2))``

	:param point: The point currently being computed
	:param mean: The apex of the peak
	:param sigma: The standard deviation of the gaussian
	:param scale: The height of the apex

	:return: a single value from a normal distribution

	:author: Sean O'Callaghan
	"""

	return scale * math.exp((-(point - mean)**2) / (2 * (sigma**2)))


def gcms_sim(
		time_list: List[float],
		mass_list: List[float],
		peak_list: List[Peak.Peak],
		) -> IntensityMatrix:
	"""
	Simulator of GCMS data.

	:param time_list: the list of scan times
	:param mass_list: the list of m/z channels
	:param peak_list: A list of peaks

	:return: A simulated Intensity Matrix object

	:author: Sean O'Callaghan
	"""

	n_mz = len(mass_list)
	n_scan = len(time_list)

	t1 = time_list[0]
	period = time_list[1] - t1

	# initialise a 2D numpy array for intensity matrix
	i_array = numpy.zeros((n_scan, n_mz), 'd')

	for peak in peak_list:
		print('-', end='')
		index = int((peak.rt - t1) / period)

		height = sum(peak.mass_spectrum.mass_spec)
		# standard deviation = area/(height * sqrt(2/pi))
		sigma: float = peak.area / (height * (math.sqrt(2 * math.pi)))
		print("width", sigma)

		for i in range(len(peak.mass_spectrum.mass_list)):
			ion_height = peak.mass_spectrum.mass_spec[i]
			ic = chromatogram(n_scan, index, sigma, ion_height)
			i_array[:, i] += ic

	im = IntensityMatrix(time_list, mass_list, i_array)

	return im

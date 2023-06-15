#############################################################################
#                                                                           #
#    PyMassSpec software for processing of mass-spectrometry data           #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                           #
#                                                                           #
#    This program is free software; you can redistribute it and/or modify   #
#    it under the terms of the GNU General Public License version 2 as      #
#    published by the Free Software Foundation.                             #
#                                                                           #
#    This program is distributed in the hope that it will be useful,        #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#    GNU General Public License for more details.                           #
#                                                                           #
#    You should have received a copy of the GNU General Public License      #
#    along with this program; if not, write to the Free Software            #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#                                                                           #
#############################################################################

# stdlib
import copy
from typing import Any, List

# 3rd party
import numpy
import pytest

# this package
from pyms.BillerBiemann import (
		BillerBiemann,
		get_maxima_indices,
		get_maxima_list,
		get_maxima_list_reduced,
		get_maxima_matrix,
		num_ions_threshold,
		rel_threshold,
		sum_maxima
		)
from pyms.IntensityMatrix import IntensityMatrix
from pyms.IonChromatogram import IonChromatogram
from pyms.Noise.Analysis import window_analyzer
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Class import Peak
from pyms.TopHat import tophat
from tests.constants import *


class TestBillerBiemann:

	def test_BillerBiemann(self, im_i: IntensityMatrix):
		im_i = copy.deepcopy(im_i)
		# Intensity matrix size (scans, masses)
		n_scan, n_mz = im_i.size

		# noise filter and baseline correct
		for ii in range(n_mz):
			ic = im_i.get_ic_at_index(ii)
			ic_smooth = savitzky_golay(ic)
			ic_bc = tophat(ic_smooth, struct="1.5m")
			im_i.set_ic_at_index(ii, ic_bc)

		# Use Biller and Biemann technique to find apexing ions at a scan
		# default is maxima over three scans and not to combine with any neighbouring
		# scan.
		peak_list = BillerBiemann(im_i)
		assert isinstance(peak_list, list)
		assert len(peak_list) == 2101
		for peak in peak_list:
			assert isinstance(peak, Peak)

		# Find apex oven 9 points and combine with neighbouring peak if two scans apex
		# next to each other.
		peak_list2 = BillerBiemann(im_i, points=9, scans=2)
		assert len(peak_list2) == 805

		assert len(peak_list2) <= len(peak_list)

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, *test_sequences, test_dict])
	def test_im_errors(self, obj: Any):
		with pytest.raises(TypeError):
			BillerBiemann(obj)

	@pytest.mark.parametrize("obj", [test_string, test_float, *test_sequences, test_dict])
	def test_scans_errors(self, obj: Any, im_i: IntensityMatrix):
		with pytest.raises(TypeError):
			BillerBiemann(im_i, scans=obj)

	@pytest.mark.parametrize("obj", [test_string, test_float, *test_sequences, test_dict])
	def test_points_errors(self, obj: Any, im_i: IntensityMatrix):
		with pytest.raises(TypeError):
			BillerBiemann(im_i, points=obj)


class Test_rel_threshold:

	def test_rel_threshold(self, peak_list: List[Peak]):
		rel_threshold(peak_list)
		rel_threshold(peak_list, 2.0)
		pl = rel_threshold(peak_list, 2)

		assert isinstance(pl, list)
		assert isinstance(pl[0], Peak)

		assert len(pl) == 805
		assert len(pl) <= len(peak_list)

		# Errors
		with pytest.raises(ValueError, match="percent' must be a number > 0"):
			rel_threshold(peak_list, percent=0)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_dict, test_int])
	def test_peak_list_errors(self, obj: Any):
		with pytest.raises(TypeError):
			rel_threshold(obj)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_dict])
	def test_percent_errors(self, obj: Any, peak_list: List[Peak]):
		with pytest.raises(TypeError):
			rel_threshold(peak_list, percent=obj)


class Test_num_ions_threshold:

	def test_num_ions_threshold(self, peak_list: List[Peak], tic: IonChromatogram):
		"""
		Filter the peak list, first by removing all intensities in a peak less
		than a given relative threshold, then by removing all peaks that have
		less than a given number of ions above a given value
		"""

		# trim by relative intensity
		pl = rel_threshold(peak_list, 2)

		# trim by threshold
		new_peak_list = num_ions_threshold(pl, 3, 10000)
		assert isinstance(new_peak_list, list)
		assert isinstance(new_peak_list[0], Peak)

		assert len(new_peak_list) == 215
		assert len(new_peak_list) <= len(peak_list)
		assert len(new_peak_list) <= len(pl)

		# With window_analyzer
		# estimate noise level from the TIC, used later to
		# discern true signal peaks
		noise_level = window_analyzer(tic)

		# trim by relative intensity
		apl = rel_threshold(peak_list, 1)

		# trim by number of ions above threshold
		peak_list = num_ions_threshold(apl, 3, noise_level)

		assert isinstance(peak_list, list)
		assert isinstance(peak_list[0], Peak)

		assert len(peak_list) in (87, 88)
		assert len(peak_list) <= len(peak_list)

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, *test_sequences, test_dict])
	def test_peak_list_errors(self, obj: Any):
		with pytest.raises(TypeError):
			num_ions_threshold(obj, n=5, cutoff=100)

	@pytest.mark.parametrize("obj", [test_string, test_float, *test_sequences, test_dict])
	def test_n_errors(self, obj: Any, peak_list: List[Peak]):
		with pytest.raises(TypeError):
			num_ions_threshold(peak_list, n=obj, cutoff=100.0)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_dict])
	def test_cutoff_errors(self, obj: Any, peak_list: List[Peak]):
		with pytest.raises(TypeError):
			num_ions_threshold(peak_list, n=5, cutoff=obj)


class Test_sum_maxima:

	def test_sum_maxima(self, im: IntensityMatrix):
		new_tic = sum_maxima(im)
		assert isinstance(new_tic, IonChromatogram)
		assert new_tic.is_tic()

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, *test_sequences, test_dict])
	def test_im_errors(self, obj: Any):
		with pytest.raises(TypeError):
			sum_maxima(obj)

	@pytest.mark.parametrize("obj", [test_string, test_float, *test_sequences, test_dict])
	def test_points_errors(self, obj: Any, im: IntensityMatrix):
		with pytest.raises(TypeError):
			sum_maxima(im, points=obj)

	@pytest.mark.parametrize("obj", [test_string, test_float, *test_sequences, test_dict])
	def test_errors_errors(self, obj: Any, im: IntensityMatrix):
		with pytest.raises(TypeError):
			sum_maxima(im, scans=obj)


class Test_get_maxima_indices:

	@pytest.mark.skip(reason="TODO")
	def test_get_maxima_indices(self):
		# TODO: main test
		pass

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, test_list_strs, test_dict])
	def test_ion_intensities_errors(self, obj: Any):
		with pytest.raises(TypeError):
			get_maxima_indices(obj)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_float, test_dict])
	def test_points_errors(self, obj: Any):
		with pytest.raises(TypeError):
			get_maxima_indices(test_list_ints, points=obj)


class Test_get_maxima_list:

	def test_get_maxima_list(self, tic: Any):
		maxima_iist = get_maxima_list(tic)
		assert isinstance(maxima_iist, list)
		assert isinstance(maxima_iist[0], list)
		assert isinstance(maxima_iist[0][0], float)
		assert maxima_iist[0][0] == 2.10800014436

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, *test_sequences, test_dict])
	def test_ic_errors(self, obj: Any):
		with pytest.raises(TypeError):
			get_maxima_list(obj)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_float, test_dict])
	def test_points_errors(self, obj: Any, tic: IonChromatogram):
		with pytest.raises(TypeError):
			get_maxima_list(tic, points=obj)


class Test_get_maxima_list_reduced:

	def test_get_maxima_list_reduced(self, tic: IonChromatogram):
		maxima_list = get_maxima_list_reduced(tic, 12.34)
		assert isinstance(maxima_list, list)
		for peak in maxima_list:
			assert isinstance(peak, tuple)
			assert len(peak) == 2
			rt, intensity = peak
			assert isinstance(rt, float)
			assert isinstance(intensity, float)
		assert maxima_list[0][0] == 10.5559998751

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, *test_sequences, test_dict])
	def test__errors(self, obj: Any):
		with pytest.raises(TypeError):
			get_maxima_list_reduced(obj, 0)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_dict])
	def test_mp_rt_errors(self, obj: Any, tic: IonChromatogram):
		with pytest.raises(TypeError):
			get_maxima_list_reduced(tic, mp_rt=obj)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_float, test_dict])
	def test_points_errors(self, obj: Any, tic: IonChromatogram):
		with pytest.raises(TypeError):
			get_maxima_list_reduced(tic, test_float, points=obj)

	@pytest.mark.parametrize("obj", [test_string, *test_sequences, test_dict])
	def test_window_errors(self, obj: Any, tic: IonChromatogram):
		with pytest.raises(TypeError):
			get_maxima_list_reduced(tic, test_float, window=obj)


class Test_get_maxima_matrix:

	def test_get_maxima_matrix(self, im: IntensityMatrix):
		maxima_matrix = get_maxima_matrix(im)
		assert isinstance(maxima_matrix, numpy.ndarray)
		# TODO: value check

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, *test_sequences, test_dict])
	def test_im_errors(self, obj: Any):
		with pytest.raises(TypeError):
			get_maxima_matrix(obj)

	@pytest.mark.parametrize("obj", [test_string, test_float, *test_sequences, test_dict])
	def test_points_errors(self, obj: Any, im: IntensityMatrix):
		with pytest.raises(TypeError):
			get_maxima_matrix(im, points=obj)

	@pytest.mark.parametrize("obj", [test_string, test_float, *test_sequences, test_dict])
	def test_scans_errors(self, obj: Any, im: IntensityMatrix):
		with pytest.raises(TypeError):
			get_maxima_matrix(im, scans=obj)

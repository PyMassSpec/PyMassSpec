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
from typing import Any, List, Type

# 3rd party
import pytest
from domdf_python_tools.paths import PathPlus

# this package
from pyms.GCMS.Class import GCMS_data
from pyms.IntensityMatrix import IntensityMatrix
from pyms.Peak.Class import Peak
from pyms.Peak.List import composite_peak, fill_peaks, is_peak_list, sele_peaks_by_rt
from pyms.Peak.List.IO import load_peaks, store_peaks
from pyms.Spectrum import MassSpectrum
from tests.constants import *


def test_composite_peak(filtered_peak_list: List[Peak]):
	print('.', end='')
	composite_peak_list = filtered_peak_list[10:20]
	print('.', end='')
	peak = composite_peak(composite_peak_list)
	print('.', end='')
	assert isinstance(peak, Peak)

	uid = peak.UID
	assert uid == "96-69-62-47.73"

	assert peak.get_third_highest_mz() == 57
	assert peak.bounds is None
	assert peak.get_int_of_ion(100) == 3.603215488138507
	assert peak.rt == 47.727200388899995
	assert peak.top_ions(10)[0] == 115

	# area = peak_sum_area(im_i, peak)
	# peak.area = area
	# assert peak.area == area

	assert isinstance(peak.mass_spectrum, MassSpectrum)
	assert isinstance(peak.mass_spectrum.mass_spec, list)
	peak.null_mass(73)
	index_73 = peak.mass_spectrum.mass_list.index(73)
	assert peak.mass_spectrum.mass_spec[index_73] == 0

	peak.crop_mass(100, 200)
	assert peak.UID != uid

	# Errors
	for obj in [test_dict, *test_sequences, *test_numbers, test_string]:
		with pytest.raises(TypeError):
			composite_peak(obj)  # type: ignore[arg-type]


def test_composite_peak_outliers(filtered_peak_list: List[Peak]):
	composite_peak_list = filtered_peak_list[10:13]
	peak = composite_peak(composite_peak_list, ignore_outliers=True)
	assert isinstance(peak, Peak)

	uid = peak.UID
	assert uid == "88-86-92-39.07"

	assert peak.get_third_highest_mz() == 85
	assert peak.bounds is None
	assert peak.get_int_of_ion(100) == 7.1236965120460285
	assert peak.rt == 39.0680015087
	assert peak.top_ions(10)[0] == 161

	# area = peak_sum_area(im_i, peak)
	# peak.area = area
	# assert peak.area == area

	peak.null_mass(73)
	index_73 = peak.mass_spectrum.mass_list.index(73)
	assert peak.mass_spectrum.mass_spec[index_73] == 0

	peak.crop_mass(100, 200)
	assert peak.UID != uid


def test_fill_peaks(im_i: IntensityMatrix, peak_list: List[Peak]):
	filled_peak_list = fill_peaks(im_i, peak_list, 10.0)
	assert is_peak_list(filled_peak_list)

	# Errors
	for obj in [test_dict, *test_sequences, *test_numbers, test_string]:
		with pytest.raises(TypeError):
			fill_peaks(im_i, obj, 10.0)  # type: ignore[arg-type]
	for obj in [test_dict, *test_sequences, test_string, test_int]:
		with pytest.raises(TypeError):
			fill_peaks(im_i, peak_list, obj)  # type: ignore[arg-type]


def test_is_peak_list(peak_list: List[Peak], ms: MassSpectrum, im_i: IntensityMatrix, data: GCMS_data):
	assert is_peak_list(peak_list)
	assert not is_peak_list(test_int)
	assert not is_peak_list(test_string)
	assert not is_peak_list(test_float)
	assert not is_peak_list(test_list_strs)
	assert not is_peak_list(test_list_ints)
	assert not is_peak_list(test_tuple)
	assert not is_peak_list(test_dict)
	assert not is_peak_list(ms)
	assert not is_peak_list(im_i)
	assert not is_peak_list(data)


def test_sele_peaks_by_rt(filtered_peak_list: List[Peak]):
	selected_peaks = sele_peaks_by_rt(filtered_peak_list, ("12m", "13m"))
	assert is_peak_list(selected_peaks)
	assert len(selected_peaks) == 18
	peak = selected_peaks[0]
	assert isinstance(peak, Peak)

	uid = peak.UID
	assert uid == "68-54-26-722.30"

	assert peak.get_third_highest_mz() == 50
	assert peak.bounds == (0, 683, 0)
	assert peak.get_int_of_ion(100) == 0.0
	assert peak.rt == 722.299976349
	assert peak.top_ions(10)[0] == 133

	peak.null_mass(73)
	index_73 = peak.mass_spectrum.mass_list.index(73)
	assert peak.mass_spectrum.mass_spec[index_73] == 0

	peak.crop_mass(100, 200)
	assert peak.UID != uid

	with pytest.raises(TypeError, match="lower/upper retention time limits must be strings"):
		sele_peaks_by_rt(filtered_peak_list, [1.2, 3.4])  # type: ignore[list-item]
	with pytest.raises(ValueError, match="lower retention time limit must be less than upper"):
		sele_peaks_by_rt(filtered_peak_list, ["50s", "10s"])

	# Errors
	for obj in [test_dict, *test_sequences, *test_numbers, test_string]:
		with pytest.raises(TypeError, match="'peaks' must be a Sequence of Peak objects"):
			sele_peaks_by_rt(obj, ("12m", "13m"))  # type: ignore[arg-type]
	for obj in [test_dict, *test_numbers, test_string]:
		with pytest.raises(TypeError, match="'rt_range' must be a Sequence"):
			sele_peaks_by_rt(filtered_peak_list, obj)  # type: ignore[arg-type]
	for obj in [*test_sequences]:
		with pytest.raises(ValueError, match="'rt_range' must have exactly two elements"):
			sele_peaks_by_rt(filtered_peak_list, obj)  # type: ignore[arg-type]


@pytest.fixture()
def peak_list_filename(filtered_peak_list: List[Peak], tmp_pathplus: PathPlus) -> PathPlus:
	filename = tmp_pathplus / "filtered_peak_list.dat"
	store_peaks(filtered_peak_list, filename)
	return filename


class TestStoreLoadPeaks:

	@pytest.mark.parametrize("obj", [test_dict, *test_sequences, *test_numbers, test_string])
	def test_store_filename_errors(self, tmp_pathplus: PathPlus, obj: Any):
		with pytest.raises(TypeError):
			store_peaks(obj, tmp_pathplus / test_string)

	@pytest.mark.parametrize("obj", [test_dict, *test_sequences, *test_numbers])
	def test_store_peak_list_errors(self, filtered_peak_list: List[Peak], obj: Any):
		with pytest.raises(TypeError):
			store_peaks(filtered_peak_list, obj)

	def test_load_peaks(self, filtered_peak_list: List[Peak], peak_list_filename: PathPlus):
		loaded_peak_list = load_peaks(peak_list_filename)

		assert loaded_peak_list == filtered_peak_list

	@pytest.mark.parametrize("filename", [test_dict, *test_sequences, *test_numbers])
	def test_load_filename_errors_1(self, filename: Any):
		with pytest.raises(TypeError):
			load_peaks(filename)

	@pytest.mark.parametrize(
			"filename, expects",
			[
					(test_string, FileNotFoundError),
					("not-an-experiment.expr", IOError),
					("test_list_ints.dat", IOError),
					("test_empty_list.dat", IOError),
					],
			)
	def test_load_filename_errors_2(self, filename: str, expects: Type[Exception], pyms_datadir: PathPlus):
		with pytest.raises(expects):
			load_peaks(pyms_datadir / filename)

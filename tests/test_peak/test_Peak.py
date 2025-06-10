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
from typing import Any

# 3rd party
import pytest
from domdf_python_tools.paths import PathPlus

# this package
from pyms.IntensityMatrix import IntensityMatrix
from pyms.Peak import Peak
from pyms.Peak.Class import ICPeak
from pyms.Peak.Function import peak_sum_area, top_ions_v1, top_ions_v2
from pyms.Spectrum import MassSpectrum
from pyms.Utils.Utils import _pickle_load_path, is_number
from tests.constants import *


def test_Peak(im_i: IntensityMatrix, peak: Peak):
	assert isinstance(peak, Peak)

	# Get the scan of a known TIC peak (at RT 31.17 minutes)
	# get the index of the scan nearest to 31.17 minutes (converted to seconds)
	scan_i = im_i.get_index_at_time(31.17 * 60.0)
	# get the MassSpectrum Object
	ms = im_i.get_ms_at_index(scan_i)

	# create a Peak object
	Peak(31.17)
	Peak(31.17, ms, outlier=True)
	Peak(31.17, ms, minutes=True)

	# Errors
	for obj in [test_string, *test_lists, test_dict]:

		with pytest.raises(TypeError):
			Peak(obj, ms, minutes=True)  # type: ignore[arg-type]
		with pytest.raises(TypeError):
			Peak(test_float, obj, minutes=False)  # type: ignore[arg-type]

	Peak(test_float, test_int)
	ICPeak(test_float, test_int)
	Peak(test_float, test_float)
	ICPeak(test_float, test_float)


def test_equality(peak: Peak):
	assert peak == Peak(peak.rt, peak.mass_spectrum)
	assert peak != Peak(50, peak.mass_spectrum)


@pytest.mark.parametrize(
		"val", [test_list_ints, test_list_strs, test_tuple, test_string, test_int, test_float, test_dict]
		)
def test_inequality(peak: Peak, val: Any):
	assert peak != val


def test_area(im_i: IntensityMatrix, peak: Peak):
	peak = copy.deepcopy(peak)

	# determine and set area
	area = peak_sum_area(im_i, peak)
	assert isinstance(area, float)
	peak.area = area

	assert peak.area == area
	assert isinstance(peak.area, float)

	scan_i = im_i.get_index_at_time(31.17 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)

	for obj in [test_string, test_dict, test_list_strs, test_list_ints]:
		with pytest.raises(TypeError):
			Peak(test_float, ms).area = obj  # type: ignore[assignment]
	with pytest.raises(ValueError, match="'Peak.area' must be a positive number"):
		Peak(test_float, ms).area = -1


def test_bounds(peak: Peak):
	peak = copy.copy(peak)

	# Setter
	peak.bounds = (11, 12, 13)

	for obj in [test_string, *test_numbers, test_dict, ['a', 'b', 'c'], test_tuple]:
		with pytest.raises(TypeError):
			peak.bounds = obj  # type: ignore[assignment]

	for obj in [*test_lists, (1, 2), [1, 2, 3, 4]]:
		with pytest.raises(ValueError, match="'Peak.bounds' must have exactly 3 elements"):
			peak.bounds = obj  # type: ignore[assignment]

	# Getter
	assert peak.bounds == (11, 12, 13)
	assert isinstance(peak.bounds, tuple)
	peak2 = Peak(test_float)
	peak2.bounds = [11, 12, 13]
	assert peak2.bounds == (11, 12, 13)
	assert isinstance(peak2.bounds, tuple)

	# set_bounds
	peak3 = Peak(test_float)
	peak3.set_bounds(11, 12, 13)
	assert peak3.bounds == (11, 12, 13)
	assert isinstance(peak3.bounds, tuple)

	for obj in [*test_sequences, test_string, test_dict, test_float]:
		print(obj)

		with pytest.raises(TypeError):
			peak3.set_bounds(obj, 12, 13)  # type: ignore[arg-type]
		with pytest.raises(TypeError):
			peak3.set_bounds(11, obj, 13)  # type: ignore[arg-type]
		with pytest.raises(TypeError):
			peak3.set_bounds(11, 12, obj)  # type: ignore[arg-type]


def test_crop_mass(peak: Peak):
	peak = copy.deepcopy(peak)
	peak2 = copy.deepcopy(peak)

	uid = peak.UID
	peak.crop_mass(100, 200)
	assert peak.UID != uid
	assert len(peak.mass_spectrum) == 101
	assert min(peak.mass_spectrum.mass_list) == 100
	assert max(peak.mass_spectrum.mass_list) == 200

	# Errors
	for obj in [test_string, *test_lists, test_dict]:
		with pytest.raises(TypeError, match="'mass_min' and 'mass_max' must be numbers"):
			peak2.crop_mass(obj, 450)  # type: ignore[arg-type]
		with pytest.raises(TypeError, match="'mass_min' and 'mass_max' must be numbers"):
			peak2.crop_mass(450, obj)  # type: ignore[arg-type]

	with pytest.raises(ValueError, match="'mass_min' must be less than 'mass_max'"):
		peak2.crop_mass(100, 0)
	with pytest.raises(ValueError, match="'mass_min' is less than the smallest mass: 50"):
		peak2.crop_mass(10, 450)
	with pytest.raises(ValueError, match="'mass_max' is greater than the largest mass: 499"):
		peak2.crop_mass(60, 500)
	with pytest.warns(Warning):
		peak2.crop_mass(60, 65)


def test_get_int_of_ion(peak: Peak):
	assert peak.get_int_of_ion(100) == 3888.0
	assert peak.get_int_of_ion(200) == 0.0
	assert is_number(peak.get_int_of_ion(100))

	with pytest.raises(IndexError):
		peak.get_int_of_ion(1)
	with pytest.raises(IndexError):
		peak.get_int_of_ion(1000000)


def test_ion_area(peak: Peak):
	peak = copy.deepcopy(peak)

	assert peak.get_ion_area(1) is None

	peak.set_ion_area(1, 1234)
	peak.set_ion_area(2, 1234.56)

	assert is_number(peak.get_ion_area(1))
	assert is_number(peak.get_ion_area(2))
	assert peak.get_ion_area(1) == 1234

	# Errors
	for obj in [test_dict, *test_sequences, test_float, test_string]:
		with pytest.raises(TypeError):
			peak.set_ion_area(obj, test_int)  # type: ignore[arg-type]
	for obj in [test_dict, *test_sequences, test_string]:
		with pytest.raises(TypeError):
			peak.set_ion_area(1, obj)  # type: ignore[arg-type]


def test_ion_areas(peak: Peak):
	peak = copy.deepcopy(peak)

	with pytest.raises(ValueError, match="no ion areas set"):
		peak.ion_areas

	peak.ion_areas = {1: 1234, 2: 1234, 3: 1234}

	for obj in [*test_numbers, test_string, test_list_strs, test_list_ints, tuple]:
		with pytest.raises(TypeError):
			peak.ion_areas = obj  # type: ignore[assignment]

	assert peak.ion_areas == {1: 1234, 2: 1234, 3: 1234}


def test_get_third_highest_mz(peak: Peak):
	assert peak.get_third_highest_mz() == 59
	assert isinstance(peak.get_third_highest_mz(), int)

	with pytest.raises(ValueError, match="Mass spectrum is unset."):
		Peak(test_float).get_third_highest_mz()


def test_ic_mass():
	peak = ICPeak(12.34, 55)
	uid = peak.UID
	assert is_number(peak.ic_mass)
	assert peak.ic_mass == 55
	peak.ic_mass = 12
	assert peak.UID != uid
	assert peak.ic_mass == 12

	peak.ic_mass = 1234
	assert peak.ic_mass == 1234

	# Errors
	for obj in [*test_sequences, test_string, test_dict]:
		with pytest.raises(TypeError):
			peak.ic_mass = obj  # type: ignore[assignment]


def test_mass_spectrum(peak: Peak, im_i: IntensityMatrix):
	scan_i = im_i.get_index_at_time(31.17 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)

	assert isinstance(peak.mass_spectrum, MassSpectrum)
	assert peak.mass_spectrum == ms

	peak = Peak(test_float)
	assert peak.mass_spectrum == MassSpectrum([], [])
	assert not peak.mass_spectrum
	peak.mass_spectrum = ms
	assert peak.mass_spectrum == ms

	peak = Peak(test_float)
	assert peak.mass_spectrum == MassSpectrum([], [])
	assert not peak.mass_spectrum

	peak.mass_spectrum = ms
	assert isinstance(peak.mass_spectrum, MassSpectrum)
	assert isinstance(peak.mass_spectrum.mass_spec, list)

	for obj in [test_string, *test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			peak.mass_spectrum = obj  # type: ignore[assignment]


def test_null_mass(peak: Peak):
	peak = copy.deepcopy(peak)
	uid = peak.UID

	peak.null_mass(73)
	peak.null_mass(147.0)

	index_73 = peak.mass_spectrum.mass_list.index(73)
	assert peak.mass_spectrum.mass_spec[index_73] == 0
	index_147 = peak.mass_spectrum.mass_list.index(147)
	assert peak.mass_spectrum.mass_spec[index_147] == 0

	assert peak.UID != uid

	# Errors
	with pytest.raises(ValueError, match="Mass spectrum is unset."):
		Peak(test_float).null_mass(1)
	for obj in [test_string, *test_lists, test_dict]:
		with pytest.raises(TypeError):
			Peak(test_float, peak.mass_spectrum).null_mass(obj)  # type: ignore[arg-type]
	with pytest.raises(IndexError):
		Peak(test_float, peak.mass_spectrum).null_mass(1)
	with pytest.raises(IndexError):
		Peak(test_float, peak.mass_spectrum).null_mass(10000)


def test_rt(peak: Peak):
	assert isinstance(peak.rt, float)
	assert peak.rt == 12.34


def test_UID(peak: Peak):
	# Get the peak's unique ID
	# Consists of the two most abundant ions and their ratio,
	# and the retention time (in the format set by minutes=True or False)
	assert isinstance(peak.UID, str)
	assert peak.UID == "131-73-42-12.34"

	assert isinstance(Peak(test_float).UID, str)


def test_another_peak(im_i: IntensityMatrix, peak: Peak):
	# A different peak
	scan_i = im_i.get_index_at_time(31.44 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)
	peak2 = Peak(31.44, ms, minutes=True)
	assert peak2.rt == 1886.4
	assert peak2.UID == "207-68-42-1886.40"
	assert peak.UID != peak2.UID


def test_outlier(peak: Peak):
	assert isinstance(peak.is_outlier, bool)
	assert peak.is_outlier is False

	assert Peak(12.34, outlier=True).is_outlier is True


def test_top_ions(peak: Peak):
	with pytest.warns(DeprecationWarning):
		assert isinstance(top_ions_v1(peak, 10), list)
	with pytest.warns(DeprecationWarning):
		assert len(top_ions_v1(peak, 10)) == 10
	with pytest.warns(DeprecationWarning):
		assert len(top_ions_v1(peak)) == 5
	with pytest.warns(DeprecationWarning):
		assert top_ions_v1(peak, 10)[0] == 55

	for obj in [test_string, *test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError), pytest.warns(DeprecationWarning):
			top_ions_v1(obj)

	for obj in [test_string, test_float, test_dict, *test_lists]:
		with pytest.raises(TypeError), pytest.warns(DeprecationWarning):
			top_ions_v1(peak, obj)

	with pytest.warns(DeprecationWarning):
		assert isinstance(top_ions_v2(peak, 10), list)
	with pytest.warns(DeprecationWarning):
		assert len(top_ions_v2(peak, 10)) == 10
	with pytest.warns(DeprecationWarning):
		assert len(top_ions_v2(peak)) == 5
	with pytest.warns(DeprecationWarning):
		assert top_ions_v2(peak, 10)[0] == 55

	for obj in [test_string, *test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError), pytest.warns(DeprecationWarning):
			top_ions_v2(obj)

	for obj in [test_string, test_float, test_dict, *test_lists]:
		with pytest.raises(TypeError), pytest.warns(DeprecationWarning):
			top_ions_v2(peak, obj)

	assert isinstance(peak.top_ions(10), list)
	assert len(peak.top_ions(10)) == 10
	assert len(peak.top_ions()) == 5
	assert peak.top_ions(10)[0] == 55

	for obj in [test_string, test_float, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			peak.top_ions(obj)  # type: ignore[arg-type]


# Inherited Methods from pymsBaseClass


def test_dump(peak: Peak, tmp_pathplus: PathPlus):
	peak.dump(tmp_pathplus / "Peak_dump.dat")

	# Errors
	for obj in [test_list_strs, test_dict, test_list_ints, test_tuple, *test_numbers]:
		with pytest.raises(TypeError):
			peak.dump(obj)  # type: ignore[arg-type]

	# Read and check values
	assert (tmp_pathplus / "Peak_dump.dat").exists()
	loaded_peak = _pickle_load_path(tmp_pathplus / "Peak_dump.dat")
	assert loaded_peak == peak

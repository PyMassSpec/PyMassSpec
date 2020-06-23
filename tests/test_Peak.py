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
import pickle
from numbers import Number

# 3rd party
import deprecation  # type: ignore
import pytest  # type: ignore

# pyms
from pyms.Peak import Peak
from pyms.Peak.Function import peak_sum_area, top_ions_v1, top_ions_v2
from pyms.Spectrum import MassSpectrum

# tests
from .constants import *


def test_Peak(im_i, peak):
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
			Peak(obj, ms, minutes=True)

		with pytest.raises(TypeError):
			Peak(test_float, obj, minutes=False)

	Peak(test_float, test_int)
	Peak(test_float, test_float)


def test_equality(peak):
	assert peak == Peak(peak.rt, peak.mass_spectrum)
	assert peak != Peak(50, peak.mass_spectrum)


@pytest.mark.parametrize(
		"val",
		[test_list_ints, test_list_strs, test_tuple, test_string, test_int, test_float, test_dict])
def test_inequality(peak, val):
	assert peak != val


def test_area(im_i, peak):
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
			Peak(test_float, ms).area = obj
	with pytest.raises(ValueError):
		Peak(test_float, ms).area = -1


def test_bounds(peak):
	peak = copy.copy(peak)

	# Setter
	peak.bounds = (11, 12, 13)

	for obj in [test_string, *test_numbers, test_dict, ["a", "b", "c"], test_tuple]:
		with pytest.raises(TypeError):
			peak.bounds = obj

	for obj in [*test_lists, (1, 2), [1, 2, 3, 4]]:
		with pytest.raises(ValueError):
			peak.bounds = obj

	# Getter
	assert peak.bounds == (11, 12, 13)
	assert isinstance(peak.bounds, tuple)
	peak2 = Peak(test_float)
	peak2.bounds = [11, 12, 13]
	assert peak2.bounds == [11, 12, 13]
	assert isinstance(peak2.bounds, list)

	# set_bounds
	peak3 = Peak(test_float)
	peak3.set_bounds(11, 12, 13)
	assert peak3.bounds == (11, 12, 13)
	assert isinstance(peak3.bounds, tuple)

	for obj in [*test_sequences, test_string, test_dict, test_float]:
		with pytest.raises(TypeError):
			print(obj)
			peak3.set_bounds(obj, 12, 13)
		with pytest.raises(TypeError):
			peak3.set_bounds(11, obj, 13)
		with pytest.raises(TypeError):
			peak3.set_bounds(11, 12, obj)


def test_crop_mass(peak):
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
		with pytest.raises(TypeError):
			peak2.crop_mass(obj, 450)
		with pytest.raises(TypeError):
			peak2.crop_mass(450, obj)

	with pytest.raises(ValueError):
		peak2.crop_mass(100, 0)
	with pytest.raises(ValueError):
		peak2.crop_mass(10, 450)
	with pytest.raises(ValueError):
		peak2.crop_mass(60, 500)
	with pytest.warns(Warning):
		peak2.crop_mass(60, 65)


@deprecation.fail_if_not_removed
def test_get_area(peak):
	with pytest.warns(DeprecationWarning):
		peak.get_area()


@deprecation.fail_if_not_removed
def test_get_ic_mass(peak):
	with pytest.warns(DeprecationWarning):
		peak.get_ic_mass()


def test_get_int_of_ion(peak):
	assert peak.get_int_of_ion(100) == 3888.0
	assert peak.get_int_of_ion(200) == 0.0
	assert isinstance(peak.get_int_of_ion(100), Number)

	with pytest.raises(IndexError):
		peak.get_int_of_ion(1)
	with pytest.raises(IndexError):
		peak.get_int_of_ion(1000000)


def test_ion_area(peak):
	peak = copy.deepcopy(peak)

	assert peak.get_ion_area(1) is None

	peak.set_ion_area(1, 1234)
	peak.set_ion_area(2, 1234.56)

	assert isinstance(peak.get_ion_area(1), Number)
	assert isinstance(peak.get_ion_area(2), Number)
	assert peak.get_ion_area(1) == 1234

	# Errors
	for obj in [test_dict, *test_sequences, test_float, test_string]:
		with pytest.raises(TypeError):
			peak.set_ion_area(obj, test_int)
	for obj in [test_dict, *test_sequences, test_string]:
		with pytest.raises(TypeError):
			peak.set_ion_area(1, obj)


def test_ion_areas(peak):
	peak = copy.deepcopy(peak)

	with pytest.raises(ValueError):
		peak.ion_areas

	peak.ion_areas = {1: 1234, 2: 1234, 3: 1234}

	with pytest.warns(DeprecationWarning):
		peak.set_ion_areas(peak.ion_areas)

	with pytest.warns(DeprecationWarning):
		peak.get_ion_areas()

	for obj in [*test_numbers, test_string, test_list_strs, test_list_ints, tuple]:
		with pytest.raises(TypeError):
			peak.ion_areas = obj

	assert peak.ion_areas == {1: 1234, 2: 1234, 3: 1234}


@deprecation.fail_if_not_removed
def test_get_mass_spectrum(peak):
	with pytest.warns(DeprecationWarning):
		peak.get_mass_spectrum()


@deprecation.fail_if_not_removed
def test_get_pt_bounds(peak):
	with pytest.warns(DeprecationWarning):
		peak.get_pt_bounds()


@deprecation.fail_if_not_removed
def test_get_rt(peak):
	with pytest.warns(DeprecationWarning):
		peak.get_rt()


def test_get_third_highest_mz(peak):
	assert peak.get_third_highest_mz() == 59
	assert isinstance(peak.get_third_highest_mz(), int)

	assert Peak(test_float, test_float).get_third_highest_mz() is None

	# with pytest.raises(AttributeError):
	assert Peak(test_float).get_third_highest_mz() is None


@deprecation.fail_if_not_removed
def test_get_UID(peak):
	with pytest.warns(DeprecationWarning):
		peak.get_UID()


def test_ic_mass():
	peak = Peak(12.34, 55)
	uid = peak.UID
	assert isinstance(peak.ic_mass, Number)
	assert peak.ic_mass == 55
	peak.ic_mass = 12
	assert peak.mass_spectrum is None
	assert peak.UID != uid
	assert peak.ic_mass == 12

	peak.ic_mass = 1234
	assert peak.ic_mass == 1234

	# Errors
	for obj in [*test_sequences, test_string, test_dict]:
		with pytest.raises(TypeError):
			peak.ic_mass = obj


def test_mass_spectrum(peak, im_i):
	scan_i = im_i.get_index_at_time(31.17 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)

	assert isinstance(peak.mass_spectrum, MassSpectrum)
	assert peak.mass_spectrum == ms
	assert peak.ic_mass is None

	peak = Peak(test_float)
	assert peak.mass_spectrum is None
	peak.mass_spectrum = ms
	assert peak.mass_spectrum == ms
	assert peak.ic_mass is None

	peak = Peak(test_float)
	assert peak.mass_spectrum is None
	peak.mass_spectrum = ms
	assert isinstance(peak.mass_spectrum, MassSpectrum)
	assert isinstance(peak.mass_spectrum.mass_spec, list)

	for obj in [test_string, *test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			peak.mass_spectrum = obj


def test_null_mass(peak):
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
	with pytest.raises(NameError):
		Peak(test_float).null_mass(1)
	for obj in [test_string, *test_lists, test_dict]:
		with pytest.raises(TypeError):
			Peak(test_float, peak.mass_spectrum).null_mass(obj)
	with pytest.raises(IndexError):
		Peak(test_float, peak.mass_spectrum).null_mass(1)
	with pytest.raises(IndexError):
		Peak(test_float, peak.mass_spectrum).null_mass(10000)


def test_rt(peak):
	assert isinstance(peak.rt, float)
	assert peak.rt == 12.34


def test_set_area(peak):
	peak = copy.copy(peak)

	with pytest.warns(DeprecationWarning):
		peak.set_area(1234)


def test_set_ic_mass(peak):
	peak = copy.deepcopy(peak)
	with pytest.warns(DeprecationWarning):
		peak.set_ic_mass(12)


def test_set_mass_spectrum(peak, im_i):
	peak = copy.deepcopy(peak)
	with pytest.warns(DeprecationWarning):
		peak.set_mass_spectrum(im_i.get_ms_at_index(123))


def test_set_pt_bounds(peak):
	peak = copy.copy(peak)

	with pytest.warns(DeprecationWarning):
		peak.set_pt_bounds((1, 2, 3))


def test_UID(peak):
	# Get the peak's unique ID
	# Consists of the two most abundant ions and their ratio,
	# and the retention time (in the format set by minutes=True or False)
	assert isinstance(peak.UID, str)
	assert peak.UID == '131-73-42-12.34'

	assert isinstance(Peak(test_float).UID, str)


def test_another_peak(im_i, peak):
	# A different peak
	scan_i = im_i.get_index_at_time(31.44 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)
	peak2 = Peak(31.44, ms, minutes=True)
	assert peak2.rt == 1886.4
	assert peak2.UID == '207-68-42-1886.40'
	assert peak.UID != peak2.UID


def test_outlier(peak):
	assert isinstance(peak.is_outlier, bool)
	assert peak.is_outlier is False

	assert Peak(12.34, outlier=True).is_outlier is True


def test_top_ions(peak):
	with pytest.warns(DeprecationWarning):
		assert isinstance(top_ions_v1(peak, 10), list)
	with pytest.warns(DeprecationWarning):
		assert len(top_ions_v1(peak, 10)) == 10
	with pytest.warns(DeprecationWarning):
		assert len(top_ions_v1(peak)) == 5
	with pytest.warns(DeprecationWarning):
		assert top_ions_v1(peak, 10)[0] == 55

	for obj in [test_string, *test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			with pytest.warns(DeprecationWarning):
				top_ions_v1(obj)

	for obj in [test_string, test_float, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			with pytest.warns(DeprecationWarning):
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
		with pytest.raises(TypeError):
			with pytest.warns(DeprecationWarning):
				top_ions_v2(obj)

	for obj in [test_string, test_float, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			with pytest.warns(DeprecationWarning):
				top_ions_v2(peak, obj)

	assert isinstance(peak.top_ions(10), list)
	assert len(peak.top_ions(10)) == 10
	assert len(peak.top_ions()) == 5
	assert peak.top_ions(10)[0] == 55

	for obj in [test_string, test_float, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			peak.top_ions(obj)


# Inherited Methods from pymsBaseClass

def test_dump(peak, outputdir):
	peak.dump(outputdir / "Peak_dump.dat")

	# Errors
	for obj in [test_list_strs, test_dict, test_list_ints, test_tuple, *test_numbers]:
		with pytest.raises(TypeError):
			peak.dump(obj)

	# Read and check values
	assert (outputdir / "Peak_dump.dat").exists()
	loaded_peak = pickle.load((outputdir / "Peak_dump.dat").open("rb"))
	assert loaded_peak == peak

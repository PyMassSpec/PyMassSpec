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
import pickle
from copy import deepcopy

# 3rd party
import pytest
from pytest_regressions.file_regression import FileRegressionFixture

# this package
from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IonChromatogram import IonChromatogram
from pyms.Spectrum import Scan
from tests.constants import *


@pytest.fixture(scope="module")
def andi(pyms_datadir):
	print("data")
	return ANDI_reader(pyms_datadir / "gc01_0812_066.cdf")


# @pytest.fixture(scope="module")
# def im_andi(data):
# 	# build an intensity matrix object from the data
# 	return build_intensity_matrix(data)
#
#
# @pytest.fixture(scope="module")
# def tic_andi(data):
# 	# get the TIC
# 	return deepcopy(data.tic)
#
#
# @pytest.fixture(scope="module")
# def im_i_andi(data):
# 	# build an intensity matrix object from the data
# 	return build_intensity_matrix_i(data)
#
#
# @pytest.fixture(scope="function")
# def peak_list(im_i):
# 	im_i = deepcopy(im_i)
#
# 	# Intensity matrix size (scans, masses)
# 	n_scan, n_mz = im_i.size
#
# 	# noise filter and baseline correct
# 	for ii in range(n_mz):
# 		ic = im_i.get_ic_at_index(ii)
# 		ic_smooth = savitzky_golay(ic)
# 		ic_bc = tophat(ic_smooth, struct="1.5m")
# 		im_i.set_ic_at_index(ii, ic_bc)
#
# 	# Use Biller and Biemann technique to find apexing ions at a scan
# 	# default is maxima over three scans and not to combine with any neighbouring
# 	# scan.
# 	peak_list = BillerBiemann(im_i, points=9, scans=2)
# 	return peak_list
#
#
# @pytest.fixture(scope="function")
# def filtered_peak_list(im_i, peak_list):
# 	# peak_list = deepcopy(peak_list)
# 	# do peak detection on pre-trimmed data
# 	# trim by relative intensity
# 	apl = rel_threshold(peak_list, 2, copy_peaks=False)
#
# 	# trim by threshold
# 	new_peak_list = num_ions_threshold(apl, 3, 30, copy_peaks=False)
#
# 	# ignore TMS ions and set mass range
# 	for peak in new_peak_list:
# 		peak.crop_mass(50, 400)
# 		peak.null_mass(73)
# 		peak.null_mass(147)
#
# 		# find area
# 		area = peak_sum_area(im_i, peak)
# 		peak.area = area
# 		area_dict = peak_top_ion_areas(im_i, peak)
# 		peak.ion_areas = area_dict
#
# 	return new_peak_list
#
#
# @pytest.fixture(scope="session")
# def peak(im_i):
# 	scan_i = im_i.get_index_at_time(31.17 * 60.0)
# 	ms = im_i.get_ms_at_index(scan_i)
# 	return Peak(12.34, ms)
#
#
# @pytest.fixture(scope="session")
# def ms(im_i):
# 	return deepcopy(im_i.get_ms_at_index(0))
#
#
# @pytest.fixture(scope="session")
# def scan(data):
# 	# return deepcopy(im_i.get_scan_at_index(0))
# 	return deepcopy(data.scan_list[0])
#
#
# @pytest.fixture(scope="function")
# def expr(filtered_peak_list):
# 	# create an experiment
# 	return Experiment("ELEY_1_SUBTRACT", filtered_peak_list)


def test_ANDI_reader(pyms_datadir):
	# Errors
	for obj in [*test_numbers, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			ANDI_reader(obj)  # type: ignore

	with pytest.raises(FileNotFoundError):
		ANDI_reader(test_string)  # type: ignore


# def test_ANDI_OpenChrom_reader(pyms_datadir):
# todo


def test_GCMS_data(andi):
	assert isinstance(andi, GCMS_data)

	GCMS_data(andi.time_list, andi.scan_list)

	# Errors
	for obj in [test_string, *test_numbers, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(obj, andi.scan_list)  # type: ignore

	for obj in [test_string, *test_numbers, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(andi.time_list, obj)  # type: ignore


def test_len(andi):
	assert len(andi) == 9865


def test_equality(andi):
	assert andi == GCMS_data(andi.time_list, andi.scan_list)
	assert andi != GCMS_data(list(range(len(andi.scan_list))), andi.scan_list)
	assert andi != test_string
	assert andi != test_int
	assert andi != test_float
	assert andi != test_list_ints
	assert andi != test_list_strs
	assert andi != test_tuple
	assert andi != test_dict


def test_info(capsys, andi):
	andi.info()
	captured = capsys.readouterr()
	assert captured.out == """ Data retention time range: 5.093 min -- 66.795 min
 Time step: 0.375 s (std=0.000 s)
 Number of scans: 9865
 Minimum m/z measured: 50.000
 Maximum m/z measured: 599.900
 Mean number of m/z values per scan: 56
 Median number of m/z values per scan: 40
"""


def test_scan_list(andi):
	# raw scans
	scans = andi.scan_list

	assert isinstance(scans, list)
	assert isinstance(scans[0], Scan)
	assert len(scans[0]) == 622
	assert isinstance(scans[0].mass_list, list)
	# 1st mass value for 1st scan
	assert isinstance(scans[0].mass_list[0], float)
	assert scans[0].mass_list[0] == 50.099998474121094

	assert isinstance(scans[0].intensity_list, list)
	# 1st intensity value for 1st scan
	assert isinstance(scans[0].intensity_list[0], float)
	assert scans[0].intensity_list[0] == 22128.0

	# minimum mass found in 1st scan
	assert isinstance(scans[0].min_mass, float)
	assert scans[0].min_mass == 50.099998474121094

	# maximum mass found in 1st scan
	assert isinstance(scans[0].max_mass, float)
	assert scans[0].max_mass == 599.4000244140625


def test_tic(andi):
	tic = andi.tic
	assert isinstance(tic, IonChromatogram)
	# number of scans in TIC
	assert len(tic) == 9865
	assert len(tic) == len(andi.time_list)

	# start time of TIC
	assert isinstance(tic.get_time_at_index(0), float)
	assert tic.get_time_at_index(0) == 305.582
	assert isinstance(tic.get_index_at_time(305.6), int)
	assert tic.get_index_at_time(305.6) == 0
	assert tic.get_index_at_time(306) == 1

	assert isinstance(tic.get_intensity_at_index(44), float)
	assert tic.get_intensity_at_index(44) == 21685482.0

	assert isinstance(tic.time_list, list)
	assert tic.time_list[0] == 305.582

	assert isinstance(tic.time_step, float)

	assert isinstance(tic.is_tic(), bool)
	assert tic.is_tic()


def test_trim(andi):
	# time
	trimmed = deepcopy(andi)
	trimmed.trim("6.5m", "21m")

	assert trimmed.min_mass == 50.099998474121094
	assert trimmed.max_mass == 542.0

	time = trimmed.time_list
	assert len(time) == 2319
	assert time[0] == 390.404
	assert trimmed.get_index_at_time(400.0) == 26

	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 50.099998474121094

	# Scans
	trimmed = deepcopy(andi)
	trimmed.trim(10, 2000)

	assert trimmed.min_mass == 50.0
	assert trimmed.max_mass == 599.9000244140625

	time = trimmed.time_list
	assert len(time) == 1992
	assert time[0] == 308.96000000000004
	assert trimmed.get_index_at_time(1000.0) == 1841

	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 50.099998474121094

	trimmed.trim(end=1000)
	assert len(trimmed.time_list) == 1001

	trimmed.trim(begin=2)
	assert len(trimmed.time_list) == 1000

	# Errors
	with pytest.raises(SyntaxError):
		trimmed.trim()

	for obj in [*test_sequences, test_dict]:
		with pytest.raises(TypeError):
			trimmed.trim(begin=obj)
		with pytest.raises(TypeError):
			trimmed.trim(end=obj)


@pytest.mark.parametrize("filename", [
		"andi_gcms_data.I.csv",
		"andi_gcms_data.mz.csv",
		])
def test_write(andi, tmp_pathplus, file_regression: FileRegressionFixture, filename):
	andi.write(tmp_pathplus / "andi_gcms_data")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			andi.write(obj)

	# Read file and check values
	assert (tmp_pathplus / filename).exists()
	file_regression.check((tmp_pathplus / filename).read_text(), encoding="UTF-8", extension=".csv")


def test_write_intensities_stream(andi, tmp_pathplus, file_regression: FileRegressionFixture):
	filename = "andi_intensity_stream.csv"
	andi.write_intensities_stream(tmp_pathplus / filename)

	# Errors
	for obj in [test_list_strs, test_dict, test_list_ints, test_tuple, *test_numbers]:
		with pytest.raises(TypeError):
			andi.write_intensities_stream(obj)

	# Read file and check values
	assert (tmp_pathplus / filename).exists()
	file_regression.check((tmp_pathplus / filename).read_text(), encoding="UTF-8", extension=".csv")


# Inherited Methods from pymsBaseClass


def test_dump(andi, tmp_pathplus):
	andi.dump(tmp_pathplus / "ANDI_dump.dat")

	# Errors
	for obj in [test_list_strs, test_dict, test_list_ints, test_tuple, *test_numbers]:
		with pytest.raises(TypeError):
			andi.dump(obj)

	# Read and check values
	assert (tmp_pathplus / "ANDI_dump.dat").exists()
	loaded_data = pickle.load((tmp_pathplus / "ANDI_dump.dat").open("rb"))
	assert loaded_data == andi
	assert len(loaded_data) == len(andi)


# Inherited Methods from TimeListMixin


def test_time_list(andi):
	time = andi.time_list
	assert isinstance(time, list)
	# number of retention times
	assert len(time) == 9865
	# retention time of 1st scan:
	assert isinstance(time[0], float)
	assert time[0] == 305.582


# Inherited Methods from MaxMinMassMixin


def test_max_mass(andi):
	# maximum mass found in all data
	assert isinstance(andi.max_mass, float)
	assert andi.max_mass == 599.9000244140625


def test_min_mass(andi):
	assert isinstance(andi.min_mass, float)
	# minimum mass found in all data
	assert andi.min_mass == 50.0


# Inherited Methods from GetIndexTimeMixin


def test_get_index_at_time(andi):
	# index of 400sec in time_list
	assert isinstance(andi.get_index_at_time(400.0), int)
	assert andi.get_index_at_time(400.0) == 252

	# Errors
	for obj in [test_dict, *test_lists, test_string, test_tuple]:
		with pytest.raises(TypeError):
			andi.get_index_at_time(obj)
	with pytest.raises(IndexError):
		andi.get_index_at_time(0)
	with pytest.raises(IndexError):
		andi.get_index_at_time(100000)


def test_get_time_at_index(andi):
	assert isinstance(andi.get_time_at_index(400), float)
	assert andi.get_time_at_index(400) == 455.71

	# Errors
	for obj in [test_dict, *test_lists, test_string, test_tuple]:
		with pytest.raises(TypeError):
			andi.get_time_at_index(obj)
	with pytest.raises(IndexError):
		andi.get_time_at_index(-1)
	with pytest.raises(IndexError):
		andi.get_time_at_index(1000000)


# Test GCMS.Function

# def test_diff(data):
# TODO

# def test_ic_window_points(data):
# todo

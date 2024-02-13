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
import os
from copy import deepcopy
from typing import cast

# 3rd party
import numpy
import pytest
from coincidence.regressions import AdvancedDataRegressionFixture, AdvancedFileRegressionFixture, _representer_for
from domdf_python_tools.paths import PathPlus
from pytest_regressions.data_regression import RegressionYamlDumper

# this package
from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IonChromatogram import IonChromatogram
from pyms.Spectrum import Scan
from pyms.Utils.Utils import _pickle_load_path

# this package
from .constants import *


def test_JCAMP_reader_errors():
	# Errors
	for obj in [*test_numbers, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			JCAMP_reader(obj)  # type: ignore[arg-type]

	with pytest.raises(FileNotFoundError):
		JCAMP_reader(test_string)


# def test_JCAMP_OpenChrom_reader(pyms_datadir):
# todo


def test_GCMS_data(data: GCMS_data):
	assert isinstance(data, GCMS_data)

	GCMS_data(data.time_list, data.scan_list)

	# Errors
	for obj in [test_string, *test_numbers, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(obj, data.scan_list)  # type: ignore[arg-type]

	for obj in [test_string, *test_numbers, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(data.time_list, obj)  # type: ignore[arg-type]


def test_len(data: GCMS_data):
	assert len(data) == 2103


def test_equality(data: GCMS_data):
	assert data == GCMS_data(data.time_list, data.scan_list)
	assert data != GCMS_data(list(range(len(data.scan_list))), data.scan_list)
	assert data != test_string
	assert data != test_int
	assert data != test_float
	assert data != test_list_ints
	assert data != test_list_strs
	assert data != test_tuple
	assert data != test_dict


def test_info(capsys, data: GCMS_data):
	data.info()
	captured = capsys.readouterr()
	assert captured.out.splitlines() == [
			" Data retention time range: 0.018 min -- 37.013 min",
			" Time step: 1.056 s (std=0.000 s)",
			" Number of scans: 2103",
			" Minimum m/z measured: 50.252",
			" Maximum m/z measured: 499.623",
			" Mean number of m/z values per scan: 99",
			" Median number of m/z values per scan: 98",
			]


def test_scan_list(data: GCMS_data):
	# raw scans
	scans = data.scan_list

	assert isinstance(scans, list)
	assert isinstance(scans[0], Scan)
	assert len(scans[0]) == 101
	assert isinstance(scans[0].mass_list, list)
	# 1st mass value for 1st scan
	assert isinstance(scans[0].mass_list[0], float)
	assert scans[0].mass_list[0] == 52.0131

	assert isinstance(scans[0].intensity_list, list)
	# 1st intensity value for 1st scan
	assert isinstance(scans[0].intensity_list[0], float)
	assert scans[0].intensity_list[0] == 5356.0

	# minimum mass found in 1st scan
	assert isinstance(scans[0].min_mass, float)
	assert scans[0].min_mass == 52.0131

	# maximum mass found in 1st scan
	assert isinstance(scans[0].max_mass, float)
	assert scans[0].max_mass == 477.6667


def test_tic(data: GCMS_data):
	tic = data.tic
	assert isinstance(tic, IonChromatogram)
	# number of scans in TIC
	assert len(tic) == 2103
	assert len(tic) == len(data.time_list)

	# start time of TIC
	assert isinstance(tic.get_time_at_index(0), float)
	assert tic.get_time_at_index(0) == 1.05200003833
	assert isinstance(tic.get_index_at_time(1.05203833), int)
	assert tic.get_index_at_time(1.05203833) == 0
	assert tic.get_index_at_time(2) == 1

	assert isinstance(tic.get_intensity_at_index(44), float)
	assert tic.get_intensity_at_index(44) == 3113615.0

	assert isinstance(tic.time_list, list)
	assert tic.time_list[0] == 1.05200003833

	assert isinstance(tic.time_step, float)

	assert isinstance(tic.is_tic(), bool)
	assert tic.is_tic()


def test_trim(data: GCMS_data):
	# time
	trimmed = deepcopy(data)
	trimmed.trim("6.5m", "21m")

	assert trimmed.min_mass == 50.2516
	assert trimmed.max_mass == 499.6226

	time = trimmed.time_list
	assert len(time) == 825
	assert time[0] == 390.715999603
	assert trimmed.get_index_at_time(400.0) == 9

	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 51.0066

	# Scans
	trimmed = deepcopy(data)
	trimmed.trim(1000, 2000)

	assert trimmed.min_mass == 50.2516
	assert trimmed.max_mass == 499.6226

	time = trimmed.time_list
	assert len(time) == 1002
	assert time[0] == 1055.99601746
	assert trimmed.get_index_at_time(1500.0) == 420

	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 50.8808

	trimmed.trim(end=1000)
	assert len(trimmed.time_list) == 1001

	trimmed.trim(begin=2)
	assert len(trimmed.time_list) == 1000

	# Errors
	with pytest.raises(SyntaxError):
		trimmed.trim()

	for obj in [*test_sequences, test_dict]:
		with pytest.raises(TypeError):
			trimmed.trim(begin=obj)  # type: ignore[type-var]
		with pytest.raises(TypeError):
			trimmed.trim(end=obj)  # type: ignore[type-var]


@pytest.mark.parametrize("filename", [
		"jcamp_gcms_data.I.csv",
		"jcamp_gcms_data.mz.csv",
		])
def test_write(
		data: GCMS_data,
		tmp_pathplus: PathPlus,
		advanced_file_regression: AdvancedFileRegressionFixture,
		filename: str,
		):
	data.write(tmp_pathplus / "jcamp_gcms_data")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			data.write(obj)  # type: ignore[arg-type]

	# Read file and check values
	assert (tmp_pathplus / filename).exists()
	advanced_file_regression.check_file(tmp_pathplus / filename)


def test_write_intensities_stream(
		data: GCMS_data,
		tmp_pathplus: PathPlus,
		advanced_file_regression: AdvancedFileRegressionFixture,
		):
	filename = "jcamp_intensity_stream.csv"
	data.write_intensities_stream(tmp_pathplus / "jcamp_intensity_stream.csv")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			data.write_intensities_stream(obj)  # type: ignore[arg-type]

	# Read file and check values
	assert (tmp_pathplus / filename).exists()
	advanced_file_regression.check_file(tmp_pathplus / filename)


# Inherited Methods from pymsBaseClass


def test_dump(data: GCMS_data, tmp_pathplus: PathPlus):
	data.dump(tmp_pathplus / "JCAMP_dump.dat")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			data.dump(obj)  # type: ignore[arg-type]

	# Read and check values
	assert (tmp_pathplus / "JCAMP_dump.dat").exists()
	loaded_data = cast(GCMS_data, _pickle_load_path(tmp_pathplus / "JCAMP_dump.dat"))
	assert loaded_data == data
	assert len(loaded_data) == len(data)


# Inherited Methods from TimeListMixin


def test_time_list(data: GCMS_data):
	time = data.time_list
	assert isinstance(time, list)
	# number of retention times
	assert len(time) == 2103
	# retention time of 1st scan:
	assert isinstance(time[0], float)
	assert time[0] == 1.05200003833


# Inherited Methods from MaxMinMassMixin


def test_max_mass(data: GCMS_data):
	# maximum mass found in all data
	assert isinstance(data.max_mass, float)
	assert data.max_mass == 499.6226


def test_min_mass(data: GCMS_data):
	assert isinstance(data.min_mass, float)
	# minimum mass found in all data
	assert data.min_mass == 50.2516


# Inherited Methods from GetIndexTimeMixin


def test_get_index_at_time(data: GCMS_data):
	# index of 400sec in time_list
	assert isinstance(data.get_index_at_time(400.0), int)
	assert data.get_index_at_time(400.0) == 378

	# Errors
	for obj in [test_dict, *test_lists, test_string, test_tuple]:
		with pytest.raises(TypeError):
			data.get_index_at_time(obj)  # type: ignore[arg-type]
	with pytest.raises(IndexError):
		data.get_index_at_time(0)
	with pytest.raises(IndexError):
		data.get_index_at_time(1000000)


def test_get_time_at_index(data: GCMS_data):
	assert isinstance(data.get_time_at_index(400), float)
	assert data.get_time_at_index(400) == 423.45199585

	# Errors
	for obj in [test_dict, *test_lists, test_string, test_tuple]:
		with pytest.raises(TypeError):
			data.get_time_at_index(obj)  # type: ignore[arg-type]
	with pytest.raises(IndexError):
		data.get_time_at_index(-1)
	with pytest.raises(IndexError):
		data.get_time_at_index(1000000)


# Test GCMS.Function

# def test_diff(data):
# TODO

# def test_ic_window_points(data):
# todo


@_representer_for(Scan)
def _represent_scan(dumper: RegressionYamlDumper, scan: Scan):  # noqa: MAN002
	return dumper.represent_data({
			"mass_list": scan.mass_list,
			"intensity_list": scan.intensity_list,
			})


@_representer_for(
		numpy.int64,
		numpy.int32,
		numpy.float64,
		)
def _represent_numpy(dumper: RegressionYamlDumper, data: int):  # noqa: MAN002
	return dumper.represent_data(int(data))


@pytest.mark.parametrize("filename", ["gc01_0812_066.jdx", "aa_4_aq_110124_056.func_0.JDX"])
def test_jcamp_reader_empty_scan(advanced_data_regression: AdvancedDataRegressionFixture, filename: str):
	data: GCMS_data = JCAMP_reader(os.path.join("tests", "data", filename))
	assert data.min_mass is not None
	assert data.max_mass is not None
	advanced_data_regression.check({"time_list": data.time_list, "scan_list": data.scan_list})

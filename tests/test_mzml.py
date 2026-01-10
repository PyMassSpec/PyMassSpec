# stdlib
from pathlib import Path

# 3rd party
import numpy
import pytest
from coincidence.regressions import AdvancedDataRegressionFixture

# this package
from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.IO.MZML import mzML_reader
from pyms.IonChromatogram import IonChromatogram
from pyms.Spectrum import Scan


@pytest.fixture(scope="module")
def mzml_data(pyms_datadir: Path) -> GCMS_data:
	mzml_file = pyms_datadir / "example.mzML"
	data = mzML_reader(mzml_file)
	return data


def test_mzml_reader(
		pyms_datadir: Path,
		advanced_data_regression: AdvancedDataRegressionFixture,
		):

	# read the raw data
	# mzml_file = pyms_datadir / "1min.mzML"
	mzml_file = pyms_datadir / "example.mzML"
	data = mzML_reader(mzml_file)

	# raw data operations
	print("minimum mass found in all data: ", data.min_mass)
	print("maximum mass found in all data: ", data.max_mass)

	# time
	time = data.time_list
	print("number of retention times: ", len(time))
	print("retention time of 1st scan: ", time[0], "sec")
	print("index of 2sec in time_list: ", data.get_index_at_time(2.0))

	# TIC
	tic = data.tic
	print("number of scans in TIC: ", len(tic))
	print("start time of TIC: ", tic.get_time_at_index(0), "sec")

	# raw scans
	scans = data.scan_list

	print("number of masses in 1st scan: ", len(scans[0]))
	print("1st mass value for 1st scan: ", scans[0].mass_list[0])
	print("1st intensity value for 1st scan: ", scans[0].intensity_list[0])

	print("minimum mass found in 1st scan: ", scans[0].min_mass)
	print("maximum mass found in 1st scan: ", scans[0].max_mass)

	advanced_data_regression.check({
			"min_mass": data.min_mass,
			"max_mass": data.max_mass,
			"n_retention_times": len(time),
			"1st_rt": time[0],
			"2sec_index": data.get_index_at_time(2.0),
			"n_scans_tic": len(tic),
			"tic_start_time": tic.get_time_at_index(0),
			"n_masses_1st_scan": len(scans[0]),
			"1st_mass_1st_scan": scans[0].mass_list[0],
			"1st_intensity_1st_scan": scans[0].intensity_list[0],
			"min_mass_1st_scan": scans[0].min_mass,
			"max_mass_1st_scan": scans[0].max_mass,
			})


def test_masses(mzml_data: GCMS_data):
	assert isinstance(mzml_data.min_mass, float)
	# minimum mass found in all data
	assert numpy.isclose(mzml_data.min_mass, 70.048691)
	# "maximum mass found in all data
	assert isinstance(mzml_data.max_mass, float)
	assert numpy.isclose(mzml_data.max_mass, 898.748962)


def test_times(mzml_data: GCMS_data):
	time = mzml_data.time_list
	assert isinstance(time, list)
	# number of retention times
	assert len(time) == 11
	# retention time of 1st scan:
	assert isinstance(time[0], float)
	assert numpy.isclose(time[0], 0.087954)
	# index of 2sec in time_list
	assert isinstance(mzml_data.get_index_at_time(2.0), int)
	assert mzml_data.get_index_at_time(2.0) == 7


def test_tic(mzml_data: GCMS_data):
	tic = mzml_data.tic

	assert isinstance(tic, IonChromatogram)
	# number of scans in TIC
	assert len(tic) == 11
	assert len(tic) == len(mzml_data.time_list)

	# start time of TIC
	assert isinstance(tic.get_time_at_index(0), float)
	assert numpy.isclose(tic.get_time_at_index(0), 0.087954)


def test_scans(mzml_data: GCMS_data):
	# raw scans
	scans = mzml_data.scan_list

	assert isinstance(scans, list)
	assert isinstance(scans[0], Scan)
	assert isinstance(scans[0].mass_list, list)
	# 1st mass value for 1st scan
	assert isinstance(scans[0].mass_list[0], float)
	assert numpy.isclose(scans[0].mass_list[0], 70.065781)

	assert isinstance(scans[0].intensity_list, list)
	# 1st intensity value for 1st scan
	assert isinstance(scans[0].intensity_list[0], float)
	assert numpy.isclose(scans[0].intensity_list[0], 70541.453125)

	# minimum mass found in 1st scan
	assert isinstance(scans[0].min_mass, float)
	assert numpy.isclose(scans[0].min_mass, 70.065781)

	# maximum mass found in 1st scan
	assert isinstance(scans[0].max_mass, float)
	assert numpy.isclose(scans[0].min_mass, 70.065781)

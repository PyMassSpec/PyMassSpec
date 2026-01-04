# stdlib
from pathlib import Path

# 3rd party
import pytest
from coincidence.regressions import AdvancedDataRegressionFixture

# this package
from pyms.GCMS.IO.MZML import mzML_reader


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


# def test_masses(data):
# 	assert isinstance(data.min_mass, float)
# 	# minimum mass found in all data
# 	assert data.min_mass == 50.2516
# 	#"maximum mass found in all data
# 	assert isinstance(data.max_mass, float)
# 	assert data.max_mass == 499.6226

# def test_times(data):
# 	time = data.time_list
# 	assert isinstance(time, list)
# 	#number of retention times
# 	assert len(time) == 2103
# 	#retention time of 1st scan:
# 	assert isinstance(time[0], float)
# 	assert time[0] == 1.05200003833
# 	#index of 400sec in time_list
# 	assert isinstance(data.get_index_at_time(400.0), int)
# 	assert data.get_index_at_time(400.0) == 378

# def test_tic(data):
# 	tic = data.tic
# 	# this package
# 	from pyms.IonChromatogram import IonChromatogram
# 	assert isinstance(tic, IonChromatogram)
# 	#number of scans in TIC
# 	assert len(tic) == 2103
# 	assert len(tic) == len(data.time_list)

# 	#start time of TIC
# 	assert isinstance(tic.get_time_at_index(0), float)
# 	assert tic.get_time_at_index(0) == 1.05200003833

# def test_scans(data):
# 	# raw scans
# 	scans = data.scan_list
# 	# this package
# 	from pyms.Scan import Scan

# 	assert isinstance(scans, list)
# 	assert isinstance(scans[0], Scan)
# 	assert isinstance(scans[0].mass_list, list)
# 	# 1st mass value for 1st scan
# 	assert isinstance(scans[0].mass_list[0], float)
# 	assert scans[0].mass_list[0] == 52.0131

# 	assert isinstance(scans[0].intensity_list, list)
# 	#1st intensity value for 1st scan
# 	assert isinstance(scans[0].intensity_list[0], float)
# 	assert scans[0].intensity_list[0] == 5356.0

# 	#minimum mass found in 1st scan
# 	assert isinstance(scans[0].min_mass, float)
# 	assert scans[0].min_mass == 52.0131

# 	#maximum mass found in 1st scan
# 	assert isinstance(scans[0].max_mass, float)
# 	assert scans[0].min_mass == 477.6667

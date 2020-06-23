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
import csv
import pickle
from copy import deepcopy

# 3rd party
import deprecation  # type: ignore
import numpy  # type: ignore
import pytest  # type: ignore

# pyms
from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IonChromatogram import IonChromatogram
from pyms.Spectrum import Scan

# tests
from .constants import *


def test_JCAMP_reader(datadir):
	# Errors
	for obj in [*test_numbers, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			JCAMP_reader(obj)

	with pytest.raises(FileNotFoundError):
		JCAMP_reader(test_string)


# def test_JCAMP_OpenChrom_reader(datadir):
	# todo


def test_GCMS_data(data):
	assert isinstance(data, GCMS_data)

	GCMS_data(data.time_list, data.scan_list)

	# Errors
	for obj in [test_string, *test_numbers, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(obj, data.scan_list)

	for obj in [test_string, *test_numbers, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(data.time_list, obj)


def test_len(data):
	assert len(data) == 2103


def test_equality(data):
	assert data == GCMS_data(data.time_list, data.scan_list)
	assert data != GCMS_data(list(range(len(data.scan_list))), data.scan_list)
	assert data != test_string
	assert data != test_int
	assert data != test_float
	assert data != test_list_ints
	assert data != test_list_strs
	assert data != test_tuple
	assert data != test_dict


@deprecation.fail_if_not_removed
def test_get_scan_list(data):
	with pytest.warns(DeprecationWarning):
		data.get_scan_list()


@deprecation.fail_if_not_removed
def test_get_tic(data):
	with pytest.warns(DeprecationWarning):
		data.get_tic()


def test_info(capsys, data):
	data.info()
	captured = capsys.readouterr()
	assert captured.out == """ Data retention time range: 0.018 min -- 37.013 min
 Time step: 1.056 s (std=0.000 s)
 Number of scans: 2103
 Minimum m/z measured: 50.252
 Maximum m/z measured: 499.623
 Mean number of m/z values per scan: 99
 Median number of m/z values per scan: 98
"""


def test_scan_list(data):
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


def test_tic(data):
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


def test_trim(data):
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
			trimmed.trim(begin=obj)
		with pytest.raises(TypeError):
			trimmed.trim(end=obj)


def test_write(data, outputdir):
	data.write(outputdir / "jcamp_gcms_data")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			data.write(obj)

	# Read .I.csv and check values
	assert (outputdir / "jcamp_gcms_data.I.csv").exists()
	i_csv = list(csv.reader((outputdir / "jcamp_gcms_data.I.csv").open()))
	assert [float(x) for x in i_csv[5]] == [
			2067.0, 2294.0, 2038.0, 13478.0, 33348.0, 1946.0, 2067.0, 1542.0, 1431.0,
			762.0, 497.0, 826.0, 2346.0, 1741.0, 1996.0, 935.0, 6849.0, 3422.0, 6371.0,
			1159.0, 13.0, 610.0, 7118.0, 2951.0, 2779.0, 4178.0, 5109.0, 2434.0, 673.0,
			673.0, 118576.0, 12575.0, 919.0, 4887.0, 673.0, 749.0, 2630.0, 1378.0,
			819.0, 1571.0, 673.0, 11260.0, 1521.0, 749.0, 8632.0, 643.0, 46844.0,
			2084.0, 1403.0, 1144.0, 1548.0, 673.0, 1499.0, 673.0, 902.0, 1131.0, 1055.0,
			673.0, 673.0, 5139.0, 657.0, 4978.0, 673.0, 749.0, 9107.0, 644.0, 1431.0,
			749.0, 978.0, 84336.0, 2708.0, 5867.0, 673.0, 902.0, 13.0, 1667200.0,
			342528.0, 102784.0, 1152.0, 673.0, 5583.0, 755.0, 749.0, 826.0, 673.0,
			1422.0, 673.0, 673.0, 896.0, 673.0, 90.0, 673.0, 673.0, 673.0, 673.0, 826.0,
			673.0, 673.0, 673.0, 2059.0, 673.0, 749.0, 673.0, 826.0,
			]

	assert [float(x) for x in i_csv[50]] == [
			1646.0, 3660.0, 579.0, 3975.0, 21666.0, 910.0, 749.0, 2434.0, 3471.0,
			661.0, 2776.0, 544.0, 673.0, 1342.0, 2776.0, 2631.0, 3922.0, 677.0, 5416.0,
			9322.0, 909.0, 3368.0, 343.0, 1370.0, 5831.0, 169952.0, 7492.0, 749.0,
			4389.0, 673.0, 902.0, 1499.0, 3724.0, 1422.0, 673.0, 680.0, 832.0, 50240.0,
			5015.0, 749.0, 749.0, 826.0, 673.0, 3888.0, 2220.0, 749.0, 749.0, 673.0,
			1662.0, 749.0, 749.0, 1532.0, 2925.0, 883.0, 1589.0, 749.0, 749.0, 749.0,
			673.0, 7293.0, 2063.0, 749.0, 826.0, 112340.0, 1400.0, 1291.0, 673.0,
			1775808.0, 281392.0, 127524.0, 2260.0, 673.0, 9469.0, 1572.0, 1651.0,
			673.0, 1322.0, 673.0, 749.0, 673.0, 749.0, 902.0, 826.0, 673.0, 749.0,
			1422.0, 673.0, 673.0, 673.0, 673.0, 673.0,
			]

	assert [float(x) for x in i_csv[500]] == [
			7353.0, 2778.0, 963.0, 1884.0, 1365.0, 702.0, 5346.0, 672.0, 5403.0,
			2089.0, 749.0, 3020.0, 4895.0, 370.0, 2420.0, 7680.0, 644.0, 673.0, 826.0,
			1548.0, 12477.0, 26.0, 5308.0, 1297.0, 2770.0, 1468.0, 1521.0, 673.0,
			673.0, 5839.0, 5211.0, 1397.0, 820.0, 3343.0, 661.0, 749.0, 1431.0,
			368160.0, 8394.0, 90.0, 1970.0, 673.0, 2404.0, 755.0, 1404.0, 9503.0,
			673.0, 2284.0, 749.0, 673.0, 36140.0, 2110.0, 673.0, 1547.0, 2077.0,
			749.0, 826.0, 1880.0, 1422.0, 673.0, 826.0, 2250.0, 732.0, 673.0, 673.0,
			4336.0, 279.0, 902.0, 749.0, 749.0, 826.0, 673.0, 1422.0, 673.0, 2036.0,
			673.0, 673.0, 673.0, 978.0, 673.0, 673.0, 1575.0, 673.0, 1499.0, 673.0,
			673.0, 902.0, 673.0, 749.0, 673.0, 2477.0, 673.0,
			]

	# Read .mz.csv and check values
	assert (outputdir / "jcamp_gcms_data.mz.csv").exists()
	i_csv = list(csv.reader((outputdir / "jcamp_gcms_data.mz.csv").open()))
	assert [float(x) for x in i_csv[5]] == [
			51.6986, 53.0197, 54.5925, 55.8507, 56.8573, 57.6751, 61.4497, 63.0854,
			64.1549, 67.4263, 69.7540, 70.5718, 72.0187, 73.0882, 74.7868, 75.6676,
			76.7370, 77.8065, 78.7502, 81.0150, 82.0215, 83.5314, 85.0412, 85.9849,
			86.9286, 87.9351, 89.0675, 90.8290, 91.7727, 92.7793, 95.8619, 96.8055,
			98.0008, 102.9708, 103.5999, 104.4177, 111.5896, 114.2318, 114.9238,
			116.3708, 117.3144, 119.0130, 121.4665, 127.2543, 130.9031, 132.0355,
			132.9163, 134.1116, 135.8731, 139.6477, 143.5482, 145.8759, 146.8824,
			147.8890, 148.8327, 150.2167, 151.2862, 152.9848, 157.0740, 160.9115,
			161.8552, 162.9247, 164.6862, 169.9077, 176.7650, 177.5199, 179.0298,
			180.7913, 187.1453, 191.1086, 192.3669, 193.3734, 195.2607, 196.7077,
			202.6213, 207.2767, 208.2832, 209.2898, 210.4851, 216.2100, 219.2926,
			220.2992, 222.1236, 233.2588, 237.7884, 243.8278, 245.6522, 246.9104,
			264.5254, 265.9724, 282.1404, 291.5140, 296.5469, 312.6520, 313.6586,
			314.8539, 317.5591, 337.8163, 343.6670, 371.0961, 377.8275, 398.3364,
			406.1373, 465.7766,
			]

	assert [float(x) for x in i_csv[50]] == [
			50.8178, 52.0761, 52.7052, 54.9699, 56.7944, 57.6751, 58.5559, 59.2479,
			60.5690, 61.4497, 62.7080, 64.7211, 66.7342, 70.0685, 72.3333, 72.9624,
			74.7239, 75.3530, 76.7370, 80.7633, 81.7070, 84.6638, 85.8591, 86.9286,
			88.9417, 95.7990, 96.7426, 99.9511, 102.5304, 104.2290, 106.5567, 117.3773,
			119.0130, 125.3670, 126.7510, 127.6318, 128.2609, 132.9163, 133.9857,
			135.1811, 138.8928, 142.9820, 145.2468, 147.0083, 148.0148, 150.5942,
			151.7895, 154.8092, 158.8984, 159.6533, 160.9744, 162.2956, 162.9876,
			163.8054, 165.3782, 168.2721, 170.0965, 174.1227, 175.6955, 177.0796,
			179.0298, 181.5462, 183.2448, 191.1715, 192.3669, 193.3105, 198.5950,
			207.2767, 208.2832, 209.2898, 210.2964, 211.1771, 219.2297, 222.6269,
			238.4174, 246.8475, 264.7770, 304.2220, 305.6689, 330.8961, 332.3431,
			335.0482, 344.6107, 361.0933, 372.8575, 397.3927, 417.3983, 435.3907,
			439.6058, 462.0649, 486.6630,
			]

	assert [float(x) for x in i_csv[500]] == [
			50.7549, 51.8244, 52.7052, 54.4037, 55.0958, 55.7878, 56.9202, 57.8009,
			58.6817, 59.8141, 61.6385, 62.6450, 64.4695, 65.2244, 68.0554, 68.7474,
			69.6910, 72.7108, 73.5286, 74.6610, 76.6112, 77.6807, 78.7502, 79.7567,
			83.2797, 85.2929, 86.4253, 87.7464, 89.5079, 90.8919, 92.0243, 93.0309,
			93.6600, 99.3849, 100.0769, 102.0271, 102.9079, 104.9210, 105.9276,
			108.7586, 110.5201, 114.1689, 115.3642, 116.5595, 117.4402, 118.7614,
			122.5989, 123.4167, 125.6186, 127.0656, 134.0486, 135.3698, 137.5088,
			144.4918, 147.8890, 151.2233, 154.8092, 155.8787, 160.5341, 162.9876,
			165.4411, 167.5171, 168.3350, 169.2786, 175.1293, 176.4504, 207.0880,
			208.2203, 211.3029, 215.8325, 224.8917, 234.8945, 235.9010, 248.8606,
			264.5883, 275.4089, 277.8624, 281.2596, 284.5939, 298.4343, 305.4802,
			315.9863, 332.4689, 338.6342, 341.9055, 361.4078, 379.2745, 392.4228,
			406.5148, 448.6649, 466.8461, 485.7822,
			]


def test_write_intensities_stream(data, outputdir):
	data.write_intensities_stream(outputdir / "jcamp_intensity_stream.csv")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			data.write_intensities_stream(obj)

	# Read and check values
	assert (outputdir / "jcamp_intensity_stream.csv").exists()
	intensity_stream = list((outputdir / "jcamp_intensity_stream.csv").open().readlines())
	assert intensity_stream[5] == "1381.0000\n"
	assert intensity_stream[50] == "673.0000\n"
	assert intensity_stream[500] == "673.0000\n"


# Inherited Methods from pymsBaseClass

def test_dump(data, outputdir):
	data.dump(outputdir / "JCAMP_dump.dat")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			data.dump(obj)

	# Read and check values
	assert (outputdir / "JCAMP_dump.dat").exists()
	loaded_data = pickle.load((outputdir / "JCAMP_dump.dat").open("rb"))
	assert loaded_data == data
	assert len(loaded_data) == len(data)


# Inherited Methods from TimeListMixin

def test_time_list(data):
	time = data.time_list
	assert isinstance(time, list)
	# number of retention times
	assert len(time) == 2103
	# retention time of 1st scan:
	assert isinstance(time[0], float)
	assert time[0] == 1.05200003833


@deprecation.fail_if_not_removed
def test_get_time_list(data):
	with pytest.warns(DeprecationWarning):
		data.get_time_list()


# Inherited Methods from MaxMinMassMixin

@deprecation.fail_if_not_removed
def test_get_max_mass(data):
	with pytest.warns(DeprecationWarning):
		data.get_max_mass()


@deprecation.fail_if_not_removed
def test_get_min_mass(data):
	with pytest.warns(DeprecationWarning):
		data.get_min_mass()


def test_max_mass(data):
	# maximum mass found in all data
	assert isinstance(data.max_mass, float)
	assert data.max_mass == 499.6226


def test_min_mass(data):
	assert isinstance(data.min_mass, float)
	# minimum mass found in all data
	assert data.min_mass == 50.2516


# Inherited Methods from GetIndexTimeMixin

def test_get_index_at_time(data):
	# index of 400sec in time_list
	assert isinstance(data.get_index_at_time(400.0), int)
	assert data.get_index_at_time(400.0) == 378

	# Errors
	for obj in [test_dict, *test_lists, test_string, test_tuple]:
		with pytest.raises(TypeError):
			data.get_index_at_time(obj)
	with pytest.raises(IndexError):
		data.get_index_at_time(0)
	with pytest.raises(IndexError):
		data.get_index_at_time(1000000)


def test_get_time_at_index(data):
	assert isinstance(data.get_time_at_index(400), float)
	assert data.get_time_at_index(400) == 423.45199585

	# Errors
	for obj in [test_dict, *test_lists, test_string, test_tuple]:
		with pytest.raises(TypeError):
			data.get_time_at_index(obj)
	with pytest.raises(IndexError):
		data.get_time_at_index(-1)
	with pytest.raises(IndexError):
		data.get_time_at_index(1000000)

# Test GCMS.Function

# def test_diff(data):
# TODO


# def test_ic_window_points(data):
# todo

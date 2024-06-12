"""
Functions for reading ANDI-MS data files.
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
import os
import pathlib

# 3rd party
from domdf_python_tools.typing import PathLike
from netCDF4 import Dataset  # type: ignore[import]

try:
	# 3rd party
	from mpi4py import MPI  # type: ignore[import]  # noqa: F401
except ModuleNotFoundError:
	pass

# this package
from pyms.GCMS.Class import GCMS_data
from pyms.Spectrum import Scan

__all__ = ["ANDI_reader"]

# netCDF dimension names
__POINT_NUMBER = "point_number"
__SCAN_NUMBER = "scan_number"

# the keys used to create and retrieve certain data from the NetCDF file
__MASS_STRING = "mass_values"
__INTENSITY_STRING = "intensity_values"
__TIME_STRING = "scan_acquisition_time"
__POINT_COUNT = "point_count"


def ANDI_reader(file_name: PathLike) -> GCMS_data:
	"""
	A reader for ANDI-MS NetCDF files.

	:param file_name: The path of the ANDI-MS file

	:return: GC-MS data object

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic, Dominic Davis-Foster
	"""

	if not isinstance(file_name, (str, pathlib.Path)):
		raise TypeError("'file_name' must be a string or a pathlib.Path object")

	if not os.path.isfile(file_name):
		# netCDF4 1.6.0 has stopped raising FileNotFoundError
		# and instead creates an empty file, for some reason.
		raise FileNotFoundError(2, "No such file or directory", file_name)

	rootgrp = Dataset(file_name, "r+", format="NETCDF3_CLASSIC")
	# TODO: find out if netCDF4 throws specific errors that we can use here

	print(f" -> Reading netCDF file '{file_name}'")

	scan_list = []
	mass = rootgrp.variables[__MASS_STRING][:]
	intensity = rootgrp.variables[__INTENSITY_STRING][:]

	scan_lengths = rootgrp.variables["point_count"]  # The number of data points in each scan

	mass_values = mass.tolist()
	intensity_values = intensity.tolist()

	if len(mass_values) != len(intensity_values):
		raise ValueError("The lengths of the mass and intensity lists differ!")

	offset = 0
	for idx, length in enumerate(scan_lengths):
		mass_list = mass_values[offset:offset + length]
		assert len(mass_values[offset:offset + length]) == length
		intensity_list = intensity_values[offset:offset + length]
		assert len(intensity_values[offset:offset + length]) == length
		scan_list.append(Scan(mass_list, intensity_list))
		offset += length

	assert offset == len(mass_values)

	time = rootgrp.variables[__TIME_STRING][:]
	time_list = time.tolist()

	# sanity check
	if len(time_list) != len(scan_list):
		raise ValueError("number of time points does not equal the number of scans")

	rootgrp.close()
	return GCMS_data(time_list, scan_list)


#
# def ANDI_writer(file_name: str, im: IntensityMatrix):
# 	"""
# 	A writer for ANDI-MS NetCDF files
#
# 	:param file_name: The name of the ANDI-MS file
# 	:param im: The IntensityMatrix
#
# 	:author: Andrew Isaac
#
# 	.. TODO:: finish this
# 	"""
#
# 	# netCDF header info for compatability
# 	# attributes
# 	# dataset_completeness   0 CHAR     6 C1+C2
# 	# dataset_origin         4 CHAR    16 Santa Clara, CA
# 	# experiment_date_time_stamp   7 CHAR    20 20081218044500+1100
# 	# experiment_title       6 CHAR     7 mix ma
# 	# experiment_type       10 CHAR    25 Centroided Mass Spectrum
# 	# external_file_ref_0    9 CHAR     8 MA_5C.M
# 	# languages              3 CHAR     8 English
# 	# ms_template_revision   1 CHAR     6 1.0.1
# 	# netcdf_file_date_time_stamp   5 CHAR    20 20090114001531+1100
# 	# netcdf_revision        2 CHAR     6 2.3.2
# 	# number_of_times_calibrated  12 INT      1 0
# 	# number_of_times_processed  11 INT      1 1
# 	# operator_name          8 CHAR    12 Dave and Su
# 	# raw_data_intensity_format  25 CHAR     6 Float
# 	# raw_data_mass_format  23 CHAR     6 Float
# 	# raw_data_time_format  24 CHAR     6 Short
# 	# sample_state          13 CHAR    12 Other State
# 	# test_detector_type    18 CHAR    20 Electron Multiplier
# 	# test_ionization_mode  16 CHAR    16 Electron Impact
# 	# test_ionization_polarity  17 CHAR    18 Positive Polarity
# 	# test_ms_inlet         15 CHAR    17 Capillary Direct
# 	# test_resolution_type  19 CHAR    20 Constant Resolution
# 	# test_scan_direction   21 CHAR     3 Up
# 	# test_scan_function    20 CHAR    10 Mass Scan
# 	# test_scan_law         22 CHAR     7 Linear
# 	# test_separation_type  14 CHAR    18 No Chromatography
#
# 	# dimensions
# 	# _128_byte_string       6    128
# 	# _16_byte_string        3     16
# 	# _255_byte_string       7    255
# 	# _2_byte_string         0      2
# 	# _32_byte_string        4     32
# 	# _4_byte_string         1      4
# 	# _64_byte_string        5     64
# 	# _8_byte_string         2      8
# 	# error_number          10      1
# 	# instrument_number     12      1
# 	# point_number           9 554826   X
# 	# range                  8      2
# 	# scan_number           11   9865
#
# 	# variables
# 	# a_d_coaddition_factor   2 SHORT      0 scan_number(9865)
# 	# a_d_sampling_rate      1 DOUBLE     0 scan_number(9865)
# 	# actual_scan_number     7 INT        0 scan_number(9865)
# 	# error_log              0 CHAR       0 error_number(1), _64_byte_string(64)
# 	# flag_count            15 INT        0 scan_number(9865)
# 	# instrument_app_version  27 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_comments   28 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_fw_version  25 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_id         20 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_mfr        21 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_model      22 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_name       19 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_os_version  26 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_serial_no  23 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# instrument_sw_version  24 CHAR       0 instrument_number(1),
# 	# _32_byte_string(32)
# 	# intensity_values      18 FLOAT      3 point_number(554826)
# 	# inter_scan_time        5 DOUBLE     0 scan_number(9865)
# 	# mass_range_max        10 DOUBLE     0 scan_number(9865)
# 	# mass_range_min         9 DOUBLE     0 scan_number(9865)
# 	# mass_values           16 FLOAT      2 point_number(554826)
# 	# point_count           14 INT        0 scan_number(9865)
# 	# resolution             6 DOUBLE     0 scan_number(9865)
# 	# scan_acquisition_time   3 DOUBLE     0 scan_number(9865)
# 	# scan_duration          4 DOUBLE     0 scan_number(9865)
# 	# scan_index            13 INT        0 scan_number(9865)
# 	# time_range_max        12 DOUBLE     0 scan_number(9865)
# 	# time_range_min        11 DOUBLE     0 scan_number(9865)
# 	# time_values           17 FLOAT      2 point_number(554826)
# 	# total_intensity        8 DOUBLE     1 scan_number(9865)
#
# 	# variable information
# 	# intensity_values attributes
#
# 	# name                 idx type   len value
# 	# -------------------- --- ----   --- -----
# 	# add_offset             1 DOUBLE   1 0.0
# 	# scale_factor           2 DOUBLE   1 1.0
# 	# units                  0 CHAR    26 Arbitrary Intensity Units
#
# 	# mass_values attributes
#
# 	# name                 idx type   len value
# 	# -------------------- --- ----   --- -----
# 	# scale_factor           1 DOUBLE   1 1.0
# 	# units                  0 CHAR     4 M/Z
#
# 	# time_values attributes
#
# 	# name                 idx type   len value
# 	# -------------------- --- ----   --- -----
# 	# scale_factor           1 DOUBLE   1 1.0
# 	# units                  0 CHAR     8 Seconds
#
# 	# total_intensity attributes
#
# 	# name                 idx type   len value
# 	# -------------------- --- ----   --- -----
# 	# units                  0 CHAR    26 Arbitrary Intensity Units
#
# 	if not isinstance(file_name, str):
# 		raise TypeError("'file_name' must be a string")
# 	try:
# 		# Open netCDF file in overwrite mode, creating it if inexistent.
# 		nc = CDF(file_name, NC.WRITE | NC.TRUNC | NC.CREATE)
# 		# Automatically set define and data modes.
# 		nc.automode()
# 	except CDFError:
# 		raise IOError(f"Cannot create file '{file_name}'")
#
# 	mass_list = im.mass_list
# 	time_list = im.time_list
#
# 	# direct access, don't modify
# 	intensity_matrix = im.intensity_array
#
# 	# compress by ignoring zero intensities
# 	# included for consistency with imported netCDF format
# 	mass_values = []
# 	intensity_values = []
# 	point_count_values = []
# 	for row in range(len(intensity_matrix)):
# 		pc = 0  # point count
# 		for col in range(len(intensity_matrix[0])):  # all rows same len
# 			if (intensity_matrix[row][col] > 0):
# 				mass_values.append(mass_list[col])
# 				intensity_values.append(intensity_matrix[row][col])
# 				pc += 1
# 		point_count_values.append(pc)
#
# 	# sanity checks
# 	if not len(time_list) == len(point_count_values):
# 		raise ValueError("number of time points does not equal the number of scans")
#
# 	# create dimensions
# 	# total number of data points
# 	dim_point_number = nc.def_dim(__POINT_NUMBER, len(mass_values))
# 	# number of scans
# 	dim_scan_number = nc.def_dim(__SCAN_NUMBER, len(point_count_values))
#
# 	# create variables
# 	# points
# 	var_mass_values = nc.def_var(__MASS_STRING, NC.FLOAT, dim_point_number)
# 	var_intensity_values = nc.def_var(__INTENSITY_STRING, NC.FLOAT, dim_point_number)
# 	# scans
# 	var_time_list = nc.def_var(__TIME_STRING, NC.DOUBLE, dim_scan_number)
# 	var_point_count_values = nc.def_var(__POINT_COUNT, NC.INT, dim_scan_number)
#
# 	# populate variables
# 	# points
# 	var_mass_values[:] = mass_values
# 	var_intensity_values[:] = intensity_values
# 	# scans
# 	var_time_list[:] = time_list
# 	var_point_count_values[:] = point_count_values
#
# 	# close file
# 	nc.close()

# 3rd party
import pytest

# this package
from pyms.DPA.Alignment import Alignment

# this package
from .constants import *


def test_Alignment(filtered_peak_list):
	#Alignment(Experiment("ELEY_1_SUBTRACT", filtered_peak_list))
	#Alignment(None)

	for obj in [test_string, *test_numbers, *test_lists, test_dict]:
		with pytest.raises(TypeError):
			Alignment(obj)


# def test_experiment():
# 	# Experiment code
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# Peak List
# 	# With list of peak objects
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# Store expr
# 	# With Experiment
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# def test_GCMS():
# 	# Time list
# 	# With list of numbers
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# scan list
# 	# With list of scans
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	#__set_time
# 	for ii in range(len(time_list) - 1):
# 		t1 = time_list[ii]
# 		t2 = time_list[ii + 1]
# 		if not t2 > t1:
# 			error("problem with retention times detected")
#
# 	# get_index_at_time
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
# 	# time out of bunds (time < self.__min_rt) or (time > self.__max_rt):
#
# 	# trim
# 	# With defaults
# 	# With values
# 	# With only End
# 	# With only Begin
# 	# With invalid begin
# 	# With invalid end
# 	# With end < beginning
# 	# With first_scan < 0:
# 	# With last_scan > N-1
#
# def test_Scan():
# 	# Mass List
# 	# With list of numbers
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# intensity list
# 	# With list of numbers
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
#
#
# def test_build_intensity_matrix():
# 	# data
# 	# With GCMS_data
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# Bin interval
# 	# With bin_interval == 0
# 	# with bin_interval > 0:
#
# 	# bin_left
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# bin_right
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
#
# def test_build_intensity_matrix_i():
# 	# data
# 	# With GCMS_data
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# bin_left
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# bin_right
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# # TODO: __fill_bins
#
#
#
# def test_ic_window_points():
# 	# window_sele
# 	# With string
# 	# With int
# 	# With float
# 	# With list
# 	# With dict
#
# 	# With half_window and window_sele % 2 == 0
#
#

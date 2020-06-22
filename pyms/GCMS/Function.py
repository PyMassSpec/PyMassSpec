"""
Provides conversion and information functions for GC-MS data objects
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
import math
import sys

# this package
from typing import Union

from pyms import IonChromatogram
from pyms.GCMS.Class import GCMS_data
from pyms.Utils.Math import rmsd
from pyms.Utils.Time import time_str_secs


def diff(data1: GCMS_data, data2: GCMS_data):
	"""
	Compares two GCMS_data objects

	:param data1: GCMS data set 1
	:type data1: pyms.GCMS.Class.GCMS_data
	:param data2: GCMS data set 2
	:type data2: pyms.GCMS.Class.GCMS_data

	:author: Qiao Wang
	:author: Andrew Isaac
	:author: Vladimir Likic
	"""

	# get time attributes
	time_list1 = data1.get_time_list()
	time_list2 = data2.get_time_list()

	# First, check if two data sets have the same number of retention times.
	if len(time_list1) != len(time_list2):
		print(" The number of retention time points differ.")
		print(f"	First data set: {len(time_list1):d} time points")
		print(f"	Second data set: {len(time_list2):d} time points")
		print(" Data sets are different.")
		return
	else:
		time_rmsd = rmsd(time_list1, time_list2)
		print(" Data sets have the same number of time points.")
		print(f"   Time RMSD: {time_rmsd:.2e}")

	# Second, check if each scan has the same number of m/z intensities
	print(" Checking for consistency in scan lengths ...", end='')
	sys.stdout.flush()

	scan_list1 = data1.get_scan_list()
	scan_list2 = data2.get_scan_list()

	if not len(scan_list1) == len(scan_list2):
		# since the number of rention times are the same, this indicated
		# some unexpected problem with data
		raise ValueError("inconsistency in data detected")

	for ii in range(len(scan_list1)):
		scan1 = scan_list1[ii]
		scan2 = scan_list2[ii]
		mass_list1 = scan1.get_mass_list()
		mass_list2 = scan2.get_mass_list()
		if len(mass_list1) != len(mass_list2):
			print(f"\n Different number of points detected in scan no. {ii:d}")
			print(" Data sets are different.")
			return

	print("OK")

	# Third, if here, calculate the max RMSD for m/z and intensities
	print(" Calculating maximum RMSD for m/z values and intensities ...", end='')
	sys.stdout.flush()

	max_mass_rmsd = 0.0
	max_intensity_rmsd = 0.0

	for ii in range(len(scan_list1)):
		scan1 = scan_list1[ii]
		scan2 = scan_list2[ii]
		mass_list1 = scan1.get_mass_list()
		mass_list2 = scan2.get_mass_list()
		intensity_list1 = scan1.get_intensity_list()
		intensity_list2 = scan2.get_intensity_list()
		mass_rmsd = rmsd(mass_list1, mass_list2)
		if mass_rmsd > max_mass_rmsd:
			max_mass_rmsd = mass_rmsd
		intensity_rmsd = rmsd(intensity_list1, intensity_list2)
		if intensity_rmsd > max_intensity_rmsd:
			max_intensity_rmsd = intensity_rmsd

	print(f"\n   Max m/z RMSD: {max_mass_rmsd:.2e}")
	print(f"   Max intensity RMSD: {max_intensity_rmsd:.2e}")


def ic_window_points(ic: IonChromatogram, window_sele: Union[int, str], half_window: bool = False):
	"""
	Converts the window selection parameter into points based on the
	time step in an ion chromatogram.

	:param ic: ion chromatogram object relevant for the conversion
	:type ic: pyms.IO.Class.IonChromatogram
	:param window_sele: The window selection parameter. This can be an
		integer or time string. If integer, taken as the number of points.
		If a string, must of the form "<NUMBER>s" or "<NUMBER>m",
		specifying a time in seconds or minutes, respectively
	:type window_sele: int or str
	:param half_window: Specifies whether to return half-window
	:type half_window: bool, optional

	:author: Vladimir Likic
	"""

	if not isinstance(window_sele, (int, str)):
		raise TypeError("'window_sele' must be either an integer or a string")

	if isinstance(window_sele, int):
		if half_window:
			if window_sele % 2 == 0:
				raise ValueError("window must be an odd number of points")
			else:
				points = int(math.floor(window_sele * 0.5))
		else:
			points = window_sele
	else:
		time = time_str_secs(window_sele)
		time_step = ic.time_step

		if half_window:
			time = time * 0.5

		points = int(math.floor(time / time_step))

	if half_window:
		if points < 1:
			raise ValueError(f"window too small (half window={points:d})")
	else:
		if points < 2:
			raise ValueError(f"window too small (window={points})")

	return points

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
from typing import Any, Type

# 3rd party
import pytest

# this package
from pyms.IntensityMatrix import IntensityMatrix
from pyms.Spectrum import Scan

# this package
from .constants import *


def test_scan(scan: Scan):
	assert isinstance(scan, Scan)
	assert isinstance(scan.mass_list, list)
	assert isinstance(scan.intensity_list, list)


@pytest.mark.parametrize(
		"obj, expects",
		[
				(test_list_ints, ValueError),
				(test_string, ValueError),
				(test_list_strs, ValueError),
				(test_int, TypeError),
				(test_dict, TypeError),
				],
		)
def test_errors(scan: Scan, obj: Any, expects: Type[Exception]):
	with pytest.raises(expects):
		Scan(obj, scan.intensity_list)

	with pytest.raises(expects):
		Scan(scan.mass_list, obj)


def test_len(scan: Scan):
	assert len(scan) == 101


def test_equality(im: IntensityMatrix, scan: Scan):
	assert scan == Scan(scan.mass_list, scan.intensity_list)
	assert scan != im.get_scan_at_index(1234)


@pytest.mark.parametrize("val", [test_list_ints, test_list_strs, test_tuple, test_string, test_int, test_float])
def test_inequality(scan: Scan, val: Any):
	assert scan != val


@pytest.mark.parametrize(
		"index, mass, intensity",
		[(5, 60.9465, 1381.0), (50, 138.8299, 673.0), (100, 477.6667, 1728.0)],
		)
def test_scan_values(scan: Scan, index: int, mass: float, intensity: float):
	assert scan.mass_list[index] == mass
	assert scan.intensity_list[index] == intensity


#
# def test_zero_length():
# 	# TODO: finish
# 	scan = Scan([], [])

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
from copy import copy, deepcopy
from timeit import timeit

# this package
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Peak.Class import Peak

data = JCAMP_reader(os.path.join("data", "ELEY_1_SUBTRACT.JDX"))
im_i = build_intensity_matrix_i(data)
scan_i = im_i.get_index_at_time(31.17 * 60.0)
ms = im_i.get_ms_at_index(scan_i)
peak = Peak(12.34, ms)


def copy_peak():
	return copy(peak)


def deepcopy_peak():
	return deepcopy(peak)


print(timeit(copy_peak))
print(timeit(deepcopy_peak))


def copy_ms():
	return copy(ms)


def deepcopy_ms():
	return deepcopy(ms)


print(timeit(copy_ms))
print(timeit(deepcopy_ms))

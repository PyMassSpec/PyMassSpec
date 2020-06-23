#############################################################################
#                                                                           #
#    PyMassSpec software for processing of mass-spectrometry data           #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                           s#
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

# 3rd party
import pytest  # type: ignore

# pyms
from pyms.Noise.Analysis import window_analyzer

# tests
from .constants import *


def test_window_anlyzer(tic):
	noise_estimate = window_analyzer(tic, rand_seed=test_int)
	assert noise_estimate == 22524.833209785025

	assert isinstance(noise_estimate, float)
	assert isinstance(window_analyzer(tic), float)
	assert isinstance(window_analyzer(tic, rand_seed=test_string), float)
	assert isinstance(window_analyzer(tic, rand_seed=test_float), float)

	for obj in [test_string, *test_numbers, *test_lists, test_dict]:
		with pytest.raises(TypeError):
			window_analyzer(obj)

	for obj in [*test_lists, test_dict]:
		with pytest.raises(TypeError):
			window_analyzer(tic, rand_seed=obj)
	for obj in [test_float, *test_lists, test_dict]:
		with pytest.raises(TypeError):
			window_analyzer(tic, window=obj)
	for obj in [test_string, test_float, *test_lists, test_dict]:
		with pytest.raises(TypeError):
			window_analyzer(tic, n_windows=obj)

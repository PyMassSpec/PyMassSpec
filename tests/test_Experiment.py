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
import copy

# 3rd party
import pytest  # type: ignore
import deprecation  # type: ignore

# pyms
from pyms.Experiment import Experiment, load_expr, read_expr_list, store_expr
from pyms.Peak.Class import Peak
from pyms.Utils.Utils import is_sequence_of

# tests
from .constants import *


def test_Experiment(expr, filtered_peak_list):
	assert isinstance(expr, Experiment)
	# Errors
	for obj in [*test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			Experiment(obj, filtered_peak_list)

	for obj in [test_string, test_int, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			Experiment(test_string, obj)


def test_equality(expr):
	assert expr == Experiment(expr.expr_code, expr.peak_list)
	assert expr != Experiment(test_string, expr.peak_list)
	assert expr != test_float
	assert expr != test_int
	assert expr != test_string
	assert expr != test_list_ints
	assert expr != test_list_strs
	assert expr != test_dict
	assert expr != test_tuple


def test_len(expr):
	assert len(expr) == 641


@deprecation.fail_if_not_removed
def test_get_expr_code(expr):
	with pytest.warns(DeprecationWarning):
		expr.get_expr_code()


def test_expr_code(expr):
	assert isinstance(expr.expr_code, str)
	assert expr.expr_code == "ELEY_1_SUBTRACT"


@deprecation.fail_if_not_removed
def test_get_peak_list(expr):
	with pytest.warns(DeprecationWarning):
		expr.get_peak_list()


def test_peak_list(expr, filtered_peak_list):
	assert isinstance(expr.peak_list, list)
	assert is_sequence_of(filtered_peak_list, Peak)
	assert expr.peak_list == filtered_peak_list


def test_sele_rt_range(expr, filtered_peak_list):
	expr = copy.copy(expr)

	expr.sele_rt_range(["6.5m", "21m"])
	assert expr.peak_list != filtered_peak_list

	for obj in [test_string, *test_numbers, test_dict]:
		with pytest.raises(TypeError):
			expr.sele_rt_range(obj)

	for obj in [*test_lists]:
		with pytest.raises(ValueError):
			expr.sele_rt_range(obj)


@deprecation.fail_if_not_removed
def test_store_expr(expr, outputdir):
	with pytest.warns(DeprecationWarning):
		store_expr(str(outputdir / "ELEY_1_SUBTRACT_DEPRECATION.expr"), expr)

	for obj in [*test_numbers, test_dict, *test_lists]:
		with pytest.warns(DeprecationWarning):
			with pytest.raises(TypeError):
				store_expr(obj, expr)

	for obj in [*test_numbers, test_string, test_dict, *test_lists]:
		with pytest.warns(DeprecationWarning):
			with pytest.raises(TypeError):
				store_expr(test_string, obj)


@pytest.fixture(scope="function")
def expr_filename(expr, outputdir):
	filename = outputdir / "ELEY_1_SUBTRACT.expr"
	expr.store(filename)
	return filename


def test_store_errors(expr):

	for obj in [*test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			expr.store(obj)


def test_load_expr(filtered_peak_list, datadir, expr_filename):
	expr = load_expr(expr_filename)
	assert isinstance(expr, Experiment)

	assert isinstance(expr.expr_code, str)
	assert expr.expr_code == "ELEY_1_SUBTRACT"

	assert isinstance(expr.peak_list, list)
	assert is_sequence_of(expr.peak_list, Peak)

	assert expr.peak_list == filtered_peak_list
	expr.sele_rt_range(["6.5m", "21m"])

	# Errors
	for obj in [*test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			load_expr(obj)

	with pytest.raises(IOError):
		load_expr(datadir / "non-existent.expr")
	with pytest.raises(IOError):
		load_expr(datadir / "not-an-experiment.expr")


def test_read_expr_list(filtered_peak_list, datadir, expr_filename):
	expr_list = read_expr_list(datadir / "read_expr_list.txt")
	assert isinstance(expr_list, list)
	assert is_sequence_of(expr_list, Experiment)

	expr = expr_list[0]
	assert isinstance(expr.expr_code, str)
	assert expr.expr_code == "ELEY_1_SUBTRACT"

	assert isinstance(expr.peak_list, list)
	assert is_sequence_of(expr.peak_list, Peak)
	assert expr.peak_list == filtered_peak_list

	expr.sele_rt_range(["6.5m", "21m"])

	# Errors
	for obj in [*test_numbers, test_dict, *test_lists]:
		with pytest.raises(TypeError):
			read_expr_list(obj)

	with pytest.raises(IOError):
		read_expr_list("non-existent.expr")
	with pytest.raises((IOError, UnicodeDecodeError)):
		read_expr_list("not-an-experiment.expr")
	with pytest.raises(IOError):
		read_expr_list("__init__.py")

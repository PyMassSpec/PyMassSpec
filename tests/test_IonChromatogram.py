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
from typing import Any, cast

# 3rd party
import numpy
import pytest
from domdf_python_tools.paths import PathPlus

# this package
from pyms.IntensityMatrix import IntensityMatrix
from pyms.IonChromatogram import IonChromatogram
from pyms.Utils.Utils import _pickle_load_path, is_number

# this package
from .constants import *


class TestIonChromatogram:

	def test_success(self, im: IntensityMatrix, tic: IonChromatogram):
		# get the first ion chromatogram of the IntensityMatrix
		ic = im.get_ic_at_index(0)
		assert isinstance(ic, IonChromatogram)
		assert not ic.is_tic()

		# get the ion chromatogram for m/z = 73
		ic = im.get_ic_at_mass(73)
		assert isinstance(ic, IonChromatogram)
		assert not ic.is_tic()

		assert isinstance(tic, IonChromatogram)
		assert tic.is_tic()

		with pytest.raises(ValueError, match="'intensity_list' and 'time_list' differ in length"):
			IonChromatogram(tic.intensity_array, test_list_ints)

	@pytest.mark.parametrize("value", [test_string, *test_numbers, test_list_strs, test_tuple, test_dict])
	def test_errors_intensity_list(self, value: Any, tic: IonChromatogram):
		with pytest.raises(TypeError):
			IonChromatogram(value, tic.time_list)

	@pytest.mark.parametrize("value", [test_string, *test_numbers, test_list_strs, test_dict])
	def test_errors_time_list(self, value: Any, tic: IonChromatogram):
		with pytest.raises(TypeError):
			IonChromatogram(tic.intensity_array, value)

	@pytest.mark.parametrize("value", [test_string, *test_sequences, test_dict])
	def test_errors_mass(self, value: Any, tic: IonChromatogram):
		with pytest.raises(TypeError):
			IonChromatogram(tic.intensity_array, tic.time_list, mass=value)


def test_len(tic: IonChromatogram):
	assert len(tic) == 2103


def test_subtract_ic(im: IntensityMatrix):
	ic1 = im.get_ic_at_index(0)
	assert isinstance(ic1, IonChromatogram)

	ic2 = im.get_ic_at_index(1)

	ic3 = ic1 - ic2
	assert isinstance(ic3, IonChromatogram)


def test_equality(tic: IonChromatogram, im: IntensityMatrix):
	assert tic == IonChromatogram(tic.intensity_array, tic.time_list)
	assert tic != im.get_ic_at_index(0)
	assert tic != test_string
	assert tic != test_int
	assert tic != test_float
	assert tic != test_list_ints
	assert tic != test_list_strs
	assert tic != test_dict
	assert tic != test_tuple


def test_get_intensity_at_index(tic: IonChromatogram):
	assert isinstance(tic.get_intensity_at_index(test_int), float)
	assert tic.get_intensity_at_index(test_int) == 421170.0

	# Errors
	for obj in [test_string, test_float, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			tic.get_intensity_at_index(obj)  # type: ignore[arg-type]

	with pytest.raises(IndexError):
		tic.get_intensity_at_index(-1)
	with pytest.raises(IndexError):
		tic.get_intensity_at_index(10000000)


def test_mass(tic: IonChromatogram, im: IntensityMatrix):
	with pytest.warns(Warning):
		tic.mass

	ic = im.get_ic_at_index(0)
	assert is_number(ic.mass)
	assert ic.mass == 50.2516


def test_intensity_array(tic: IonChromatogram, im: IntensityMatrix):
	tic = copy.deepcopy(tic)

	assert isinstance(tic.intensity_array, numpy.ndarray)
	assert all(
			numpy.equal(
					IonChromatogram(tic.intensity_array, tic.time_list).intensity_array,
					tic.intensity_array,
					)
			)

	ic = im.get_ic_at_index(0)
	tic.intensity_array = ic.intensity_array
	assert all(numpy.equal(tic.intensity_array, ic.intensity_array))

	assert isinstance(tic.intensity_array, numpy.ndarray)
	assert isinstance(tic.intensity_array[0], float)
	assert tic.intensity_array[0] == 0.0
	assert tic.intensity_array[2] == 622.0


def test_time_step(tic: IonChromatogram):
	assert isinstance(tic.time_step, float)
	assert tic.time_step == 1.0560000035830972


def test_write(tic: IonChromatogram, tmp_pathplus: PathPlus):
	tic.write(tmp_pathplus / "tic.dat", minutes=False, formatting=False)

	with (tmp_pathplus / "tic.dat").open() as fp:
		for line, ii in zip(fp.readlines(), range(len(tic.time_list))):
			assert line == f"{tic.time_list[ii]} {tic.intensity_array[ii]}\n"

	tic.write(tmp_pathplus / "tic_minutes.dat", minutes=True, formatting=False)

	with (tmp_pathplus / "tic_minutes.dat").open() as fp:
		for line, ii in zip(fp.readlines(), range(len(tic.time_list))):
			assert line == f"{tic.time_list[ii] / 60.0} {tic.intensity_array[ii]}\n"

	tic.write(tmp_pathplus / "tic_formatting.dat", minutes=False)

	with (tmp_pathplus / "tic_formatting.dat").open() as fp:
		for line, ii in zip(fp.readlines(), range(len(tic.time_list))):
			assert line == f"{tic.time_list[ii]:8.4f} {tic.intensity_array[ii]:#.6e}\n"

	for obj in [test_dict, *test_sequences, *test_numbers]:
		with pytest.raises(TypeError):
			tic.write(obj)  # type: ignore[arg-type]


# Inherited Methods from pymsBaseClass


def test_dump(im_i: IntensityMatrix, tmp_pathplus: PathPlus):
	im_i.dump(tmp_pathplus / "im_i_dump.dat")

	# Errors
	for obj in [*test_sequences, test_dict, *test_numbers]:
		with pytest.raises(TypeError):
			im_i.dump(obj)  # type: ignore[arg-type]

	# Read and check values
	assert (tmp_pathplus / "im_i_dump.dat").exists()
	loaded_im_i = cast(IntensityMatrix, _pickle_load_path(tmp_pathplus / "im_i_dump.dat"))
	assert loaded_im_i == im_i
	assert len(loaded_im_i) == len(im_i)


# Inherited Methods from TimeListMixin


def test_time_list(tic: IntensityMatrix):
	assert isinstance(tic.time_list, list)
	assert isinstance(tic.time_list[0], float)
	assert tic.time_list[0] == 1.05200003833
	assert len(tic.time_list) == 2103


# Inherited Methods from IntensityArrayMixin


def test_intensity_matrix(im: IntensityMatrix):
	with pytest.warns(DeprecationWarning, match="Use 'intensity_array' attribute instead"):
		assert isinstance(im.intensity_matrix, numpy.ndarray)
		assert isinstance(im.intensity_matrix[0], numpy.ndarray)
		assert isinstance(im.intensity_matrix[0][0], float)
		assert im.intensity_matrix[0][0] == 0.0
		assert im.intensity_matrix[3][5] == 0.0
		assert im.intensity_matrix[0][0] == im.intensity_array[0][0]
		assert numpy.equal(im.intensity_matrix.all(), im.intensity_array.all())


def test_intensity_array_list(im: IntensityMatrix):
	assert isinstance(im.intensity_array_list, list)
	assert isinstance(im.intensity_array_list[0], list)
	assert isinstance(im.intensity_array_list[0][0], float)
	assert im.intensity_array_list[0][0] == 0.0
	assert im.intensity_array_list[3][5] == 0.0
	assert im.intensity_array[0][0] == im.intensity_array_list[0][0]
	assert im.intensity_array_list == im.intensity_array.tolist()


def test_matrix_list(im: IntensityMatrix):
	with pytest.warns(DeprecationWarning, match="Use 'intensity_array_list' attribute instead"):
		assert isinstance(im.matrix_list, numpy.ndarray)


# Inherited methods from GetIndexTimeMixin


def test_get_index_at_time(tic: IonChromatogram):
	assert isinstance(tic.get_index_at_time(12), int)
	assert tic.get_index_at_time(12) == 10

	# Errors
	for obj in [test_string, *test_sequences, test_dict]:
		with pytest.raises(TypeError):
			tic.get_index_at_time(obj)  # type: ignore[arg-type]

	with pytest.raises(IndexError):
		tic.get_index_at_time(-1)

	with pytest.raises(IndexError):
		tic.get_index_at_time(1000000)


def test_get_time_at_index(tic: IonChromatogram):
	assert isinstance(tic.get_time_at_index(test_int), float)
	assert tic.get_time_at_index(test_int) == 1304.15599823

	# Errors
	with pytest.raises(TypeError):
		tic.get_time_at_index(test_string)  # type: ignore[arg-type]
	with pytest.raises(TypeError):
		tic.get_time_at_index(12.34)  # type: ignore[arg-type]
	with pytest.raises(TypeError):
		tic.get_time_at_index([1, 2, 3, 4])  # type: ignore[arg-type]
	with pytest.raises(TypeError):
		tic.get_time_at_index({'a': 1, 'b': 2, 'c': 3, 'd': 4})  # type: ignore[arg-type]
	# tic.get_time_at_index(0)

	with pytest.raises(IndexError):
		tic.get_time_at_index(-1)
	with pytest.raises(IndexError):
		tic.get_time_at_index(10000000)

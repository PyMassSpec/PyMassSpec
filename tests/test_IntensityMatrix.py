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
import pathlib
import pickle
import types

# 3rd party
import deprecation  # type: ignore
import numpy  # type: ignore
import pytest
from pytest_regressions.file_regression import FileRegressionFixture

# this package
from pyms.IntensityMatrix import (
		ASCII_CSV,
		IntensityMatrix,
		build_intensity_matrix,
		build_intensity_matrix_i,
		import_leco_csv
		)
from pyms.IonChromatogram import IonChromatogram
from pyms.Spectrum import MassSpectrum
from tests.constants import *


@pytest.fixture(scope="module")
def im_leco_filename(im, tmpdir_factory):
	"""
	Create the im_leco.csv file ahead of time and return the path to it
	"""

	filename = pathlib.Path(tmpdir_factory.mktemp('im_leco')) / 'im_leco.csv'
	im.export_leco_csv(filename)
	yield filename


class TestIntensityMatrix:

	def test_creation(self, im):
		assert isinstance(im, IntensityMatrix)

		IntensityMatrix(im.time_list, im.mass_list, im.intensity_array)

	args = [
			(test_string, TypeError),
			(test_int, TypeError),
			(test_float, TypeError),
			(test_dict, TypeError),
			]

	@pytest.mark.parametrize("obj, expects", [(test_list_strs, TypeError), *args])
	def test_time_list_errors(self, obj, im, expects):
		with pytest.raises(expects):
			IntensityMatrix(obj, im.mass_list, im.intensity_array)

	@pytest.mark.parametrize("obj, expects", [
			*args,
			([test_list_ints], ValueError),
			])
	def test_mass_list_errors(self, obj, im, expects):
		with pytest.raises(TypeError):
			IntensityMatrix(im.time_list, obj, im.intensity_array)

	@pytest.mark.parametrize("obj, expects", [
			*args,
			([test_list_ints], ValueError),
			])
	def test_intensity_array_errors(self, obj, im, expects):
		with pytest.raises(expects):
			IntensityMatrix(im.time_list, im.mass_list, obj)

	# Inherited Methods from pymsBaseClass

	def test_dump(self, im_i, tmp_pathplus):
		im_i.dump(tmp_pathplus / "im_i_dump.dat")

		# Errors
		for obj in [test_list_strs, test_dict, test_list_ints, test_tuple, *test_numbers]:
			with pytest.raises(TypeError):
				im_i.dump(obj)

		# Read and check values
		assert (tmp_pathplus / "im_i_dump.dat").exists()
		loaded_im_i = pickle.load((tmp_pathplus / "im_i_dump.dat").open("rb"))
		assert loaded_im_i == im_i
		assert len(loaded_im_i) == len(im_i)

	# Inherited Methods from TimeListMixin

	def test_time_list(self, im):
		time = im.time_list
		assert isinstance(time, list)
		# number of retention times
		assert len(time) == 2103
		# retention time of 1st scan:
		assert isinstance(time[0], float)
		assert time[0] == 1.05200003833

	# Inherited Methods from MassListMixin

	def test_mass_list(self, im):
		# get the list of masses (bin centers), and print the first ten
		assert isinstance(im.mass_list, list)
		assert isinstance(im.mass_list[0], float)
		assert im.mass_list[0] == 50.2516

	# Inherited Methods from MaxMinMassMixin

	def test_min_mass(self, im):
		# start mass
		assert isinstance(im.min_mass, float)
		assert im.min_mass == 50.2516

	def test_max_mass(self, im):
		# end mass
		assert isinstance(im.max_mass, float)
		assert im.max_mass == 499.2516

	# Inherited Methods from IntensityArrayMixin

	def test_intensity_array(self, im):
		assert isinstance(im.intensity_array, numpy.ndarray)
		assert isinstance(im.intensity_array[0], numpy.ndarray)
		assert isinstance(im.intensity_array[0][0], float)
		assert im.intensity_array[0][0] == 0.0
		assert im.intensity_array[2][3] == 1216.0
		print(im.intensity_array)

	def test_intensity_matrix(self, im):
		assert isinstance(im.intensity_matrix, numpy.ndarray)
		assert isinstance(im.intensity_matrix[0], numpy.ndarray)
		assert isinstance(im.intensity_matrix[0][0], float)
		assert im.intensity_matrix[0][0] == 0.0
		assert im.intensity_matrix[2][3] == 1216.0
		assert im.intensity_matrix[0][0] == im.intensity_array[0][0]
		assert numpy.equal(im.intensity_matrix.all(), im.intensity_array.all())

	def test_intensity_array_list(self, im):
		assert isinstance(im.intensity_array_list, list)
		assert all(isinstance(x, list) for x in im.intensity_array_list)
		assert all(isinstance(x, float) for x in im.intensity_array_list[0])
		assert im.intensity_array_list[0][0] == 0.0
		assert im.intensity_array_list[2][3] == 1216.0
		assert im.intensity_array[0][0] == im.intensity_array_list[0][0]
		assert im.intensity_array_list == im.intensity_array.tolist()

	def test_matrix_list(self, im):
		assert isinstance(im.matrix_list, numpy.ndarray)

	# Inherited methods from GetIndexTimeMixin

	def test_get_index_at_time(self, im):
		assert im.get_index_at_time(test_int) == 1168
		assert im.get_index_at_time(test_float) == 11

	@pytest.mark.parametrize(
			"obj, expects",
			[
					(test_string, TypeError),
					(test_dict, TypeError),
					(test_list_ints, TypeError),
					(test_list_strs, TypeError),
					(-1, IndexError),
					(1000000, IndexError),
					]
			)
	def test_get_index_at_time_errors(self, im, obj, expects):
		with pytest.raises(expects):
			im.get_index_at_time(obj)

	def test_get_time_at_index(self, im):
		assert im.get_time_at_index(test_int) == 1304.15599823

	@pytest.mark.parametrize(
			"obj, expects",
			[
					(test_string, TypeError),
					(test_dict, TypeError),
					(test_float, TypeError),
					(test_list_ints, TypeError),
					(test_list_strs, TypeError),
					(-1, IndexError),
					(1000000, IndexError),
					]
			)
	def test_get_time_at_index_errors(self, im, obj, expects):
		with pytest.raises(expects):
			im.get_time_at_index(obj)

	def test_len(self, im):
		assert len(im) == 2103

	def test_equality(self, im):
		assert im == IntensityMatrix(im.time_list, im.mass_list, im.intensity_array)
		assert im != test_string
		assert im != test_int
		assert im != test_float
		assert im != test_tuple
		assert im != test_list_ints
		assert im != test_list_strs

	@deprecation.fail_if_not_removed
	def test_get_local_size(self, im):
		with pytest.warns(DeprecationWarning):
			im.get_local_size()

	def test_local_size(self, im):
		assert isinstance(im.local_size, tuple)
		assert isinstance(im.local_size[0], int)
		assert im.local_size[0] == 2103

	@deprecation.fail_if_not_removed
	def test_get_size(self, im):
		with pytest.warns(DeprecationWarning):
			im.get_size()

	def test_size(self, im):
		# size of intensity matrix (#scans, #bins)
		assert isinstance(im.size, tuple)
		assert isinstance(im.size[0], int)
		assert im.size == (2103, 450)

	def test_iter_ms_indices(self, im):
		iter_ms = im.iter_ms_indices()
		assert isinstance(iter_ms, types.GeneratorType)
		for index, scan in enumerate(iter_ms):
			assert scan == index

	def test_iter_ic_indices(self, im):
		iter_ic = im.iter_ic_indices()
		assert isinstance(iter_ic, types.GeneratorType)
		for index, intensity in enumerate(iter_ic):
			assert intensity == index

	def test_set_ic_at_index(self, im):
		im = copy.deepcopy(im)

		im.set_ic_at_index(123, im.get_ic_at_index(0))
		assert im.get_ic_at_index(123).time_list == im.get_ic_at_index(0).time_list
		assert all(numpy.equal(im.get_ic_at_index(123).intensity_array, im.get_ic_at_index(0).intensity_array))

		for obj in [test_dict, test_list_strs, test_list_ints, test_string, test_float]:
			with pytest.raises(TypeError):
				im.set_ic_at_index(obj, im.get_ic_at_index(0))
			with pytest.raises(TypeError):
				im.set_ic_at_index(123, obj)

	def test_get_ic_at_index(self, im):
		ic = im.get_ic_at_index(123)

		# TODO: Check values for IC

		for obj in [test_dict, test_list_strs, test_list_ints, test_string, test_float]:
			with pytest.raises(TypeError):
				im.get_ic_at_index(obj)
		with pytest.raises(IndexError):
			im.get_ic_at_index(test_int)

	def test_get_ic_at_mass(self, im):
		# TODO: im.get_ic_at_mass() # Broken
		ic = im.get_ic_at_mass(123)

		assert isinstance(ic, IonChromatogram)
		assert not ic.is_tic()
		assert len(ic) == 2103
		assert isinstance(ic.intensity_array, numpy.ndarray)
		assert ic.get_time_at_index(test_int) == 1304.15599823
		assert ic.time_list[0] == 1.05200003833
		assert ic.mass == 123.2516
		assert ic.time_step == 1.0560000035830972
		assert ic.get_index_at_time(12) == 10

		for val in [test_int, 0, test_float]:
			with pytest.raises(IndexError):
				im.get_ic_at_mass(val)
		for obj in [test_dict, test_list_strs, test_list_ints, test_string]:
			with pytest.raises(TypeError):
				im.get_ic_at_mass(obj)

	def test_get_ms_at_index(self, im):
		ms = im.get_ms_at_index(123)
		assert isinstance(ms, MassSpectrum)

		assert isinstance(ms.mass_list, list)
		assert ms.mass_list[123] == 173.2516
		assert isinstance(ms.mass_spec, list)

		scan = im.get_scan_at_index(123)
		assert ms.mass_spec[123] == 0.0
		assert ms.mass_spec[123] == scan[123]
		assert ms.mass_spec[20] == 0.0
		assert ms.mass_spec[20] == scan[20]

		for obj in [test_dict, test_list_strs, test_list_ints, test_string]:
			with pytest.raises(TypeError):
				im.get_ms_at_index(obj)

	def test_get_scan_at_index(self, im):
		scan = im.get_scan_at_index(test_int)

		assert isinstance(scan, list)

		assert scan[123] == 0.0
		assert scan[20] == 2314.0

		for obj in [test_dict, test_list_strs, test_list_ints, test_string, test_float]:
			with pytest.raises(TypeError):
				im.get_scan_at_index(obj)
		with pytest.raises(IndexError):
			im.get_scan_at_index(-1)
		with pytest.raises(IndexError):
			im.get_scan_at_index(1000000)

	def test_mass_index(self, im):
		"""
		get_mass_at_index
		get_index_of_mass
		"""

		# the index of the nearest mass to 73.3m/z
		index = im.get_index_of_mass(73.3)
		assert isinstance(index, int)
		assert index == 23

		# the nearest mass to 73.3m/z
		assert isinstance(im.get_mass_at_index(index), float)
		assert im.get_mass_at_index(index) == 73.2516

		for obj in [test_string, test_list_strs, test_list_ints, test_dict]:
			with pytest.raises(TypeError):
				im.get_index_of_mass(obj)

		for obj in [test_float, test_string, test_list_strs, test_list_ints, test_dict]:
			with pytest.raises(TypeError):
				im.get_mass_at_index(obj)

		with pytest.raises(IndexError):
			im.get_mass_at_index(-1)
		with pytest.raises(IndexError):
			im.get_mass_at_index(1000000)

	def test_crop_mass(self, im):
		im = copy.deepcopy(im)

		for obj in [test_dict, *test_lists, test_string]:
			with pytest.raises(TypeError):
				im.crop_mass(obj, 200)

		for obj in [test_dict, *test_lists, test_string]:
			with pytest.raises(TypeError):
				im.crop_mass(100, obj)

		with pytest.raises(ValueError):
			im.crop_mass(200, 100)

		im.crop_mass(100, 200)

		with pytest.raises(ValueError):
			im.crop_mass(50, 200)
		with pytest.raises(ValueError):
			im.crop_mass(150, 500)

		im.crop_mass(101.5, 149.5)

	def test_null_mass(self, im):
		im = copy.deepcopy(im)

		for obj in [test_dict, *test_lists, test_string]:
			with pytest.raises(TypeError):
				im.null_mass(obj)

		with pytest.raises(IndexError):
			im.null_mass(500)
		with pytest.raises(IndexError):
			im.null_mass(10)

		im.null_mass(120)

		# TODO: Check that the nulling worked
		print(sum(im.get_ic_at_mass(120).intensity_array))

	def test_reduce_mass_spectra(self, im):
		im = copy.deepcopy(im)
		# TODO:

		for obj in [test_dict, *test_lists, test_string]:
			with pytest.raises(TypeError):
				im.reduce_mass_spectra(obj)


class Test_export_ascii:

	def test_export_ascii(self, im, tmp_pathplus, file_regression: FileRegressionFixture):
		"""
		Export the entire IntensityMatrix as CSV. This will create
		data.im.csv, data.mz.csv, and data.rt.csv where
		these are the intensity matrix, retention time
		vector, and m/z vector in the CSV format
		"""

		im.export_ascii(tmp_pathplus / "im_ascii")

		file_regression.check((tmp_pathplus / "im_ascii.im.dat").read_text(),
								encoding="UTF-8",
								extension="_im_ascii.im.dat")
		file_regression.check((tmp_pathplus / "im_ascii.mz.dat").read_text(),
								encoding="UTF-8",
								extension="_im_ascii.mz.dat")
		file_regression.check((tmp_pathplus / "im_ascii.rt.dat").read_text(),
								encoding="UTF-8",
								extension="_im_ascii.rt.dat")

		im.export_ascii(tmp_pathplus / "im_csv", fmt=ASCII_CSV)
		file_regression.check((tmp_pathplus / "im_csv.im.csv").read_text(),
								encoding="UTF-8",
								extension="_im_csv.im.csv")
		file_regression.check((tmp_pathplus / "im_csv.mz.csv").read_text(),
								encoding="UTF-8",
								extension="_im_csv.mz.csv")
		file_regression.check((tmp_pathplus / "im_csv.rt.csv").read_text(),
								encoding="UTF-8",
								extension="_im_csv.rt.csv")

	@pytest.mark.parametrize("obj", [test_dict, *test_lists, *test_numbers])
	def test_errors(self, obj, im, tmp_pathplus):
		with pytest.raises(TypeError):
			im.export_ascii(obj)

		with pytest.raises(ValueError):
			im.export_ascii(tmp_pathplus / "im_ascii", fmt=3)


class Test_leco_csv:

	def test_import_leco_csv(self, im, im_leco_filename):
		imported_im = import_leco_csv(im_leco_filename)
		assert isinstance(imported_im, IntensityMatrix)
		for imported, original in zip(imported_im.time_list, im.time_list):
			assert f"{imported:.3f}" == f"{original:.3f}"
		for imported, original in zip(imported_im.mass_list, im.mass_list):
			assert f"{imported:.0f}" == f"{original:.0f}"
		for imported1, original1 in zip(imported_im.intensity_array, im.intensity_array):
			for imported2, original2 in zip(imported1, original1):
				assert f"{imported2:.6e}" == f"{original2:.6e}"

		# Check size to original
		print("Output dimensions:", im.size, " Input dimensions:", imported_im.size)

	@pytest.mark.parametrize("obj", [test_dict, *test_lists, *test_numbers])
	def test_import_leco_csv_errors(self, im, im_leco_filename, obj):
		with pytest.raises(TypeError):
			import_leco_csv(obj)

	@pytest.mark.parametrize("obj", [test_dict, *test_lists, *test_numbers])
	def test_export_leco_csv_errors(self, im, im_leco_filename, obj):
		"""
		Export the entire IntensityMatrix as LECO CSV. This is
		useful for import into AnalyzerPro
		"""

		with pytest.raises(TypeError):
			im.export_leco_csv(obj)


def test_IntensityMatrix_custom(data):
	# IntensityMatrix
	# must build intensity matrix before accessing any intensity matrix methods.

	# bin interval of 0.5, eg. for double charge ions
	# intensity matrix, bin interval = 0.5, boundary +/- 0.25
	im = build_intensity_matrix(data, 0.5, 0.25, 0.25)
	assert isinstance(im, IntensityMatrix)

	# size of intensity matrix (#scans, #bins)
	assert isinstance(im.size, tuple)
	assert im.size == (2103, 900)

	# start mass
	assert isinstance(im.min_mass, float)
	assert im.min_mass == 50.2516

	# end mass
	assert isinstance(im.max_mass, float)
	assert im.max_mass == 499.7516

	# the index of the nearest mass to 73.3m/z
	index = im.get_index_of_mass(73.3)
	assert isinstance(index, int)
	assert index == 46

	# the nearest mass to 73.3m/z
	assert isinstance(im.get_mass_at_index(index), float)
	assert im.get_mass_at_index(index) == 73.2516

	# get the list of masses (bin centers), and print the first ten
	masses = im.mass_list
	assert isinstance(masses, list)
	assert masses[0] == 50.2516


def test_build_intensity_matrix(data):
	# todo

	for obj in [test_dict, *test_lists, test_string, *test_numbers]:
		with pytest.raises(TypeError):
			build_intensity_matrix(obj)  # type: ignore
	for obj in [test_dict, *test_lists, test_string]:
		with pytest.raises(TypeError):
			build_intensity_matrix(data, bin_interval=obj)  # type: ignore
	for obj in [test_dict, *test_lists, test_string]:
		with pytest.raises(TypeError):
			build_intensity_matrix(data, bin_left=obj)  # type: ignore
	for obj in [test_dict, *test_lists, test_string]:
		with pytest.raises(TypeError):
			build_intensity_matrix(data, bin_right=obj)  # type: ignore
	for obj in [test_dict, *test_lists, test_string]:
		with pytest.raises(TypeError):
			build_intensity_matrix(data, min_mass=obj)  # type: ignore
	with pytest.raises(ValueError):
		build_intensity_matrix(data, bin_interval=0)  # type: ignore


def test_build_intensity_matrix_i(data, im_i):
	assert isinstance(im_i, IntensityMatrix)

	# size of intensity matrix (#scans, #bins)
	assert isinstance(im_i.size, tuple)
	assert im_i.size == (2103, 450)

	# start mass
	assert isinstance(im_i.min_mass, int)
	assert im_i.min_mass == 50

	# end mass
	assert isinstance(im_i.max_mass, int)
	assert im_i.max_mass == 499

	# the index of the nearest mass to 73.3m/z
	index = im_i.get_index_of_mass(73.3)
	assert isinstance(index, int)
	assert index == 23

	# the nearest mass to 73.3m/z
	assert isinstance(im_i.get_mass_at_index(index), int)
	assert im_i.get_mass_at_index(index) == 73

	# get the list of masses (bin centers), and print the first ten
	masses = im_i.mass_list
	assert isinstance(masses, list)
	assert masses[0] == 50

	for obj in [test_dict, *test_lists, test_string, *test_numbers]:
		with pytest.raises(TypeError):
			build_intensity_matrix_i(obj)  # type: ignore
	for obj in [test_dict, *test_lists, test_string]:
		with pytest.raises(TypeError):
			build_intensity_matrix_i(data, bin_left=obj)  # type: ignore
	for obj in [test_dict, *test_lists, test_string]:
		with pytest.raises(TypeError):
			build_intensity_matrix_i(data, bin_right=obj)  # type: ignore


# TODO; Saving data
# # save the intensity matrix values to a file
# mat = im.matrix_list
# print("saving intensity matrix intensity values...")
# save_data("output/im.dat", mat)

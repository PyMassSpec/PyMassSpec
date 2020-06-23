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
import requests

# 3rd party
import pytest  # type: ignore
import numpy  # type: ignore

# pyms
from pyms.Spectrum import MassSpectrum
# tests
from .constants import *


def test_MassSpectrum(ms):
	assert isinstance(ms, MassSpectrum)


@pytest.mark.parametrize("obj, expects", [
		(test_list_ints, ValueError),
		(test_string, ValueError),
		(test_list_strs, ValueError),
		(test_int, TypeError),
		(test_dict, TypeError),
		])
def test_errors(ms, obj, expects):
	with pytest.raises(expects):
		MassSpectrum(obj, ms.intensity_list)

	with pytest.raises(expects):
		MassSpectrum(ms.mass_list, obj)


def test_len(ms):
	assert len(ms) == 450


def test_equality(im, ms):
	assert ms != im.get_ms_at_index(1234)
	assert ms == MassSpectrum(ms.mass_list, ms.mass_spec)
	assert ms != test_list_ints
	assert ms != test_list_strs
	assert ms != test_tuple
	assert ms != test_string
	assert ms != test_int
	assert ms != test_float


def test_mass_spec(ms):
	ms = copy.deepcopy(ms)
	assert ms.mass_spec[5] == 4192.0
	assert ms.mass_spec[50] == 3459.0
	assert ms.mass_spec[100] == 0.0

	ms.mass_spec[0] = 123
	assert ms.mass_spec[0] == 123

	assert ms.mass_spec == ms.intensity_list

	ms.mass_spec = list(range(len(ms.mass_spec)))
	assert ms.mass_spec == list(range(len(ms.mass_spec)))

	# for type in [test_list_ints, test_tuple]:
	# 	with pytest.raises(ValueError):
	# 		ms.mass_spec = type
	# 	with pytest.raises(ValueError):
	# 		ms.intensity_list = types


@pytest.mark.parametrize("obj, expects", [
		(test_string, ValueError),
		(test_dict, TypeError),
		(test_list_strs, ValueError),
		])
def test_mass_spec_errors(ms, obj, expects):
	with pytest.raises(expects):
		ms.mass_spec = obj

	with pytest.raises(expects):
		ms.intensity_list = obj


def test_mass_list(ms):
	assert ms.mass_list[5] == 55
	assert ms.mass_list[50] == 100
	assert ms.mass_list[100] == 150

	ms.mass_list = list(range(len(ms.mass_list)))
	assert list(ms.mass_list) == list(range(len(ms.mass_list)))

	# for obj in [test_list_ints, test_tuple]:
	# 	with pytest.raises(ValueError):
	# 		ms.mass_list = obj


@pytest.mark.parametrize("obj, expects", [
		(test_string, ValueError),
		(test_dict, TypeError),
		(test_list_strs, ValueError),
		])
def test_mass_list_errors(ms, obj, expects):
	with pytest.raises(expects):
		ms.mass_list = obj


def test_from_jcamp():
	nist_data_dir = pathlib.Path("nist_jdx_files")

	if not nist_data_dir.exists():
		nist_data_dir.mkdir(parents=True)

	# Compounds from nist
	for cas in [
			"122-39-4", "71-43-2", "85-98-3",
			"107-10-8", "50-37-3", "57-13-6",
			"77-92-9", "118-96-7"
			]:
		print(f"Testing CAS {cas}")
		jcamp_file = nist_data_dir / f"{cas}.jdx"

		if not jcamp_file.exists():
			r = requests.get(
					f"https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C{cas.replace('-', '')}&Index=0&Type=Mass")
			jcamp_file.write_bytes(r.content)

		MassSpectrum.from_jcamp(jcamp_file)

	# TODO: test jdx files from other sources


def test_from_mz_int_pairs():
	# Diphenylamine
	mz_int_pairs = [
			(27, 138),	(28, 210),	(32, 59),	(37, 70),	(38, 273),
			(39, 895),	(40, 141),	(41, 82),	(50, 710),	(51, 2151),
			(52, 434),	(53, 49),	(57, 41),	(59, 121),	(61, 73),
			(62, 229),	(63, 703),	(64, 490),	(65, 1106), (66, 932),
			(67, 68),	(70, 159),	(71, 266),	(72, 297),	(73, 44),
			(74, 263),	(75, 233),	(76, 330),	(77, 1636),	(78, 294),
			(84, 1732),	(87, 70),	(88, 86),	(89, 311),	(90, 155),
			(91, 219),	(92, 160),	(93, 107),	(101, 65),	(102, 111),
			(103, 99),	(104, 188),	(113, 107),	(114, 120),	(115, 686),
			(116, 150),	(117, 91),	(126, 46),	(127, 137),	(128, 201),
			(129, 73),	(130, 69),	(139, 447),	(140, 364),	(141, 584),
			(142, 279),	(143, 182),	(152, 37),	(153, 60),	(154, 286),
			(166, 718),
			(167, 3770),
			(168, 6825),
			(169, 9999),
			(170, 1210),
			(171, 85),
			]

	ms = MassSpectrum.from_mz_int_pairs(mz_int_pairs)

	assert isinstance(ms, MassSpectrum)
	assert len(ms) == len(mz_int_pairs)
	assert ms.intensity_list[6] == 141
	assert ms.intensity_list[30] == 1732
	assert ms.mass_list[-1] == 171
	assert ms.mass_list[2] == 32

	# Errors
	for obj in [test_string, test_int, test_list_strs, test_dict, test_list_ints, test_tuple, (["abc", "123"])]:
		with pytest.raises(TypeError):
			MassSpectrum.from_mz_int_pairs(obj)

	for obj in [[(1, 2, 3)], ([1, 2, 3],), [(1,)], ([1],), [("abc", "123")]]:
		with pytest.raises(ValueError):
			MassSpectrum.from_mz_int_pairs(obj)

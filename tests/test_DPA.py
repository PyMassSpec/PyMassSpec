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
import csv
import math
import operator
import pathlib
import tempfile
from typing import Any, Iterator, List

# 3rd party
import numpy
import pytest
from coincidence.regressions import AdvancedFileRegressionFixture
from domdf_python_tools.paths import PathPlus

# this package
from pyms.BillerBiemann import BillerBiemann, num_ions_threshold, rel_threshold
from pyms.DPA.Alignment import Alignment, exprl2alignment
from pyms.DPA.PairwiseAlignment import PairwiseAlignment, align_with_tree
from pyms.Experiment import Experiment, load_expr
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas
from pyms.Peak.List.Function import composite_peak
from pyms.Peak.List.IO import store_peaks
from pyms.TopHat import tophat
from pyms.Utils.Utils import is_number

# this package
from .constants import *

eley_codes = [
		"ELEY_1_SUBTRACT",
		"ELEY_2_SUBTRACT",
		"ELEY_3_SUBTRACT",
		"ELEY_4_SUBTRACT",
		"ELEY_5_SUBTRACT",
		]
geco_codes = [
		"GECO_1_SUBTRACT",
		"GECO_2_SUBTRACT",
		"GECO_3_SUBTRACT",
		"GECO_4_SUBTRACT",
		"GECO_5_SUBTRACT",
		]

# within replicates alignment parameters
Dw = 2.5  # rt modulation [s]
Gw = 0.30  # gap penalty


@pytest.fixture(scope="module")
def expr_list(pyms_datadir: PathPlus) -> Iterator[List[Experiment]]:

	with tempfile.TemporaryDirectory() as tmpdir:
		outputdir = pathlib.Path(tmpdir)

		# Create experiment files
		for jcamp_file in eley_codes:

			im = build_intensity_matrix_i(JCAMP_reader(pyms_datadir / f"{jcamp_file}.JDX"))

			# Intensity matrix size (scans, masses)
			n_scan, n_mz = im.size

			# noise filter and baseline correct
			for ii in range(n_mz):
				ic = im.get_ic_at_index(ii)
				ic_smooth = savitzky_golay(ic)
				ic_bc = tophat(ic_smooth, struct="1.5m")
				im.set_ic_at_index(ii, ic_bc)

			peak_list = BillerBiemann(im, points=9, scans=2)

			print('#')
			apl = rel_threshold(peak_list, 2)
			new_peak_list = num_ions_threshold(apl, 3, 3000)
			print('#')

			# ignore TMS ions and set mass range
			for peak in new_peak_list:
				peak.crop_mass(50, 400)
				peak.null_mass(73)
				peak.null_mass(147)

				# find area
				area = peak_sum_area(im, peak)
				peak.area = area
				area_dict = peak_top_ion_areas(im, peak)
				peak.ion_areas = area_dict

			expr = Experiment(jcamp_file, new_peak_list)

			# set time range for all experiments
			expr.sele_rt_range(["6.5m", "21m"])

			print('#')
			expr.dump(outputdir / f"{jcamp_file}.expr")
			print('#')

		# Load experiments
		expr_list = []
		for expr_code in eley_codes:
			expr = load_expr(outputdir / f"{expr_code}.expr")
			assert isinstance(expr, Experiment)
			expr_list.append(expr)

		yield expr_list


def test_expr_inequality(expr_list: List[Experiment]):
	assert expr_list[0] != expr_list[1]


@pytest.fixture(scope="module")
def F1(expr_list: List[Experiment]) -> List[Alignment]:
	# do the alignment
	print("Aligning ELEY SUBTRACT")
	F1 = exprl2alignment(expr_list)
	assert isinstance(F1, list)

	return F1


@pytest.fixture(scope="module")
def T1(F1: List[Alignment]) -> PairwiseAlignment:
	T1 = PairwiseAlignment(F1, Dw, Gw)
	assert isinstance(T1, PairwiseAlignment)

	return T1


@pytest.fixture(scope="module")
def A1(T1: PairwiseAlignment) -> Alignment:
	A1 = align_with_tree(T1, min_peaks=2)
	assert isinstance(A1, Alignment)

	assert isinstance(len(A1), int)
	assert len(A1) == 232
	A1.filter_min_peaks(5)
	assert len(A1) == 50

	return A1


class Test_alignment_Errors:

	@pytest.mark.parametrize("obj", [test_string, test_int, *test_sequences, test_dict])
	def test_alignments_errors(self, F1: List[Alignment], obj: Any):
		with pytest.raises(TypeError):
			PairwiseAlignment(F1, obj, Gw)
		with pytest.raises(TypeError):
			PairwiseAlignment(F1, Dw, obj)

	@pytest.mark.parametrize("obj", [*test_numbers, test_string, *test_sequences, test_dict])
	def test_expr_errors(self, obj: Any):
		with pytest.raises(TypeError):
			exprl2alignment(obj)
		with pytest.raises(TypeError):
			PairwiseAlignment(obj, Dw, Gw)
		with pytest.raises(TypeError):
			Alignment(obj)

	@pytest.mark.parametrize("obj", [test_float, test_string, test_dict, test_list_strs, test_list_ints])
	def test_min_peaks_errors(self, A1: Alignment, obj: Any):
		with pytest.raises(TypeError):
			A1.filter_min_peaks(obj)

	@pytest.mark.parametrize("obj", [*test_numbers, test_dict, *test_lists])
	def test_file_name_errors(self, A1: Alignment, obj: Any, tmp_pathplus: PathPlus):
		with pytest.raises(TypeError):
			A1.write_csv(obj, tmp_pathplus / "alignment_area.csv")
		with pytest.raises(TypeError):
			A1.write_csv(tmp_pathplus / "alignment_rt.csv", obj)
		with pytest.raises(TypeError):
			A1.write_common_ion_csv(obj, A1.common_ion())
		with pytest.raises(TypeError):
			A1.write_ion_areas_csv(obj)

	@pytest.mark.parametrize("obj", [*test_numbers, test_dict, test_list_strs, test_string])
	def test_top_ion_list_errors(self, A1: Alignment, obj: Any, tmp_pathplus: PathPlus):
		with pytest.raises(TypeError):
			A1.write_common_ion_csv(tmp_pathplus / "alignent_ion_area.csv", obj)


def test_write_csv(A1: Alignment, tmp_pathplus: PathPlus):
	A1.write_csv(tmp_pathplus / "alignment_rt.csv", tmp_pathplus / "alignment_area.csv")

	# Read alignment_rt.csv and alignment_area.csv and check values
	assert (tmp_pathplus / "alignment_rt.csv").exists()
	assert (tmp_pathplus / "alignment_area.csv").exists()

	with (tmp_pathplus / "alignment_rt.csv").open() as fp:
		rt_csv = list(csv.reader(fp))

	with (tmp_pathplus / "alignment_area.csv").open() as fp:
		area_csv = list(csv.reader(fp))

	assert rt_csv[0][0:2] == area_csv[0][0:2] == ["UID", "RTavg"]
	assert rt_csv[0][2:] == area_csv[0][2:] == A1.expr_code

	for peak_idx in range(len(A1.peakpos[0])):  # loop through peak lists (rows)

		new_peak_list = []

		for align_idx in range(len(A1.peakpos)):
			peak = A1.peakpos[align_idx][peak_idx]

			if peak is not None:

				if peak.rt is None or numpy.isnan(peak.rt):
					assert rt_csv[peak_idx + 1][align_idx + 2] == "NA"
				else:
					assert rt_csv[peak_idx + 1][align_idx + 2] == f"{peak.rt / 60:.3f}"

				if peak.area is None or numpy.isnan(peak.area):
					assert area_csv[peak_idx + 1][align_idx + 2] == "NA"
				else:
					assert area_csv[peak_idx + 1][align_idx + 2] == f"{peak.area:.0f}"

				new_peak_list.append(peak)

		compo_peak = composite_peak(new_peak_list)
		assert compo_peak is not None

		assert rt_csv[peak_idx + 1][0] == area_csv[peak_idx + 1][0] == compo_peak.UID

		assert rt_csv[peak_idx + 1][1] == area_csv[peak_idx + 1][1] == f"{float(compo_peak.rt / 60):.3f}"

	A1.write_csv(
			tmp_pathplus / "alignment_rt_seconds.csv",
			tmp_pathplus / "alignment_area_seconds.csv",
			minutes=False,
			)

	# Read alignment_rt_seconds.csv and alignment_area_seconds.csv and check values
	assert (tmp_pathplus / "alignment_rt_seconds.csv").exists()
	assert (tmp_pathplus / "alignment_area_seconds.csv").exists()

	with (tmp_pathplus / "alignment_rt_seconds.csv").open() as fp:
		rt_csv = list(csv.reader(fp))

	with (tmp_pathplus / "alignment_area_seconds.csv").open() as fp:
		area_csv = list(csv.reader(fp))

	assert rt_csv[0][0:2] == area_csv[0][0:2] == ["UID", "RTavg"]
	assert rt_csv[0][2:] == area_csv[0][2:] == A1.expr_code

	for peak_idx in range(len(A1.peakpos[0])):  # loop through peak lists (rows)

		new_peak_list = []

		for align_idx in range(len(A1.peakpos)):
			peak = A1.peakpos[align_idx][peak_idx]

			if peak is not None:

				if peak.rt is None or numpy.isnan(peak.rt):
					assert rt_csv[peak_idx + 1][align_idx + 2] == "NA"
				else:
					assert rt_csv[peak_idx + 1][align_idx + 2] == f"{peak.rt:.3f}"

				if peak.area is None or numpy.isnan(peak.area):
					assert area_csv[peak_idx + 1][align_idx + 2] == "NA"
				else:
					assert area_csv[peak_idx + 1][align_idx + 2] == f"{peak.area:.0f}"

				new_peak_list.append(peak)

		compo_peak = composite_peak(new_peak_list)

		assert compo_peak is not None

		assert rt_csv[peak_idx + 1][0] == area_csv[peak_idx + 1][0] == compo_peak.UID

		assert rt_csv[peak_idx + 1][1] == area_csv[peak_idx + 1][1] == f"{float(compo_peak.rt):.3f}"


def test_write_ion_areas_csv(A1: Alignment, tmp_pathplus: PathPlus):
	A1.write_ion_areas_csv(tmp_pathplus / "alignment_ion_areas.csv")
	A1.write_ion_areas_csv(tmp_pathplus / "alignment_ion_areas_seconds.csv", minutes=False)

	# Read alignment_ion_areas.csv and check values
	assert (tmp_pathplus / "alignment_ion_areas.csv").exists()

	with (tmp_pathplus / "alignment_ion_areas.csv").open() as fp:
		ion_csv = list(csv.reader(fp, delimiter='|'))

	with (tmp_pathplus / "alignment_ion_areas_seconds.csv").open() as fp:
		seconds_ion_csv = list(csv.reader(fp, delimiter='|'))

	assert ion_csv[0][0:2] == seconds_ion_csv[0][0:2] == ["UID", "RTavg"]
	assert ion_csv[0][2:] == seconds_ion_csv[0][2:] == A1.expr_code

	for peak_idx in range(len(A1.peakpos[0])):  # loop through peak lists (rows)

		new_peak_list = []

		for align_idx in range(len(A1.peakpos)):
			peak = A1.peakpos[align_idx][peak_idx]

			if peak is not None:
				ia = peak.ion_areas
				ia.update((mass, math.floor(intensity)) for mass, intensity in ia.items())
				sorted_ia = sorted(ia.items(), key=operator.itemgetter(1), reverse=True)

				assert ion_csv[peak_idx + 1][align_idx + 2] == str(sorted_ia)
				assert seconds_ion_csv[peak_idx + 1][align_idx + 2] == str(sorted_ia)

				new_peak_list.append(peak)

		compo_peak = composite_peak(new_peak_list)
		assert compo_peak is not None

		assert ion_csv[peak_idx + 1][0] == seconds_ion_csv[peak_idx + 1][0] == compo_peak.UID

		assert ion_csv[peak_idx + 1][1] == f"{float(compo_peak.rt / 60):.3f}"
		assert seconds_ion_csv[peak_idx + 1][1] == f"{float(compo_peak.rt):.3f}"


def test_write_common_ion_csv(
		A1: Alignment,
		tmp_pathplus: PathPlus,
		advanced_file_regression: AdvancedFileRegressionFixture,
		):
	common_ion = A1.common_ion()
	assert isinstance(common_ion, list)
	assert is_number(common_ion[0])
	assert common_ion[0] == 77

	# read the csv and check values
	A1.write_common_ion_csv(tmp_pathplus / "alignment_common_ion.csv", A1.common_ion())
	advanced_file_regression.check_file(
			tmp_pathplus / "alignment_common_ion.csv",
			extension="_alignment_common_ion.csv",
			)

	A1.write_common_ion_csv(tmp_pathplus / "alignment_common_ion_seconds.csv", A1.common_ion(), minutes=False)
	advanced_file_regression.check_file(
			tmp_pathplus / "alignment_common_ion_seconds.csv",
			extension="_alignment_common_ion_seconds.csv",
			)


# TODO: read the csv and check values


def test_align_2_alignments(A1: Alignment, pyms_datadir: PathPlus, tmp_pathplus: PathPlus):
	expr_list = []

	for jcamp_file in geco_codes:
		im = build_intensity_matrix_i(JCAMP_reader(pyms_datadir / f"{jcamp_file}.JDX"))

		# Intensity matrix size (scans, masses)
		n_scan, n_mz = im.size

		# noise filter and baseline correct
		for ii in range(n_mz):
			ic = im.get_ic_at_index(ii)
			ic_smooth = savitzky_golay(ic)
			ic_bc = tophat(ic_smooth, struct="1.5m")
			im.set_ic_at_index(ii, ic_bc)

		peak_list = BillerBiemann(im, points=9, scans=2)
		apl = rel_threshold(peak_list, 2)
		new_peak_list = num_ions_threshold(apl, 3, 3000)

		# ignore TMS ions and set mass range
		for peak in new_peak_list:
			peak.crop_mass(50, 400)
			peak.null_mass(73)
			peak.null_mass(147)

			# find area
			area = peak_sum_area(im, peak)
			peak.area = area
			area_dict = peak_top_ion_areas(im, peak)
			peak.ion_areas = area_dict

		expr = Experiment(jcamp_file, new_peak_list)

		# set time range for all experiments
		expr.sele_rt_range(["6.5m", "21m"])

		expr_list.append(expr)

	F2 = exprl2alignment(expr_list)
	T2 = PairwiseAlignment(F2, Dw, Gw)
	A2 = align_with_tree(T2, min_peaks=2)

	# top_ion_list = A2.common_ion()
	# A2.write_common_ion_csv(tmp_pathplus/'area1.csv', top_ion_list)

	# between replicates alignment parameters
	Db = 10.0  # rt modulation
	Gb = 0.30  # gap penalty

	print("Aligning input {1,2}")
	T9 = PairwiseAlignment([A1, A2], Db, Gb)
	A9 = align_with_tree(T9)

	A9.write_csv(tmp_pathplus / "rt.csv", tmp_pathplus / "area.csv")

	aligned_peaks = list(filter(None, A9.aligned_peaks()))
	store_peaks(aligned_peaks, tmp_pathplus / "peaks.bin")

	top_ion_list = A9.common_ion()
	A9.write_common_ion_csv(tmp_pathplus / "area.csv", top_ion_list)


# def test_alignment_compare():
# todo

# def test_dp()
# todo

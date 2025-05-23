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
from pathlib import Path
from typing import List

# 3rd party
import pytest

# this package
from pyms.BillerBiemann import BillerBiemann, num_ions_threshold, rel_threshold
from pyms.Experiment import Experiment
from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import IntensityMatrix, build_intensity_matrix, build_intensity_matrix_i
from pyms.IonChromatogram import IonChromatogram
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Class import Peak
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas
from pyms.Spectrum import MassSpectrum, Scan
from pyms.TopHat import tophat

pytest_plugins = ("coincidence", )


@pytest.fixture(scope="session")
def pyms_datadir() -> Path:
	return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def data(pyms_datadir: Path) -> GCMS_data:
	print("data")
	return JCAMP_reader(pyms_datadir / "ELEY_1_SUBTRACT.JDX")


@pytest.fixture(scope="session")
def im(data: GCMS_data) -> IntensityMatrix:
	# build an intensity matrix object from the data
	return build_intensity_matrix(data)


@pytest.fixture(scope="session")
def tic(data: GCMS_data) -> IonChromatogram:
	# get the TIC
	return copy.deepcopy(data.tic)


@pytest.fixture(scope="session")
def im_i(data: GCMS_data) -> IntensityMatrix:
	# build an intensity matrix object from the data
	return build_intensity_matrix_i(data)


@pytest.fixture(scope="session")  # noqa: PT005
def _peak_list(im_i: IntensityMatrix) -> List[Peak]:  # noqa: PT005
	im_i = copy.deepcopy(im_i)

	# Intensity matrix size (scans, masses)
	n_scan, n_mz = im_i.size

	# noise filter and baseline correct
	for ii in range(n_mz):
		ic = im_i.get_ic_at_index(ii)
		ic_smooth = savitzky_golay(ic)
		ic_bc = tophat(ic_smooth, struct="1.5m")
		im_i.set_ic_at_index(ii, ic_bc)

	# Use Biller and Biemann technique to find apexing ions at a scan
	# default is maxima over three scans and not to combine with any neighbouring
	# scan.
	peak_list = BillerBiemann(im_i, points=9, scans=2)
	return peak_list


@pytest.fixture()
def peak_list(_peak_list: List[Peak]) -> List[Peak]:
	return copy.deepcopy(_peak_list)


@pytest.fixture(scope="session")  # noqa: PT005
def _filtered_peak_list(im_i: IntensityMatrix, _peak_list: List[Peak]) -> List[Peak]:  # noqa: PT005
	peak_list = copy.deepcopy(_peak_list)
	# do peak detection on pre-trimmed data
	# trim by relative intensity
	apl = rel_threshold(peak_list, 2, copy_peaks=False)

	# trim by threshold
	new_peak_list = num_ions_threshold(apl, 3, 3000, copy_peaks=False)

	# ignore TMS ions and set mass range
	for peak in new_peak_list:
		peak.crop_mass(50, 400)
		peak.null_mass(73)
		peak.null_mass(147)

		# find area
		area = peak_sum_area(im_i, peak)
		peak.area = area
		area_dict = peak_top_ion_areas(im_i, peak)
		peak.ion_areas = area_dict

	return new_peak_list


@pytest.fixture()
def filtered_peak_list(_filtered_peak_list: List[Peak]) -> List[Peak]:
	return copy.deepcopy(_filtered_peak_list)


@pytest.fixture(scope="session")
def peak(im_i: IntensityMatrix) -> Peak:
	scan_i = im_i.get_index_at_time(31.17 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)
	return Peak(12.34, ms)


@pytest.fixture(scope="session")
def ms(im_i: IntensityMatrix) -> MassSpectrum:
	return copy.deepcopy(im_i.get_ms_at_index(0))


@pytest.fixture(scope="session")
def scan(data: GCMS_data) -> Scan:
	# return deepcopy(im_i.get_scan_at_index(0))
	return copy.deepcopy(data.scan_list[0])


@pytest.fixture()
def expr(filtered_peak_list: List[Peak]) -> Experiment:
	# create an experiment
	return Experiment("ELEY_1_SUBTRACT", filtered_peak_list)

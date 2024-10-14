"""
Functions to fill missing peak objects.
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
import pathlib
from typing import List, Optional

# 3rd party
import pandas  # type: ignore[import-untyped]
from domdf_python_tools.typing import PathLike
from enum_tools import IntEnum

# this package
from pyms.BillerBiemann import get_maxima_list_reduced
from pyms.Gapfill.Class import MissingPeak, Sample
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import ion_area
from pyms.TopHat import tophat
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Utils import is_path

__all__ = [
		"file2dataframe",
		"missing_peak_finder",
		"mp_finder",
		"write_filled_csv",
		"write_filled_rt_csv",
		"MZML",
		"NETCDF",
		"MissingPeakFiletype",
		]


class MissingPeakFiletype(IntEnum):
	"""
	Flag to indicate the filetype for :func:`pyms.Gapfill.Function.missing_peak_finder`.

	.. versionadded:: 2.3.0
	"""

	MZML = 1
	NETCDF = 2


MZML = MissingPeakFiletype.MZML
NETCDF = MissingPeakFiletype.NETCDF


def file2dataframe(file_name: PathLike) -> pandas.DataFrame:
	"""
	Convert a .csv file to a pandas DataFrame.

	:param file_name: CSV file to read.

	:authors: Jairus Bowne, Sean O'Callaghan, Dominic Davis-Foster (pathlib support)

	.. versionadded:: 2.3.0
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	return pandas.read_csv(
			file_name,
			delimiter=',',
			quotechar='"',
			header=0,
			)


def missing_peak_finder(
		sample: Sample,
		file_name: str,
		points: int = 3,
		null_ions: Optional[List] = None,
		crop_ions: Optional[List] = None,
		threshold: int = 1000,
		rt_window: float = 1,
		filetype: MissingPeakFiletype = MZML,
		) -> None:
	r"""
	Integrates raw data around missing peak locations to fill ``NA``\s in the data matrix.

	:param sample: The sample object containing missing peaks
	:param file_name: Name of the raw data file
	:param points: Peak finding - Peak if maxima over 'points' number of scans.
	:param null_ions: Ions to be deleted in the matrix.
	:default null_ions: ``[73, 147]``
	:param crop_ions: Range of Ions to be considered.
	:default crop_ions: ``[50, 540]``
	:param threshold: Minimum intensity of IonChromatogram allowable to fill.
	:param rt_window: Window in seconds around average RT to look for.
	:param filetype:

	:author: Sean O'Callaghan
	"""

	if not null_ions:
		null_ions = [73, 147]
	if not crop_ions:
		crop_ions = [50, 540]

	# TODO: some error checks on null and crop ions

	# TODO: a for root,files,dirs in os.path.walk(): loop
	print("Sample:", sample.name, "File:", file_name)

	if filetype == NETCDF:
		# this package
		from pyms.GCMS.IO.ANDI import ANDI_reader
		data = ANDI_reader(file_name)

	elif filetype == MZML:
		# this package
		from pyms.GCMS.IO.MZML import mzML_reader
		data = mzML_reader(file_name)

	else:
		print("file type not valid")

	# build integer intensity matrix
	im = build_intensity_matrix_i(data)

	for null_ion in null_ions:
		im.null_mass(null_ion)

	im.crop_mass(crop_ions[0], crop_ions[1])

	# get the size of the intensity matrix
	n_scan, n_mz = im.size

	# smooth data
	for ii in range(n_mz):
		ic = im.get_ic_at_index(ii)
		ic1 = savitzky_golay(ic, points)
		ic_smooth = savitzky_golay(ic1, points)
		ic_base = tophat(ic_smooth, struct="1.5m")
		im.set_ic_at_index(ii, ic_base)

	for mp in sample.missing_peaks:

		mp_rt = mp.rt
		common_ion = mp.common_ion
		qual_ion_1 = float(mp.qual_ion1)
		qual_ion_2 = float(mp.qual_ion2)

		ci_ion_chrom = im.get_ic_at_mass(common_ion)
		print("ci = ", common_ion)
		qi1_ion_chrom = im.get_ic_at_mass(qual_ion_1)
		print("qi1 = ", qual_ion_1)
		qi2_ion_chrom = im.get_ic_at_mass(qual_ion_2)
		print("qi2 = ", qual_ion_2)
		######
		# Integrate the CI around that particular RT
		#######

		# Convert time to points
		# How long between scans?

		points_1 = ci_ion_chrom.get_index_at_time(float(mp_rt))
		points_2 = ci_ion_chrom.get_index_at_time(float(mp_rt) - rt_window)
		print("rt_window = ", points_1 - points_2)

		rt_window_points = points_1 - points_2

		maxima_list = get_maxima_list_reduced(ci_ion_chrom, mp_rt, rt_window_points)

		large_peaks = []

		for rt, intens in maxima_list:
			if intens > threshold:
				q1_index = qi1_ion_chrom.get_index_at_time(rt)
				q2_index = qi2_ion_chrom.get_index_at_time(rt)

				q1_intensity = qi1_ion_chrom.get_intensity_at_index(q1_index)
				q2_intensity = qi2_ion_chrom.get_intensity_at_index(q2_index)

				if q1_intensity > threshold / 2 and q2_intensity > threshold / 2:
					large_peaks.append([rt, intens])

		print(f"found {len(large_peaks):d} peaks above threshold")

		areas = []
		for peak in large_peaks:
			apex = ci_ion_chrom.get_index_at_time(peak[0])
			ia = ci_ion_chrom.intensity_array.tolist()
			area, left, right, l_share, r_share = ion_area(ia, apex, 0)
			areas.append(area)

		########################

		areas.sort()
		if len(areas) > 0:
			biggest_area = areas[-1]
			mp.common_ion_area = biggest_area
			# mp.exact_rt = f"{float(mp_rt) / 60.0:.3f}"
			mp.exact_rt = float(mp_rt) / 60.0
			print("found area:", biggest_area, "at rt:", mp_rt)
		else:
			print("Missing peak at rt = ", mp_rt)
			mp.common_ion_area = None


def mp_finder(input_matrix: List) -> List[Sample]:
	r"""
	Finds the ``'NA'``\s in the transformed ``area_ci.csv`` file and makes
	:class:`pyms.Gapfill.Class.Sample` objects with them

	:param input_matrix: Data matrix derived from the ``area_ci.csv`` file.

	:rtype:

	:authors: Jairus Bowne, Sean O'Callaghan
	"""  # noqa: D400

	sample_list = []

	try:
		ci_pos = input_matrix[0].index(' "Quant Ion"')
		print("found Quant Ion position:", ci_pos)
	except ValueError:
		ci_pos = input_matrix[0].index('"Quant Ion"')

	uid_pos = input_matrix[0].index("UID")

	# Set up the sample objects
	# All entries on line 1 beyond the Qual Ion position are sample names
	for i, sample_name in enumerate(input_matrix[0][ci_pos:]):
		print(sample_name)
		sample = Sample(sample_name, i + 3)  # add 4 to allow for UID, RT,QualIon
		sample_list.append(sample)

	for line in input_matrix[1:]:
		uid = line[uid_pos]
		common_ion = line[ci_pos]

		qual_ion_1 = uid.split('-')[0]
		qual_ion_2 = uid.split('-')[1]
		rt = uid.split('-')[-1]

		for i, area in enumerate(line[ci_pos:]):
			if area == "NA":
				missing_peak = MissingPeak(common_ion, qual_ion_1, qual_ion_2, rt)
				sample_list[i].add_missing_peak(missing_peak)

	return sample_list


def write_filled_csv(
		sample_list: List[Sample],
		area_file: PathLike,
		filled_area_file: PathLike,
		) -> None:
	r"""
	Creates a new ``area_ci.csv`` file, replacing NAs with values from the sample_list objects where possible.

	:param sample_list:
	:param area_file: The file ``'area_ci.csv'`` from PyMassSpec output.
	:param filled_area_file: the new output file which has ``'NA'``\s values replaced.

	:authors: Jairus Bowne, Sean O'Callaghan, Dominic Davis-Foster
	"""

	if not is_path(filled_area_file):
		raise TypeError("'filled_area_file' must be a string or a pathlib.Path object")

	filled_area_file = prepare_filepath(filled_area_file)

	df = file2dataframe(area_file)

	uid_list: List[str] = df["UID"]
	rt_list: List[float] = []
	for uid in uid_list:
		rt = uid.split('-')[-1]
		rt_list.append(float(rt))

	for sample_name in df.columns[3:]:

		for sample in sample_list:
			if sample_name in sample.name:
				rt_area_dict = sample.rt_areas
				break
		else:
			raise ValueError(f"Sample {sample_name!r} not found in sample_list.")

		for i, part in enumerate(df[sample_name]):
			if part == "NA":
				try:
					df[sample_name][i] = rt_area_dict[rt_list[i]]
				except KeyError:
					pass

	df.to_csv(filled_area_file, index=False, na_rep="NA")


def write_filled_rt_csv(
		sample_list: List[Sample],
		rt_file: PathLike,
		filled_rt_file: PathLike,
		) -> None:
	r"""
	Creates a new rt.csv file, replacing ``'NA'``\s with values from the sample_list objects where possible.

	:param sample_list: A list of samples.
	:param rt_file: the file ``rt.csv`` from PyMassSpec output.
	:param filled_rt_file: the new output file which has ``NA`` values replaced.

	:authors: Jairus Bowne, Sean O'Callaghan, Dominic Davis-Foster
	"""

	if not isinstance(filled_rt_file, (str, pathlib.Path)):
		raise TypeError("'filled_rt_file' must be a string or a pathlib.Path object")

	filled_rt_file = prepare_filepath(filled_rt_file)

	df = file2dataframe(rt_file)

	uid_list: List[str] = df["UID"]
	rt_list: List[str] = []
	for uid in uid_list:
		rt = uid.split('-')[-1]
		rt_list.append(rt)

	for sample_name in df.columns[3:]:
		for sample in sample_list:
			if sample_name in sample.name:
				rt_exact_rt_dict = sample.get_mp_rt_exact_rt_dict()
				break
		else:
			raise ValueError(f"Sample {sample_name!r} not found in sample_list.")

		for i, part in enumerate(df[sample_name]):
			if part == "NA":
				try:
					df[sample_name][i] = rt_exact_rt_dict[float(rt_list[i])]
				except KeyError:
					pass

	df.to_csv(filled_rt_file, index=False, na_rep="NA")

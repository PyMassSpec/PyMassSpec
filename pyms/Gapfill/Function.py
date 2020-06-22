"""
Functions to fill missing peak objects
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
import csv
import pathlib

# 3rd party
import numpy

# this package
from pyms.Utils.Utils import is_path
from pyms.BillerBiemann import get_maxima_list_reduced
from pyms.Gapfill.Class import MissingPeak, Sample
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import ion_area
from pyms.TopHat import tophat
from pyms.Utils.IO import prepare_filepath

MZML = 1
NETCDF = 2


# .csv reader (cloned from gcqc project)
def file2matrix(file_name: Union[str, pathlib.Path) -> numpy.ndarray:
	"""
	Convert a .csv file to a numpy array

	:param file_name: Filename (.csv) to convert (area.csv, area_ci.csv)
	:type file_name: str or pathlib.Path

	:return: Data matrix
	:rtype: :class:`numpy.array`

	:author: Jairus Bowne
	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	with file_name.open() as fp:
		reader = csv.reader(fp, delimiter=",", quotechar='"')
		matrix = []
		for row in reader:
			newrow = []
			for each in row:
				try:
					each = float(each)
				except:
					pass
				newrow.append(each)
			matrix.append(newrow)

	return numpy.array(matrix)


def missing_peak_finder(sample: Sample, file_name: str, points: int = 3, null_ions: Optional[List] = None,
						crop_ions: Optional[List] = None, threshold: int = 1000, rt_window: float = 1, filetype: int = MZML):
	"""
	Integrates raw data around missing peak locations to fill NAs in the data matrix

	:param sample: The sample object containing missing peaks
	:type sample: :class:`pyms.Gapfill.Class.Sample`

	:param file_name: Name of the raw data file
	:type file_name: str
	:param points: Peak finding - Peak if maxima over 'points' number of scans (Default 3)
	:type points: int, optional
	:param  null_ions: Ions to be deleted in the matrix (Default [73, 147])
	:type null_ions: list, optional
	:param crop_ions: Range of Ions to be considered (Default [50, 540])
	:type crop_ions: list, optional
	:param threshold: Minimum intensity of IonChromatogram allowable to fill (Default 1000)
	:type threshold: int, optional
	:param  rt_window: Window in seconds around average RT to look for (Default 1)
	:type rt_window: float, optional
	:param filetype: either `MZML` (default) or `NETCDF`
	:type filetype: int, optional

	:author: Sean O'Callaghan
	"""

	if not null_ions:
		null_ions = [73, 147]
	if not crop_ions:
		crop_ions = [50, 540]

	# TODO: some error checks on null and crop ions

	# TODO: a for root,files,dirs in os.path.walk(): loop
	print("Sample:", sample.get_name(), "File:", file_name)

	if filetype.lower() == 'cdf':
		from pyms.GCMS.IO.ANDI import ANDI_reader
		data = ANDI_reader(file_name)
	elif filetype.lower() == 'mzml':
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

	for mp in sample.get_missing_peaks():

		mp_rt = mp.rt
		common_ion = mp.get_ci()
		qual_ion_1 = float(mp.get_qual_ion1())
		qual_ion_2 = float(mp.get_qual_ion2())

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

		maxima_list = get_maxima_list_reduced(
			ci_ion_chrom, mp_rt, rt_window_points
		)

		large_peaks = []

		for rt, intens in maxima_list:
			if intens > threshold:
				q1_index = qi1_ion_chrom.get_index_at_time(rt)
				q2_index = qi2_ion_chrom.get_index_at_time(rt)

				q1_intensity = qi1_ion_chrom.get_intensity_at_index(q1_index)
				q2_intensity = qi2_ion_chrom.get_intensity_at_index(q2_index)

				if q1_intensity > threshold / 2 and q2_intensity > threshold / 2:
					large_peaks.append([rt, intens])

		print(f'found {len(large_peaks):d} peaks above threshold')

		areas = []
		for peak in large_peaks:
			apex = ci_ion_chrom.get_index_at_time(peak[0])
			ia = ci_ion_chrom.get_intensity_array().tolist()
			area, left, right, l_share, r_share = ion_area(ia, apex, 0)
			areas.append(area)

		########################

		areas.sort()
		if len(areas) > 0:
			biggest_area = areas[-1]
			mp.set_ci_area(biggest_area)
			mp.set_exact_rt(f"{float(mp_rt) / 60.0:.3f}")
			print("found area:", biggest_area, "at rt:", mp_rt)
		else:
			print("Missing peak at rt = ", mp_rt)
			mp.set_ci_area('na')


def mp_finder(input_matrix: List) -> Sample:
	"""
	Finds the 'NA's in the transformed area_ci.csv file and makes
	:class:`pyms.Gapfill.Class.Sample` objects with them

	:param input_matrix: Data matrix derived from the area_ci.csv file
	:type input_matrix: list

	:return: list of Samples
	:rtype: :class:`list` of :class:`pyms.Gapfill.Class.Sample` objects

	:author: Jairus Bowne
	:author: Sean O'Callaghan
	"""

	sample_list = []

	try:
		ci_pos = input_matrix[0].index(' "Quant Ion"')
		print("found Quant Ion position:", ci_pos)
	except ValueError:
		ci_pos = input_matrix[0].index('"Quant Ion"')

	uid_pos = input_matrix[0].index('UID')

	# Set up the sample objects
	# All entries on line 1 beyond the Qual Ion position are sample names
	for i, sample_name in enumerate(input_matrix[0][ci_pos:]):
		print(sample_name)
		sample = Sample(sample_name, i + 3)  # add 4 to allow for UID, RT,QualIon
		sample_list.append(sample)

	for line in input_matrix[1:]:
		uid = line[uid_pos]
		common_ion = line[ci_pos]

		qual_ion_1 = uid.split("-")[0]
		qual_ion_2 = uid.split("-")[1]
		rt = uid.split("-")[-1]

		for i, area in enumerate(line[ci_pos:]):
			if area == 'NA':
				missing_peak = MissingPeak(common_ion, qual_ion_1, qual_ion_2, rt)
				sample_list[i].add_missing_peak(missing_peak)

	return sample_list


def transposed(lists: List) -> List[lists]:
	"""
	transposes a list of lists

	:param lists: the list of lists to be transposed
	:type lists: list

	:return: transposed list of lists
	:rtype: :class:`list` of lists

	:author: Jairus Bowne
	:author: Sean O'Callaghan
	"""

	if not lists:
		return []

	return map(lambda *row: list(row), *lists)


def write_filled_csv(sample_list: Sample, area_file: Union[str, pathlib.Path], filled_area_file: Union[str, pathlib.Path]):
	"""
	creates a new area_ci.csv file, replacing NAs with values from the sample_list objects where possible

	:param sample_list: A list of samples
	:type sample_list: :class:`list` of :class:`pyms.Gapfill.Class.Sample` objects
	:param area_file: the file 'area_ci.csv' from PyMassSpec output
	:type area_file: str or pathlib.Path
	:param filled_area_file: the new output file which has NA values replaced
	:type filled_area_file: str or pathlib.Path

	:author: Jairus Bowne
	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (pathlib support)
	"""

	if not isinstance(filled_area_file, (str, pathlib.Path)):
		raise TypeError("'filled_area_file' must be a string or a pathlib.Path object")

	if not isinstance(filled_area_file, pathlib.Path):
		filled_area_file = pathlib.Path(filled_area_file)

	if not filled_area_file.parent.is_dir():
		filled_area_file.parent.mkdir(parents=True)

	old_matrix = file2matrix(area_file)

	# Invert it to be a little more efficient
	invert_old_matrix = zip(*old_matrix)
	# print invert_old_matrix[0:5]

	uid_list = invert_old_matrix[0][1:]
	rt_list = []
	for uid in uid_list:
		rt = uid.split('-')[-1]
		rt_list.append(rt)

	# print(rt_list)

	# start setting up the output file
	invert_new_matrix = []
	for line in invert_old_matrix[0:2]:
		invert_new_matrix.append(line)

	for line in invert_old_matrix[3:]:
		sample_name = line[0]

		new_line = []
		new_line.append(sample_name)
		for sample in sample_list:
			if sample_name in sample.get_name():
				rt_area_dict = sample.get_mp_rt_area_dict()
				# print rt_area_dict

		for i, part in enumerate(line[1:]):
			# print part
			if part == 'NA':
				try:
					area = rt_area_dict[str(rt_list[i])]
					new_line.append(area)
				except KeyError:
					pass
			else:
				new_line.append(part)

		invert_new_matrix.append(new_line)

	fp_new = filled_area_file.open('w')

	#    new_matrix = numpy.empty(matrix_size)
	new_matrix = transposed(invert_new_matrix)

	for i, line in enumerate(new_matrix):
		for j, part in enumerate(line):
			fp_new.write(f"{part},")
		fp_new.write("\n")

	fp_new.close()


def write_filled_rt_csv(sample_list: Sample, rt_file: Union[str, pathlib.Path], filled_rt_file: Union[str, pathlib.Path]):
	"""
	creates a new rt.csv file, replacing NAs with values from the sample_list objects where possible

	:param sample_list: A list of samples
	:type sample_list: :class:`list` of :class:`pyms.Gapfill.Class.Sample` objects
	:param rt_file: the file 'rt.csv' from PyMassSpec output
	:type rt_file: str or pathlib.Path
	:param filled_rt_file: the new output file which has NA values replaced
	:type filled_rt_file: str or pathlib.Path

	:author: Jairus Bowne
	:author: Sean O'Callaghan
	:author: Dominic Davis-Foster (pathlib support)
	"""

	if not isinstance(filled_rt_file, (str, pathlib.Path)):
		raise TypeError("'filled_rt_file' must be a string or a pathlib.Path object")

	if not isinstance(filled_rt_file, pathlib.Path):
		filled_rt_file = pathlib.Path(filled_rt_file)

	if not filled_rt_file.parent.is_dir():
		filled_rt_file.parent.mkdir(parents=True)

	old_matrix = file2matrix(rt_file)

	# Invert it to be a little more efficent
	invert_old_matrix = zip(*old_matrix)

	uid_list = invert_old_matrix[0][1:]
	rt_list = []
	for uid in uid_list:
		rt = uid.split('-')[-1]
		rt_list.append(rt)

	# start setting up the output file
	invert_new_matrix = []
	for line in invert_old_matrix[0:1]:
		invert_new_matrix.append(line)

	for line in invert_old_matrix[2:]:
		sample_name = line[0]

		new_line = [sample_name]
		for sample in sample_list:
			if sample_name in sample.get_name():

				rt_exact_rt_dict = sample.get_mp_rt_exact_rt_dict()

		for i, part in enumerate(line[1:]):
			if part == 'NA':
				try:
					rt_new = rt_exact_rt_dict[str(rt_list[i])]
					new_line.append(rt_new)
				except KeyError:
					pass

			else:
				new_line.append(part)

		invert_new_matrix.append(new_line)

	fp_new = open(filled_rt_file, 'w')

	# new_matrix = numpy.empty(matrix_size)
	new_matrix = transposed(invert_new_matrix)

	for i, line in enumerate(new_matrix):
		for j, part in enumerate(line):
			fp_new.write(f"{str(part),}")
		fp_new.write("\n")

	fp_new.close()

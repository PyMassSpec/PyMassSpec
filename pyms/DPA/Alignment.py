"""
Classes for peak alignment by dynamic programming.
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
import copy
import math
import operator
import pathlib
from typing import Dict, List, Optional, Sequence

# 3rd party
import numpy
import pandas  # type: ignore[import]
from domdf_python_tools.typing import PathLike
from pandas import DataFrame

# this package
from pyms.Experiment import Experiment
from pyms.Peak import Peak
from pyms.Peak.List.Function import composite_peak
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Utils import _number_types, is_path, is_sequence, is_sequence_of

__all__ = ["Alignment", "exprl2alignment"]

# Ensure that the intersphinx links are correct.
DataFrame.__module__ = "pandas"


class Alignment:
	"""
	Models an alignment of peak lists.

	:param expr: The experiment to be converted into an alignment object.

	:authors: Woon Wai Keen, Qiao Wang, Vladimir Likic, Dominic Davis-Foster.
	"""

	#:
	peakpos: List[List[Peak]]

	#: List of experiment codes.
	expr_code: List[str]

	#:
	peakalgt: List[List[Peak]]

	#:
	similarity: Optional[float]

	def __init__(self, expr: Optional[Experiment]):

		if expr is None:
			self.peakpos = []
			self.peakalgt = []
			self.expr_code = []
			self.similarity = None
		elif not isinstance(expr, Experiment):
			raise TypeError("'expr' must be an 'Experiment' object")
		else:
			# for peak in expr.peak_list:
			#    if peak.area() == None or peak.area() <= 0:
			#        error("All peaks must have an area for alignment")
			self.peakpos = [copy.deepcopy(expr.peak_list)]
			self.peakalgt = numpy.transpose(self.peakpos).tolist()  # type: ignore[arg-type]
			self.expr_code = [expr.expr_code]
			self.similarity = None

	def __len__(self) -> int:
		"""
		Returns the length of the alignment, defined as the number of
		peak positions in the alignment.

		:rtype:

		:authors: Qiao Wang, Vladimir Likic
		"""  # noqa: D400

		return len(self.peakalgt)

	def aligned_peaks(self, minutes: bool = False) -> Sequence[Optional[Peak]]:
		"""
		Returns a list of Peak objects where each peak has the combined spectra
		and average retention time of all peaks that aligned.

		:param minutes: Whether retention times are in minutes.
			If :py:obj:`False`, retention time are in seconds.

		:return: A list of composite peaks based on the alignment.

		:author: Andrew Isaac

		.. TODO:: minutes currently does nothing
		"""  # noqa: D400

		# for all peaks found
		peak_list = []

		for peak_idx in range(len(self.peakpos[0])):
			# get aligned peaks, ignore missing

			new_peak_list = []
			for align_idx in range(len(self.peakpos)):
				peak = self.peakpos[align_idx][peak_idx]

				if peak is not None:
					new_peak_list.append(peak)

			# create composite
			new_peak = composite_peak(new_peak_list)
			peak_list.append(new_peak)

		return peak_list

	def common_ion(self) -> List[float]:
		"""
		Calculates a common ion among the peaks of an aligned peak.

		:return: A list of the highest intensity common ion for all aligned peaks.

		:author: Sean O'Callaghan
		"""

		# print("#peak lists =", len(self.peakpos))
		# print("#peaks =", len(self.peakpos[0]))

		# a list to contain the dictionaries
		# each dictionary contains the
		# top ions and their frequency for each peak
		# in the alignment
		list_of_top_ion_dicts: List = []
		empty_count_list: List[int] = []

		for peak_list in self.peakpos:
			# (re)initialise the peak index
			index = 0

			for peak in peak_list:
				# if the dict has not been created, create it
				# will only run on first peak list
				if len(list_of_top_ion_dicts) < len(peak_list):
					list_of_top_ion_dicts.append({})

				# copy the list entry to a local dict for code clarity
				top_ion_dict = list_of_top_ion_dicts[index]

				# count the empty peaks
				if (len(empty_count_list)) < len(peak_list):
					empty_count = 0
				else:
					empty_count = empty_count_list[index]
				# make sure the peak is present
				if peak is not None:
					# get the top 5 ions
					top_5 = peak.ion_areas.keys()

					for ion in top_5:
						if ion in top_ion_dict:
							# if we have seen it, increment the count
							top_ion_dict[ion] += 1
						else:
							# if we haven't seen it before, add it
							top_ion_dict[ion] = 1

					empty_count += 1

					# copy the dict back to the list
					list_of_top_ion_dicts[index] = top_ion_dict

					if len(empty_count_list) < len(peak_list):
						empty_count_list.append(empty_count)
					else:
						empty_count_list[index] = empty_count
				# increment for the next peak
				index += 1

		# print("length of list of dicts", len(list_of_top_ion_dicts))

		# index = 0
		# for entry in list_of_top_ion_dicts:
		#    for ion in sorted(entry, key=entry.get, reverse=True)[0:3]:
		#        print(ion,":", entry[ion],end='')
		#    print('  empty:', empty_count_list[index])
		#    index += 1

		top_ion_list = []

		for entry in list_of_top_ion_dicts:
			# This was in the Repo and commented in easyGC
			# top_ion_list.append(self.get_highest_mz_ion(entry))
			# This was in the version being used by GSMatch and easyGC
			top_ion_list.append(max(entry, key=entry.get))

		# TODO: Change over to the new version

		return top_ion_list

	def filter_min_peaks(self, min_peaks: int) -> None:
		"""
		Filters alignment positions that have less peaks than ``min_peaks``.

		This function is useful only for within state alignment.

		:param min_peaks: Minimum number of peaks required for the alignment
			position to survive filtering.

		:author: Qiao Wang
		"""

		if not isinstance(min_peaks, int):
			raise TypeError("'min_peaks' must be an integer")

		filtered_list = []

		for pos in range(len(self.peakalgt)):
			# if len(filter(None,self.peakalgt[pos])) >= min_peaks:
			if len([x for x in self.peakalgt[pos] if x is not None]) >= min_peaks:
				filtered_list.append(self.peakalgt[pos])

		self.peakalgt = filtered_list
		self.peakpos = numpy.transpose(self.peakalgt).tolist()  # type: ignore[arg-type]

	@staticmethod
	def get_highest_mz_ion(ion_dict: Dict[float, int]) -> float:
		"""
		Returns the preferred ion for quantitiation.

		Looks at the list of candidate ions, selects those which have
		highest occurrence, and selects the heaviest of those.

		:param ion_dict: a dictionary of *m/z* value: number of occurrences.

		:return ion: The ion with the highest *m/z* value.
		"""

		max_occurrences = max(ion_dict.values())
		most_freq_mzs = []

		for key, value in ion_dict.items():
			if value == max_occurrences:
				most_freq_mzs.append(key)

		return max(most_freq_mzs)

	def write_csv(
			self,
			rt_file_name: PathLike,
			area_file_name: PathLike,
			minutes: bool = True,
			) -> None:
		"""
		Writes the alignment to CSV files.

		This function writes two files: one containing the alignment of peak
		retention times and the other containing the alignment of peak areas.

		:param rt_file_name: The name for the retention time alignment file.
		:param area_file_name: The name for the areas alignment file.
		:param minutes: Whether to save retention times in minutes.
			If :py:obj:`False`, retention time will be saved in seconds.

		:authors: Woon Wai Keen, Andrew Isaac, Vladimir Likic, David Kainer,
			Dominic Davis-Foster (pathlib support)
		"""

		if not isinstance(rt_file_name, (str, pathlib.Path)):
			raise TypeError("'rt_file_name' must be a string or a pathlib.Path object")

		if not isinstance(area_file_name, (str, pathlib.Path)):
			raise TypeError("'area_file_name' must be a string or a pathlib.Path object")

		rt_file_name = prepare_filepath(rt_file_name)
		area_file_name = prepare_filepath(area_file_name)

		with rt_file_name.open('w', encoding="UTF-8") as fp1, area_file_name.open('w', encoding="UTF-8") as fp2:

			# create header
			header = ["UID", "RTavg"]
			for item in self.expr_code:
				header.append(f'"{item}"')

			# write headers
			fp1.write(','.join(header) + '\n')
			fp2.write(','.join(header) + '\n')

			# for each alignment position write alignment's peak and area
			for peak_idx in range(len(self.peakpos[0])):  # loop through peak lists (rows)

				rts = []
				areas = []
				new_peak_list = []

				for align_idx in range(len(self.peakpos)):
					peak = self.peakpos[align_idx][peak_idx]

					if peak is not None:

						if minutes:
							rt = peak.rt / 60.0
						else:
							rt = peak.rt

						rts.append(rt)
						areas.append(peak.area)
						new_peak_list.append(peak)

					else:
						rts.append(None)
						areas.append(None)

				compo_peak = composite_peak(new_peak_list)
				if compo_peak is None:
					continue

				# write to retention times file
				fp1.write(compo_peak.UID)

				if minutes:
					fp1.write(f",{float(compo_peak.rt / 60):.3f}")
				else:
					fp1.write(f",{compo_peak.rt:.3f}")

				for rt in rts:
					if rt is None or numpy.isnan(rt):
						fp1.write(",NA")
					else:
						fp1.write(f",{rt:.3f}")
				fp1.write('\n')

				# write to peak areas file
				fp2.write(compo_peak.UID)

				if minutes:
					fp2.write(f",{float(compo_peak.rt / 60):.3f}")
				else:
					fp2.write(f",{compo_peak.rt:.3f}")

				for area in areas:
					if area is None:
						fp2.write(",NA")
					else:
						fp2.write(f",{area:.0f}")
				fp2.write('\n')

	def write_common_ion_csv(
			self,
			area_file_name: PathLike,
			top_ion_list: Sequence[float],
			minutes: bool = True,
			) -> None:
		"""
		Writes the alignment to CSV files.

		This function writes two files: one containing the alignment of peak
		retention times and the other containing the alignment of peak areas.

		:param area_file_name: The name for the areas alignment file.
		:param top_ion_list: A list of the highest intensity common ion along the aligned peaks.
		:param minutes: Whether to save retention times in minutes.
			If :py:obj:`False`, retention time will be saved in seconds.

		:authors: Woon Wai Keen, Andrew Isaac, Sean O'Callaghan, Vladimir Likic,
			Dominic Davis-Foster (pathlib support)

		.. TODO:: minutes currently does nothing
		"""

		if not is_path(area_file_name):
			raise TypeError("'area_file_name' must be a string or a PathLike object")

		if not is_sequence_of(top_ion_list, _number_types):
			raise TypeError(f"'top_ion_list' must be a Sequence of numbers")

		area_file_name = prepare_filepath(area_file_name)

		with area_file_name.open('w') as fp:

			# create header
			header = ['"UID"', '"RTavg"', '"Quant Ion"']
			for item in self.expr_code:
				header.append(f'"{item}"')

			# write headers
			fp.write(','.join(header) + '\n')

			rtsums: List[float] = []
			rtcounts = []

			# The following two arrays will become list of lists
			# such that:
			# areas = [  [align1_peak1, align2_peak1, .....,alignn_peak1]
			#            [align1_peak2, ................................]
			#              .............................................
			#            [align1_peakm,....................,alignn_peakm]  ]
			areas: List[List] = []
			new_peak_lists: List[List[Peak]] = []

			for peak_list in self.peakpos:
				index = 0
				for peak in peak_list:
					# one the first iteration, populate the lists
					if len(areas) < len(peak_list):
						areas.append([])
						new_peak_lists.append([])
						rtsums.append(0)
						rtcounts.append(0)

					if peak is not None:
						rt = peak.rt

						# get the area of the common ion for the peak
						# an area of 'na' shows that while the peak was
						# aligned, the common ion was not present
						area = peak.get_ion_area(top_ion_list[index])

						areas[index].append(area)
						new_peak_lists[index].append(peak)

						# The following code to the else statement is
						# just for calculating the average rt
						rtsums[index] += rt
						rtcounts[index] += 1

					else:
						areas[index].append(None)

					index += 1

			out_strings = []
			index = 0
			# now write the strings for the file
			for area_list in areas:

				# write initial info:
				# peak unique id, peak average rt
				compo_peak = composite_peak(new_peak_lists[index])
				if compo_peak is None:
					continue

				peak_UID = compo_peak.UID
				peak_UID_string = f'"{peak_UID}"'

				rt_avg = rtsums[index] / rtcounts[index]

				out_strings.append(f"{peak_UID_string},{rt_avg / 60:.3f},{top_ion_list[index]:f}")

				for area in area_list:
					if area is not None:
						out_strings[index] += f",{area:.4f}"
					else:
						out_strings[index] += ",NA"

				index += 1

			# now write the file
			#        print("length of areas[0]", len(areas[0]))
			#        print("length of areas", len(areas))
			#        print("length of out_strings", len(out_strings))
			for row in out_strings:
				fp.write(row + '\n')

	def write_ion_areas_csv(self, ms_file_name: PathLike, minutes: bool = True) -> None:
		"""
		Write Ion Areas to CSV File.

		:param ms_file_name: The name of the file
		:param minutes: Whether to save retention times in minutes.
			If :py:obj:`False`, retention time will be saved in seconds.

		:authors: David Kainer, Dominic Davis-Foster (pathlib support)
		"""

		if not is_path(ms_file_name):
			raise TypeError("'ms_file_name' must be a string or a PathLike object")

		ms_file_name = prepare_filepath(ms_file_name)

		with ms_file_name.open('w') as fp1:

			# create header

			header = ['"UID"', '"RTavg"']
			for item in self.expr_code:
				header.append(f'"{item}"')

			# write headers
			fp1.write('|'.join(header) + '\n')

			for peak_idx in range(len(self.peakpos[0])):

				ias = []
				new_peak_list = []

				for align_idx in range(len(self.peakpos)):

					peak = self.peakpos[align_idx][peak_idx]

					if peak is not None:

						ia = peak.ion_areas
						ia.update((mass, math.floor(intensity)) for mass, intensity in ia.items())
						sorted_ia = sorted(ia.items(), key=operator.itemgetter(1), reverse=True)
						ias.append(sorted_ia)
						new_peak_list.append(peak)

				compo_peak = composite_peak(new_peak_list)
				if compo_peak is None:
					continue

				# write to ms file
				fp1.write(compo_peak.UID)

				if minutes:
					fp1.write(f"|{compo_peak.rt/60:.3f}")
				else:
					fp1.write(f"|{compo_peak.rt:.3f}")

				for intensity_array in ias:
					if intensity_array is None:
						fp1.write("|NA")
					else:
						fp1.write(f"|{intensity_array}")

				fp1.write('\n')

	def get_peak_alignment(
			self,
			minutes: bool = True,
			require_all_expr: bool = True,
			) -> DataFrame:
		"""
		Returns a Pandas dataframe of aligned retention times.

		:param minutes: Whether to return retention times in minutes.
			If :py:obj:`False`, retention time will be returned in seconds.
		:param require_all_expr: Whether the peak must be present in
			all experiments to be included in the data frame.

		:authors: Woon Wai Keen, Andrew Isaac, Vladimir Likic, Dominic Davis-Foster
		"""

		rt_table = []

		# for each alignment position write alignment's RT
		for peak_idx in range(len(self.peakpos[0])):
			rts = []
			countrt = 0

			for align_idx in range(len(self.peakpos)):
				peak = self.peakpos[align_idx][peak_idx]

				if peak is not None:

					if minutes:
						rt = peak.rt / 60.0
					else:
						rt = peak.rt

					rts.append(rt)
					countrt = countrt + 1

				else:
					rts.append(None)

			if (require_all_expr and countrt == len(self.expr_code)) or not require_all_expr:
				rt_table.append(rts)

		rt_alignment = pandas.DataFrame(rt_table, columns=self.expr_code)
		rt_alignment = rt_alignment.reindex(sorted(rt_alignment.columns), axis=1)

		return rt_alignment

	def get_ms_alignment(self, require_all_expr: bool = True) -> DataFrame:
		"""
		Returns a Pandas dataframe of mass spectra for the aligned peaks.

		:param require_all_expr: Whether the peak must be present in
			all experiments to be included in the data frame.

		:authors: Woon Wai Keen, Andrew Isaac, Vladimir Likic, Dominic Davis-Foster
		"""

		ms_table = []

		# for each alignment position write alignment's ms
		for peak_idx in range(len(self.peakpos[0])):
			specs = []
			countms = 0

			for align_idx in range(len(self.peakpos)):
				peak = self.peakpos[align_idx][peak_idx]

				if peak is not None:
					ms = peak.mass_spectrum
					specs.append(ms)
					countms = countms + 1
				else:
					specs.append(None)

			if (require_all_expr and countms == len(self.expr_code)) or not require_all_expr:
				ms_table.append(specs)

		ms_alignment = pandas.DataFrame(ms_table, columns=self.expr_code)
		ms_alignment = ms_alignment.reindex(sorted(ms_alignment.columns), axis=1)

		return ms_alignment

	def get_peaks_alignment(self, require_all_expr: bool = True) -> DataFrame:
		"""
		Returns a Pandas dataframe of Peak objects for the aligned peaks.

		:param require_all_expr: Whether the peak must be present in
			all experiments to be included in the data frame.

		:authors: Woon Wai Keen, Andrew Isaac, Vladimir Likic, Dominic Davis-Foster
		"""

		peaks_table = []

		# for each alignment position write alignment's ms
		for peak_idx in range(len(self.peakpos[0])):
			peaks = []
			count_peaks = 0

			for align_idx in range(len(self.peakpos)):
				peak = self.peakpos[align_idx][peak_idx]

				if peak is not None:
					peaks.append(peak)
					count_peaks = count_peaks + 1
				else:
					peaks.append(None)

			if (require_all_expr and count_peaks == len(self.expr_code)) or not require_all_expr:
				peaks_table.append(peaks)

		peak_alignment = pandas.DataFrame(peaks_table, columns=self.expr_code)
		peak_alignment = peak_alignment.reindex(sorted(peak_alignment.columns), axis=1)

		return peak_alignment

	def get_area_alignment(self, require_all_expr: bool = True) -> DataFrame:
		"""
		Returns a Pandas dataframe containing the peak areas of the aligned peaks.

		:param require_all_expr: Whether the peak must be present in all experiments
			to be included in the data frame.

		:authors: Woon Wai Keen, Andrew Isaac, Vladimir Likic, Dominic Davis-Foster
		"""

		areas_table = []

		# for each alignment position write alignment's ms
		for peak_idx in range(len(self.peakpos[0])):
			areas = []
			count_areas = 0

			for align_idx in range(len(self.peakpos)):
				peak = self.peakpos[align_idx][peak_idx]

				if peak is not None:
					area = peak.area
					areas.append(area)
					count_areas = count_areas + 1
				else:
					areas.append(None)

			if (require_all_expr and count_areas == len(self.expr_code)) or not require_all_expr:
				areas_table.append(areas)

		area_alignment = pandas.DataFrame(areas_table, columns=self.expr_code)
		area_alignment = area_alignment.reindex(sorted(area_alignment.columns), axis=1)

		return area_alignment


def exprl2alignment(expr_list: List[Experiment]) -> List[Alignment]:
	"""
	Converts a list of experiments into a list of alignments.

	:param expr_list: The list of experiments to be converted into an alignment objects.

	:return: A list of alignment objects for the experiments.

	:author: Vladimir Likic
	"""

	if not is_sequence(expr_list):
		raise TypeError("'expr_list' must be a Sequence")

	alignments = []

	for item in expr_list:
		if not isinstance(item, Experiment):
			raise TypeError("list items must be 'Experiment' instances")

		alignments.append(Alignment(item))

	return alignments

"""
Class to model GC-MS data.
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
import pathlib
from statistics import mean, median, stdev
from typing import List, Optional, Sequence, TypeVar, cast

# 3rd party
import numpy
from domdf_python_tools.doctools import prettify_docstrings
from domdf_python_tools.typing import PathLike

# this package
from pyms.Base import pymsBaseClass
from pyms.IonChromatogram import IonChromatogram
from pyms.Mixins import GetIndexTimeMixin, MaxMinMassMixin, TimeListMixin
from pyms.Spectrum import MassSpectrum, Scan
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Time import time_str_secs
from pyms.Utils.Utils import _number_types, is_path, is_sequence_of, signedinteger

__all__ = ["GCMS_data", "IntStr"]

MassSpectrum = MassSpectrum  # For legacy imports. Stops PyCharm complaining TODO: Remove eventually.

IntStr = TypeVar("IntStr", int, str)


@prettify_docstrings
class GCMS_data(pymsBaseClass, TimeListMixin, MaxMinMassMixin, GetIndexTimeMixin):
	"""
	Generic object for GC-MS data.

	Contains the raw data as a list of scans and a list of times.

	:param time_list: Scan retention times.
	:param scan_list:

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic,
		Dominic Davis-Foster (type assertions and properties)
	"""

	def __init__(self, time_list: Sequence[float], scan_list: Sequence[Scan]):
		if not is_sequence_of(time_list, _number_types):
			raise TypeError("'time_list' must be a Sequence of numbers")

		if not is_sequence_of(scan_list, Scan):
			raise TypeError("'scan_list' must be a Sequence of Scan objects")

		self._time_list = list(time_list)
		self._scan_list = list(scan_list)
		self._set_time()
		self._set_min_max_mass()
		self._calc_tic()

	def __eq__(self, other) -> bool:  # noqa: MAN001
		"""
		Return whether this GCMS_data object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return self.scan_list == other.scan_list and self.time_list == other.time_list

		return NotImplemented

	def __len__(self) -> int:
		"""
		Returns the length of the data object, defined as the number of scans.

		:author: Vladimir Likic
		"""

		return len(self._scan_list)

	def __repr__(self) -> str:
		return f"<GCMS_data({self.min_rt} - {self.max_rt} seconds, time step {self.time_step}, {len(self)} scans)>"

	def __str__(self) -> str:
		return self.__repr__()

	def _calc_tic(self) -> None:
		"""
		Calculate the total ion chromatogram.

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
		"""

		intensities = []
		for scan in self._scan_list:
			intensities.append(sum(scan.intensity_list))
		ia = numpy.array(intensities)
		rt = copy.deepcopy(self._time_list)
		tic = IonChromatogram(ia, rt)

		self._tic = tic

	def _set_time(self) -> None:
		"""
		Sets time-related properties of the data.

		:author: Vladimir Likic
		"""

		# calculate the time step, its spread, and along the way
		# check that retention times are increasing
		time_diff_list = []

		for index, t1 in enumerate(self._time_list):
			if index == len(self._time_list) - 1:
				break
			t2 = self._time_list[index + 1]
			if not t2 > t1:
				raise ValueError("Retention times are not in ascending order!")
			time_diff = t2 - t1
			time_diff_list.append(time_diff)

		time_step = mean(time_diff_list)
		time_step_std = stdev(time_diff_list)

		self._time_step = time_step
		self._time_step_std = time_step_std
		self._min_rt = min(self._time_list)
		self._max_rt = max(self._time_list)

	def _set_min_max_mass(self) -> None:
		"""
		Sets the min and max mass value.

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
		"""

		min_mass = None
		max_mass = None

		for scan in self._scan_list:

			if min_mass is None:
				min_mass = scan.min_mass
			elif scan.min_mass is not None:
				min_mass = min(scan.min_mass, min_mass)

			if max_mass is None:
				max_mass = scan.max_mass
			elif scan.max_mass is not None:
				max_mass = max(scan.max_mass, max_mass)

		self._min_mass = min_mass
		self._max_mass = max_mass

	def info(self, print_scan_n: bool = False) -> None:
		"""
		Prints some information about the data.

		:param print_scan_n: If set to :py:obj:`True` will print the number of *m/z* values in each scan.

		:author: Vladimir Likic
		"""

		# print the summary of simply attributes
		print(f" Data retention time range: {self._min_rt / 60.0:.3f} min -- {self._max_rt / 60:.3f} min")
		print(f" Time step: {self._time_step:.3f} s (std={self._time_step_std:.3f} s)")
		print(f" Number of scans: {len(self._scan_list):d}")
		print(f" Minimum m/z measured: {self._min_mass:.3f}")
		print(f" Maximum m/z measured: {self._max_mass:.3f}")

		# calculate median number of m/z values measured per scan
		n_list = []
		for ii in range(len(self._scan_list)):
			scan = self._scan_list[ii]
			n = len(scan)
			n_list.append(n)
			if print_scan_n:
				print(n)
		mz_mean = mean(n_list)
		mz_median = median(n_list)
		print(f" Mean number of m/z values per scan: {mz_mean:.0f}")
		print(f" Median number of m/z values per scan: {mz_median:.0f}")

	@property
	def scan_list(self) -> List[Scan]:
		"""
		Return a list of the scan objects.

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
		"""

		return copy.deepcopy(self._scan_list)

	@property
	def time_list(self) -> List[float]:
		"""
		Return a copy of the time list.
		"""

		return self._time_list[:]

	@property
	def tic(self) -> IonChromatogram:
		"""
		Returns the total ion chromatogram.

		:author: Andrew Isaac
		"""

		return self._tic

	@property
	def min_rt(self) -> float:
		"""
		Returns the minimum retention time for the data in seconds.
		"""

		return self._min_rt

	@property
	def max_rt(self) -> float:
		"""
		Returns the maximum retention time for the data in seconds.
		"""

		return self._max_rt

	@property
	def time_step(self) -> float:
		"""
		Returns the time step of the data.
		"""

		return self._time_step

	@property
	def time_step_std(self) -> float:
		"""
		Returns the standard deviation of the time step of the data.
		"""

		return self._time_step_std

	def trim(
			self,
			begin: Optional[IntStr] = None,
			end: Optional[IntStr] = None,
			) -> None:
		"""
		Trims data in the time domain.

		The arguments ``begin`` and ``end`` can be either integers (in which case
		they are taken as the first/last scan number for trimming) or strings
		in which case they are treated as time strings and converted to scan
		numbers.

		At present both ``begin`` and ``end`` must be of the same type, either both
		scan numbers or time strings.

		At least one of ``begin`` and ``end`` is required.

		:param begin: The start time or scan number
		:param end: The end time or scan number

		:author: Vladimir Likic
		"""

		# trim called with defaults, or silly arguments
		if begin is None and end is None:
			raise SyntaxError("At least one of 'begin' and 'end' is required")

		N = len(self._scan_list)

		# process 'begin' and 'end'
		if begin is None:
			first_scan = 0
		elif isinstance(begin, (int, signedinteger)):
			first_scan = cast(int, begin) - 1
		elif isinstance(begin, str):
			time = time_str_secs(begin)
			scan_ = self.get_index_at_time(time)
			if scan_ is None:
				raise TypeError("invalid 'begin' argument")
			first_scan = scan_ + 1
		else:
			raise TypeError("invalid 'begin' argument")

		if end is None:
			last_scan = N - 1
		elif isinstance(end, (int, signedinteger)):
			last_scan = cast(int, end)
		elif isinstance(end, str):
			time = time_str_secs(end)
			scan_ = self.get_index_at_time(time)
			if scan_ is None:
				raise TypeError("invalid 'end' argument")
			last_scan = scan_ + 1
		else:
			raise TypeError("invalid 'end' argument")

		# sanity checks
		if not last_scan > first_scan:
			raise ValueError(f"last scan={last_scan:d}, first scan={first_scan:d}")
		elif first_scan < 0:
			raise ValueError("scan number must be greater than one")
		elif last_scan > N - 1:
			raise ValueError(f"last scan={last_scan:d}, total number of scans={N:d}")

		print(f"Trimming data to between {first_scan + 1:d} and {last_scan + 1:d} scans")

		scan_list_new = []
		time_list_new = []
		for ii in range(len(self._scan_list)):
			if first_scan <= ii <= last_scan:
				scan = self._scan_list[ii]
				time = self._time_list[ii]
				scan_list_new.append(scan)
				time_list_new.append(time)

		# update info
		self._scan_list = scan_list_new
		self._time_list = time_list_new
		self._set_time()
		self._set_min_max_mass()
		self._calc_tic()

	def write(self, file_root: PathLike) -> None:
		"""
		Writes the entire raw data to two CSV files:

		- :file:`{<file_root>}.I.csv`, containing the intensities; and
		- :file:`{<file_root>}.mz.csv`, containing the corresponding *m/z* values.

		In general these are not two-dimensional matrices, because different
		scans may have different numbers of *m/z* values recorded.

		:param file_root: The root for the output file names

		:authors: Vladimir Likic, Dominic Davis-Foster (pathlib support)
		"""  # noqa: D400

		if not isinstance(file_root, (str, pathlib.Path)):
			raise TypeError("'file_root' must be a string or a pathlib.Path object")

		file_root = prepare_filepath(file_root)

		file_name1 = str(file_root) + ".I.csv"
		file_name2 = str(file_root) + ".mz.csv"

		print(f" -> Writing intensities to '{file_name1}'")
		print(f" -> Writing m/z values to '{file_name2}'")

		with open(file_name1, 'w', encoding="UTF-8") as fp1, open(file_name2, 'w', encoding="UTF-8") as fp2:

			for scan in self._scan_list:

				for index, intensity in enumerate(scan.intensity_list):
					if index == 0:
						fp1.write(f"{intensity:.4f}")
					else:
						fp1.write(f",{intensity:.4f}")
				fp1.write('\n')

				for index, mass in enumerate(scan.mass_list):
					if index == 0:
						fp2.write(f"{mass:.4f}")
					else:
						fp2.write(f",{mass:.4f}")
				fp2.write('\n')

	def write_intensities_stream(self, file_name: PathLike) -> None:
		"""
		Loop over all scans and, for each scan, write the intensities to the
		given file, one intensity per line.

		Intensities from different scans are joined without any delimiters.

		:param file_name: Output file name.

		:authors: Vladimir Likic, Dominic Davis-Foster (pathlib support)
		"""  # noqa: D400

		if not is_path(file_name):
			raise TypeError("'file_name' must be a string or a PathLike object")

		file_name = prepare_filepath(file_name)

		# n = len(self._scan_list)

		print(" -> Writing scans to a file")

		with file_name.open('w', encoding="UTF-8") as fp:

			for scan in self._scan_list:
				intensities = scan.intensity_list
				for i in intensities:
					fp.write(f"{i:8.4f}\n")

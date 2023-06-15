"""
Class to model Intensity Matrix.
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
from typing import Iterator, List, Optional, Sequence, Tuple, Union
from warnings import warn

# 3rd party
import numpy
from domdf_python_tools.typing import PathLike
from enum_tools import IntEnum, document_enum

# this package
from pyms.Base import pymsBaseClass
from pyms.GCMS.Class import GCMS_data
from pyms.IonChromatogram import BasePeakChromatogram, IonChromatogram
from pyms.Mixins import GetIndexTimeMixin, IntensityArrayMixin, MassListMixin, TimeListMixin
from pyms.Spectrum import MassSpectrum
from pyms.Utils.IO import prepare_filepath, save_data
from pyms.Utils.Utils import _number_types, is_number, is_path, is_sequence, is_sequence_of

__all__ = [
		"BaseIntensityMatrix",
		"AsciiFiletypes",
		"IntensityMatrix",
		"import_leco_csv",
		"build_intensity_matrix",
		"build_intensity_matrix_i",
		"ASCII_DAT",
		"ASCII_CSV",
		]


@document_enum
class AsciiFiletypes(IntEnum):
	"""
	Enumeration of supported ASCII filetypes for :meth:`~pyms.IntensityMatrix.IntensityMatrix.export_ascii`.

	.. versionadded:: 2.3.0
	"""

	ASCII_DAT = 1  # doc: Tab-delimited ASCII file
	ASCII_CSV = 0  # doc: Comma-separated values file


ASCII_DAT = AsciiFiletypes.ASCII_DAT
ASCII_CSV = AsciiFiletypes.ASCII_CSV


class BaseIntensityMatrix(pymsBaseClass, TimeListMixin, MassListMixin, IntensityArrayMixin, GetIndexTimeMixin):
	"""
	Base class for intensity matrices of binned raw data.

	:param time_list: Retention time values
	:param mass_list: Binned mass values
	:param intensity_array: List of lists of binned intensity values per scan

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions and properties)
	"""

	# Rows are scans, columns are masses

	_min_mass: float
	_max_mass: float

	def __init__(
			self,
			time_list: Sequence[float],
			mass_list: Sequence[float],
			intensity_array: Union[Sequence[Sequence[float]], numpy.ndarray],
			):
		# sanity check
		if not is_sequence_of(time_list, _number_types):
			raise TypeError("'time_list' must be a Sequence of numbers")

		if not is_sequence_of(mass_list, _number_types):
			raise TypeError("'mass_list' must be a Sequence of numbers")

		if not is_sequence(intensity_array) or not is_sequence_of(intensity_array[0], _number_types):
			raise TypeError("'intensity_array' must be a Sequence, of Sequences, of numbers")

		if not isinstance(intensity_array, numpy.ndarray):
			intensity_array = numpy.array(intensity_array)

		if not len(time_list) == len(intensity_array):
			raise ValueError("'time_list' is not the same length as 'intensity_array'")

		if not len(mass_list) == len(intensity_array[0]):
			raise ValueError("'mass_list' is not the same size as 'intensity_array'")

		self._time_list = list(time_list)
		self._mass_list = list(mass_list)

		self._intensity_array = intensity_array

		self._min_rt = min(time_list)
		self._max_rt = max(time_list)

		self._min_mass = min(mass_list)
		self._max_mass = max(mass_list)

	def __len__(self) -> int:
		"""
		Returns the number of scans in the intensity matrix.
		"""

		return len(self.time_list)

	def __eq__(self, other) -> bool:  # noqa: MAN001
		"""
		Return whether this intensity matrix object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return (
					self.time_list == other.time_list and self.mass_list == other.mass_list
					and numpy.array_equal(self.intensity_array, other.intensity_array)
					)

		return NotImplemented

	@property
	def size(self) -> Tuple[int, int]:
		"""
		Gets the size of intensity matrix.

		:return: Number of rows and cols

		:authors: Qiao Wang, Andrew Isaac, Luke Hodkinson, Vladimir Likic
		"""

		n_scan = len(self._intensity_array)
		n_mz = len(self._intensity_array[0])

		return n_scan, n_mz

	def iter_ms_indices(self) -> Iterator[int]:
		"""
		Iterates over row indices.
		"""

		yield from range(0, len(self._intensity_array))

	def iter_ic_indices(self) -> Iterator[int]:
		"""
		Iterate over column indices.
		"""

		yield from range(0, len(self._intensity_array[0]))

	def set_ic_at_index(self, ix: int, ic: IonChromatogram) -> None:
		"""
		Sets the intensity of the mass at index ``ix`` in each scan to a new value.

		:param ix: Index of an ion chromatogram in the intensity data matrix to be set
		:param ic: Ion chromatogram that will be copied at position ``ix`` in the data matrix

		The length of the ion chromatogram must match the appropriate
		dimension of the intensity matrix.

		:author: Vladimir Likic
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an an integer")

		if not isinstance(ic, IonChromatogram):
			raise TypeError("'ic' must be an IonChromatogram object")

		ia: numpy.ndarray = ic.intensity_array

		# check if the dimension is ok
		if len(ia) != len(self._intensity_array):
			raise ValueError("ion chromatogram incompatible with the intensity matrix")

		for row_idx, intensity in enumerate(ia):
			self._intensity_array[row_idx][ix] = intensity

	def get_ic_at_index(self, ix: int) -> IonChromatogram:
		"""
		Returns the ion chromatogram at the specified index.

		:param ix: Index of an ion chromatogram in the intensity data matrix.

		:return: Ion chromatogram at given index.

		:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an integer")

		ia = [intensities[ix] for intensities in self._intensity_array]

		ic_ia = numpy.array(ia)
		mass = self.get_mass_at_index(ix)
		rt = self._time_list[:]

		return IonChromatogram(ic_ia, rt, mass)

	def get_ms_at_index(self, ix: int) -> MassSpectrum:
		"""
		Returns a mass spectrum for a given scan index.

		:param ix: The index of the scan.

		:author: Andrew Isaac
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an an integer")

		scan = self.get_scan_at_index(ix)

		return MassSpectrum(self.mass_list, scan)

	def get_scan_at_index(self, ix: int) -> List[float]:
		"""
		Returns the spectral intensities for scan index.

		:param ix: The index of the scan

		:return: Intensity values of scan spectra

		:author: Andrew Isaac
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an an integer")

		if ix < 0 or ix >= len(self._intensity_array):
			raise IndexError("index out of range")

		return self._intensity_array[ix].tolist()

	def get_mass_at_index(self, ix: int) -> float:
		"""
		Returns binned mass at index.

		:param ix: Index of binned mass

		:return: Binned mass

		:author: Andrew Isaac
		"""

		if not isinstance(ix, int):
			raise TypeError("'ix' must be an an integer")

		if ix < 0 or ix >= len(self._mass_list):
			raise IndexError("index out of range")

		return self._mass_list[ix]

	def get_index_of_mass(self, mass: float) -> int:
		"""
		Returns the index of the nearest binned mass to the given mass.

		:param mass: Mass to lookup in list of masses

		:author: Andrew Isaac
		"""

		if not is_number(mass):
			raise TypeError("'mass' must be a number")

		best = self._max_mass
		ix = 0
		for idx, mass_ in enumerate(self._mass_list):
			tmp = abs(mass_ - mass)
			if tmp < best:
				best = tmp
				ix = idx
		return ix

	def crop_mass(self, mass_min: float, mass_max: float) -> None:
		"""
		Crops mass spectrum.

		:param mass_min: Minimum mass value
		:param mass_max: Maximum mass value

		:author: Andrew Isaac
		"""

		if not is_number(mass_min) or not is_number(mass_max):
			raise TypeError("'mass_min' and 'mass_max' must be numbers")
		if mass_min >= mass_max:
			raise ValueError("'mass_min' must be less than 'mass_max'")
		if mass_min < self._min_mass:
			raise ValueError(f"'mass_min' is less than the smallest mass: {self._min_mass:.3f}")
		if mass_max > self._max_mass:
			raise ValueError(f"'mass_max' is greater than the largest mass: {self._max_mass:.3f}")

		# pre build mass_list and list of indecies
		mass_list = self._mass_list
		new_mass_list = []
		ii_list = []
		for ii, mass in enumerate(mass_list):
			if mass_min <= mass <= mass_max:
				new_mass_list.append(mass)
				ii_list.append(ii)

		# update intensity matrix
		im: List[List[float]] = self._intensity_array.tolist()
		for spec_jj in range(len(im)):
			new_spec = []
			for ii in ii_list:
				new_spec.append(im[spec_jj][ii])
			im[spec_jj] = new_spec
		self._intensity_array = numpy.array(im)

		self._mass_list = new_mass_list
		self._min_mass = min(new_mass_list)
		self._max_mass = max(new_mass_list)

	def null_mass(self, mass: float) -> None:
		"""
		Ignore given (closest) mass in spectra.

		:param mass: Mass value to remove

		:author: Andrew Isaac
		"""

		if not is_number(mass):
			raise TypeError("'mass' must be a number")
		if mass < self._min_mass or mass > self._max_mass:
			raise IndexError(f"'mass' not in mass range: {self._min_mass:.3f} to {self._max_mass:.3f}")

		ii = self.get_index_of_mass(mass)

		im = self._intensity_array
		for spec_jj in range(len(im)):
			im[spec_jj][ii] = 0

	def reduce_mass_spectra(self, n_intensities: int = 5) -> None:
		"""
		Reduces the mass spectra by retaining the top `n_intensities`,
		discarding all other intensities.

		:param n_intensities: The number of top intensities to keep

		:author: Vladimir Likic
		"""  # noqa: D400

		if not is_number(n_intensities):
			raise TypeError("'n_intensities' must be a number")

		# loop over all mass spectral scans
		for ii, intensity_list in enumerate(self._intensity_array):

			# get the next mass spectrum as list of intensities
			# intensity_list = self._intensity_array[ii]
			n = len(intensity_list)

			# get the indices of top N intensities
			top_indices = list(range(n))
			top_indices.sort(key=lambda i: intensity_list[i], reverse=True)
			top_indices = top_indices[:n_intensities]

			# initiate new mass spectrum, and retain only top N intensities
			intensity_list_new = []

			for jj in range(n):
				intensity_list_new.append(0.0)
				if jj in top_indices:
					intensity_list_new[jj] = intensity_list[jj]

			self._intensity_array[ii] = intensity_list_new


class IntensityMatrix(BaseIntensityMatrix):
	"""
	Intensity matrix of binned raw data.

	:param time_list: Retention time values
	:param mass_list: Binned mass values
	:param intensity_array: List of lists of binned intensity values per scan

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions and properties)
	"""

	def __init__(
			self,
			time_list: Sequence[float],
			mass_list: Sequence[float],
			intensity_array: Union[Sequence[Sequence[float]], numpy.ndarray],
			):
		super().__init__(time_list, mass_list, intensity_array)

		# Try to include parallelism.
		try:
			# 3rd party
			from mpi4py import MPI  # type: ignore[import]
			comm = MPI.COMM_WORLD
			num_ranks = comm.Get_size()
			rank = comm.Get_rank()
			M, N = len(intensity_array), len(intensity_array[0])
			lrr = (rank * M / num_ranks, (rank + 1) * M / num_ranks)
			lcr = (rank * N / num_ranks, (rank + 1) * N / num_ranks)
			m, n = (lrr[1] - lrr[0], lcr[1] - lcr[0])
			self.comm = comm
			self.num_ranks = num_ranks
			self.rank = rank
			self.M = M
			self.N = N
			self.local_row_range = lrr
			self.local_col_range = lcr
			self.m = m
			self.n = n

		# If we can't import mpi4py then continue in serial.
		except ModuleNotFoundError:
			pass

	@property
	def local_size(self) -> Tuple[int, int]:
		"""
		Gets the local size of intensity matrix.

		:return: Number of rows and cols
		:rtype: int

		:author: Luke Hodkinson
		"""

		# Check for parallel.
		if hasattr(self, "comm"):
			return self.m, self.n

		# If serial call the regular routine.
		return self.size

	def iter_ms_indices(self) -> Iterator[int]:
		"""
		Iterates over the local row indices.

		:author: Luke Hodkinson
		"""

		# Check for parallel.
		if hasattr(self, "comm"):
			# At the moment we assume we break the matrix into contiguous
			# ranges. We've allowed for this to change by wrapping up the
			# iteration in this method.
			yield from range(self.local_row_range[0], self.local_row_range[1])

		else:
			# Iterate over global indices.
			return super().iter_ms_indices()

	def iter_ic_indices(self) -> Iterator[int]:
		"""
		Iterate over local column indices.

		:author: Luke Hodkinson
		"""

		# Check for parallel.
		if hasattr(self, "comm"):
			# At the moment we assume we break the matrix into contiguous
			# ranges. We've allowed for this to change by wrapping up the
			# iteration in this method.
			yield from range(self.local_col_range[0], self.local_col_range[1])

		else:
			# Iterate over global indices.
			return super().iter_ic_indices()

	def get_ic_at_mass(self, mass: Optional[float] = None) -> IonChromatogram:
		"""
		Returns the ion chromatogram for the nearest binned mass to the specified mass.

		If no mass value is given, the function returns the total ion chromatogram.

		:param mass: Mass value of an ion chromatogram

		:return: Ion chromatogram for given mass

		:authors: Andrew Isaac, Vladimir Likic
		"""

		if mass is None:
			return self.tic
		elif not is_number(mass):
			raise TypeError("'mass' must be a number")

		if mass < self._min_mass or mass > self._max_mass:
			print("min mass: ", self._min_mass, "max mass:", self._max_mass)
			raise IndexError("mass is out of range")

		ix = self.get_index_of_mass(mass)

		return self.get_ic_at_index(ix)

	@property
	def tic(self) -> IonChromatogram:
		"""
		Returns the TIC of the intensity matrix.

		.. versionadded:: 2.3.0
		"""

		intensity_list = []
		for i in range(len(self._intensity_array)):
			intensity_list.append(sum(self._intensity_array[i]))

		return IonChromatogram(numpy.array(intensity_list), self._time_list[:], None)

	def export_ascii(
			self,
			root_name: PathLike,
			fmt: AsciiFiletypes = AsciiFiletypes.ASCII_DAT,
			) -> None:
		"""
		Exports the intensity matrix, retention time vector, and m/z vector to the ascii format.

		By default, export_ascii("NAME") will create NAME.im.dat, NAME.rt.dat,
		and NAME.mz.dat where these are the intensity matrix, retention time
		vector, and m/z vector in tab delimited format.

		If ``format`` == ``<AsciiFiletypes.ASCII_CSV>``, the files will be in the CSV format, named
		NAME.im.csv, NAME.rt.csv, and NAME.mz.csv.

		:param root_name: Root name for the output files
		:param fmt: Format of the output file, either ``<AsciiFiletypes.ASCII_DAT>`` or ``<AsciiFiletypes.ASCII_CSV>``

		:authors: Milica Ng, Andrew Isaac, Vladimir Likic, Dominic Davis-Foster (pathlib support)
		"""

		if not is_path(root_name):
			raise TypeError("'root_name' must be a string or a pathlib.Path object")

		root_name = prepare_filepath(root_name, mkdirs=True)
		fmt = AsciiFiletypes(fmt)

		if fmt is AsciiFiletypes.ASCII_DAT:
			separator = ' '
			extension = "dat"
		elif fmt is AsciiFiletypes.ASCII_CSV:
			separator = ','
			extension = "csv"

		# export 2D matrix of intensities
		vals = self._intensity_array

		# save_data is annotated as taking only lists but list-like (i.e. numpy array) is ok
		save_data(
				f"{root_name}.im.{extension}",
				vals,  # type: ignore[arg-type]
				sep=separator,
				)

		# export 1D vector of m/z's, corresponding to rows of
		# the intensity matrix
		mass_list = self._mass_list
		save_data(f"{root_name}.mz.{extension}", mass_list, sep=separator)

		# export 1D vector of retention times, corresponding to
		# columns of the intensity matrix
		time_list = self._time_list
		save_data(f"{root_name}.rt.{extension}", time_list, sep=separator)

	def export_leco_csv(self, file_name: PathLike) -> None:
		"""
		Exports data in LECO CSV format.

		:param file_name: The name of the output file.

		:authors: Andrew Isaac, Vladimir Likic, Dominic Davis-Foster (pathlib support)
		"""

		if not is_path(file_name):
			raise TypeError("'file_name' must be a string or a PathLike object")

		file_name = prepare_filepath(file_name, mkdirs=False)

		if not file_name.parent.is_dir():
			file_name.parent.mkdir(parents=True)

		mass_list = self._mass_list
		time_list = self._time_list
		vals = self._intensity_array

		with file_name.open('w', encoding="UTF-8") as fp:

			# Format is text header with:
			# "Scan","Time",...
			# and the rest is "TIC" or m/z as text, i.e. "50","51"...
			# The following lines are:
			# scan_number,time,value,value,...
			# scan_number is an int, rest seem to be fixed format floats.
			# The format is 0.000000e+000

			# write header
			fp.write('"Scan","Time"')
			for ii in mass_list:
				if is_number(ii):
					fp.write(f',"{int(ii):d}"')
				else:
					raise TypeError("mass list datum not a number")
			fp.write("\r\n")  # windows CR/LF

			# write lines
			for ii, time_ in enumerate(time_list):
				fp.write(f"{ii},{time_:#.6e}")
				for jj in range(len(vals[ii])):
					if is_number(vals[ii][jj]):
						fp.write(f",{vals[ii][jj]:#.6e}")
					else:
						raise TypeError("datum not a number")
				fp.write("\r\n")

	@property
	def bpc(self) -> IonChromatogram:
		"""
		Constructs a Base Peak Chromatogram from the data.

		This represents the most intense ion for each scan.

		:authors: Dominic Davis-Foster

		.. versionadded:: 2.3.0
		"""

		return BasePeakChromatogram(
				[max(intensities) for intensities in self._intensity_array],
				self._time_list[:],
				)


def import_leco_csv(file_name: PathLike) -> IntensityMatrix:
	"""
	Imports data in LECO CSV format.

	:param file_name: Path of the file to read.

	:return: Data as an IntensityMatrix.

	:authors: Andrew Isaac, Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	data = []
	time_list = []
	mass_list = []

	# Format is text header with:
	# "Scan","Time",...
	# and the rest is "TIC" or m/z as text, i.e. "50","51"...
	# The following lines are:
	# scan_number,time,value,value,...
	# scan_number is an int, rest seem to be fixed format floats.
	# The format is 0.000000e+000

	num_mass = 0
	FIRST = True
	HEADER = True
	data_col = -1
	time_col = -1

	with file_name.open('r', encoding="UTF-8") as lines_list:

		# get each line
		for line in lines_list:
			cols = -1
			data_row = []
			if len(line.strip()) > 0:
				data_list = line.strip().split(',')
				# get each value in line
				for item in data_list:
					item = item.strip()
					item = item.strip("'\"")  # remove quotes (in header)

					# Get header
					if HEADER:
						cols += 1
						if len(item) > 0:
							if item.lower().find("time") > -1:
								time_col = cols
							try:
								value = float(item)
								# find 1st col with number as header
								if FIRST and value > 1:  # assume >1 mass
									data_col = cols
									# assume time col is previous col
									if time_col < 0:
										time_col = cols - 1
									FIRST = False
								mass_list.append(value)
								num_mass += 1
							except ValueError:
								pass
					# Get rest
					else:
						cols += 1
						if len(item) > 0:
							try:
								value = float(item)
								if cols == time_col:
									time_list.append(value)
								elif cols >= data_col:
									data_row.append(value)
							except ValueError:
								pass

				# check row length
				if not HEADER:
					if len(data_row) == num_mass:
						data.append(data_row)
					else:
						warn("ignoring row")

				HEADER = False

	# check col lengths
	if len(time_list) != len(data):
		warn("number of data rows and time list length differ")

	return IntensityMatrix(time_list, mass_list, data)


def build_intensity_matrix(
		data: GCMS_data,
		bin_interval: float = 1,
		bin_left: float = 0.5,
		bin_right: float = 0.5,
		min_mass: Optional[float] = None,
		) -> IntensityMatrix:
	"""
	Sets the full intensity matrix with flexible bins.

	The first bin is centered around ``min_mass``,
	and subsequent bins are offset by ``bin_interval``.

	:param data: Raw GCMS data
	:param bin_interval: interval between bin centres.
	:param bin_left: left bin boundary offset.
	:param bin_right: right bin boundary offset.
	:param min_mass: Minimum mass to bin (default minimum mass from data)

	:return: Binned IntensityMatrix object

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
	"""

	# this package
	from pyms.GCMS.Class import GCMS_data

	if not isinstance(data, GCMS_data):
		raise TypeError("'data' must be a GCMS_data object")

	if bin_interval <= 0:
		raise ValueError("The bin interval must be larger than zero.")

	if not is_number(bin_left):
		raise TypeError("'bin_left' must be a number.")

	if not is_number(bin_right):
		raise TypeError("'bin_right' must be a number.")

	if min_mass is None:
		min_mass = data.min_mass
	elif not is_number(min_mass):
		raise TypeError("'min_mass' must be a number.")

	max_mass = data.max_mass

	if max_mass is None:
		raise ValueError("'max_mass' cannot be None")
	if min_mass is None:
		raise ValueError("'min_mass' cannot be None")

	return _fill_bins(data, min_mass, max_mass, bin_interval, bin_left, bin_right)


def build_intensity_matrix_i(
		data: GCMS_data,
		bin_left: float = 0.3,
		bin_right: float = 0.7,
		) -> IntensityMatrix:
	"""
	Sets the full intensity matrix with integer bins.

	:param data: Raw GCMS data
	:param bin_left: left bin boundary offset.
	:param bin_right: right bin boundary offset.

	:return: Binned IntensityMatrix object

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
	"""

	# this package
	from pyms.GCMS.Class import GCMS_data

	if not isinstance(data, GCMS_data):
		raise TypeError("'data' must be a GCMS_data object")

	if not is_number(bin_left):
		raise TypeError("'bin_left' must be a number.")

	if not is_number(bin_right):
		raise TypeError("'bin_right' must be a number.")

	min_mass = data.min_mass
	max_mass = data.max_mass

	if max_mass is None:
		raise ValueError("'max_mass' cannot be None")
	if min_mass is None:
		raise ValueError("'min_mass' cannot be None")

	# Calculate integer min mass based on right boundary
	bin_right = abs(bin_right)
	min_mass = int(min_mass + 1 - bin_right)

	return _fill_bins(data, min_mass, max_mass, 1, bin_left, bin_right)


def _fill_bins(
		data: GCMS_data,
		min_mass: float,
		max_mass: float,
		bin_interval: Union[int, float],
		bin_left: float,
		bin_right: float,
		) -> IntensityMatrix:
	"""
	Fills the intensity values for all bins.

	:param data: Raw GCMS data
	:param min_mass: minimum mass value
	:param max_mass: maximum mass value
	:param bin_interval: interval between bin centres
	:param bin_left: left bin boundary offset
	:param bin_right: right bin boundary offset

	:return: Binned IntensityMatrix object

	:authors: Qiao Wang, Andrew Isaac, Moshe Olshansky, Vladimir Likic
	"""

	if not (abs(bin_left + bin_right - bin_interval) < 1.0e-6 * bin_interval):
		raise ValueError("there should be no gaps or overlap between the bins.")

	bin_left = abs(bin_left)
	# bin_right = abs(bin_right)

	# To convert to int range, ensure bounds are < 1
	bl = bin_left - int(bin_left)

	# Number of bins
	num_bins = int(float(max_mass + bl - min_mass) / bin_interval) + 1

	# initialise masses to bin centres
	mass_list = [i * bin_interval + min_mass for i in range(num_bins)]

	# Modified binning loops. I've replaced the deepcopy getting routines with
	# the alias properties. This way we can avoid performing the copies when
	# it is clear that we do not intend on modifying the contents of the arrays
	# here.
	#           - Luke Hodkinson, 18/05/2010

	# fill the bins
	intensity_matrix = []
	for scan in data.scan_list:  # use the alias, not the copy (Luke)
		intensity_list = [0.0] * num_bins
		masses = scan.mass_list  # use the alias, not the copy (Luke)
		intensities = scan.intensity_list  # use the alias, not the copy (Luke)
		for ii, mass in enumerate(masses):
			mm = int((mass + bl - min_mass) / bin_interval)
			intensity_list[mm] += intensities[ii]
		intensity_matrix.append(intensity_list)

	return IntensityMatrix(data.time_list, mass_list, intensity_matrix)


def _fill_bins_old(
		data: GCMS_data,
		min_mass: float,
		max_mass: float,
		bin_interval: float,
		bin_left: float,
		bin_right: float,
		) -> IntensityMatrix:
	"""
	Fills the intensity values for all bins.

	:param data: Raw GCMS data
	:param min_mass: minimum mass value
	:param max_mass: maximum mass value
	:param bin_interval: interval between bin centres
	:param bin_left: left bin boundary offset
	:param bin_right: right bin boundary offset

	:return: Binned IntensityMatrix object

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic
	"""

	bin_left = abs(bin_left)
	bin_right = abs(bin_right)

	# To convert to int range, ensure bounds are < 1
	bl = bin_left - int(bin_left)

	# Number of bins
	num_bins = int(float(max_mass + bl - min_mass) / bin_interval) + 1

	# initialise masses to bin centres
	mass_list = [i * bin_interval + min_mass for i in range(num_bins)]

	# fill the bins
	intensity_matrix = []
	for scan in data.scan_list:
		intensity_list = [0.0] * num_bins
		masses = scan.mass_list
		intensities = scan.intensity_list
		for mm in range(num_bins):
			for ii in range(len(scan)):
				if mass_list[mm] - bin_left <= masses[ii] < mass_list[mm] + bin_right:
					intensity_list[mm] += intensities[ii]
		intensity_matrix.append(intensity_list)

	return IntensityMatrix(data.time_list, mass_list, intensity_matrix)

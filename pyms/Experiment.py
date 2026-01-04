"""
Models a GC-MS experiment, represented by a list of signal peaks.
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
from typing import List, Sequence

# 3rd party
from domdf_python_tools.typing import PathLike

# this package
from pyms.Base import pymsBaseClass
from pyms.Peak.Class import Peak
from pyms.Peak.List.Function import is_peak_list, sele_peaks_by_rt
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Utils import _pickle_load_path, is_path, is_sequence

__all__ = ["Experiment", "read_expr_list", "load_expr"]


class Experiment(pymsBaseClass):
	"""
	Models an experiment.

	:param expr_code: A unique identifier for the experiment.
	:param peak_list:

	:author: Vladimir Likic, Andrew Isaac,  Dominic Davis-Foster (type assertions, properties and pathlib support)
	"""

	def __init__(self, expr_code: str, peak_list: Sequence[Peak]):
		if not isinstance(expr_code, str):
			raise TypeError("'expr_code' must be a string")

		if not is_peak_list(peak_list):
			raise TypeError("'peak_list' must be a list of Peak objects")

		self._expr_code = expr_code
		self._peak_list = list(peak_list)

	def __eq__(self, other) -> bool:  # noqa: MAN001
		"""
		Return whether this Experiment object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return self.peak_list == other.peak_list and self.expr_code == other.expr_code

		return NotImplemented

	def __len__(self) -> int:
		"""
		Returns the number of peaks in the Experiment.
		"""

		return len(self.peak_list)

	def __copy__(self) -> "Experiment":
		"""
		Returns a new Experiment object containing a copy of the data in this object.
		"""

		return Experiment(copy.copy(self._expr_code), copy.copy(self.peak_list))

	def __deepcopy__(self, memodict={}) -> "Experiment":  # noqa: MAN001,PRM002
		"""
		Returns a new Experiment object containing a copy of the data in this object.
		"""

		return self.__copy__()

	@property
	def expr_code(self) -> str:
		"""
		Returns the expr_code of the experiment.
		"""

		return self._expr_code

	@property
	def peak_list(self) -> List[Peak]:
		"""
		Returns the peak list.
		"""

		return self._peak_list

	def sele_rt_range(self, rt_range: Sequence[str]) -> None:
		"""
		Discards all peaks which have the retention time outside the specified range.

		:param rt_range: Min, max retention time given as a sequence ``[rt_min, rt_max]``.
		"""

		if not is_sequence(rt_range):
			raise TypeError("'rt_range' must be a Sequence")

		peaks_sele = sele_peaks_by_rt(self._peak_list, rt_range)
		self._peak_list = peaks_sele


def read_expr_list(file_name: PathLike) -> List[Experiment]:
	"""
	Reads the set of experiment files and returns a list of :class:`pyms.Experiment.Experiment` objects.

	:param file_name: The name of the file which lists experiment dump file names, one file per line.

	:return: A list of Experiment instances.

	:author: Vladimir Likic
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	with file_name.open(encoding="UTF-8") as fp:
		exprfiles = fp.readlines()

	expr_list = []

	for exprfile in exprfiles:

		exprfile = exprfile.strip()
		expr = load_expr(exprfile)

		expr_list.append(expr)

	return expr_list


def load_expr(file_name: PathLike) -> Experiment:
	"""
	Loads an experiment saved with :meth:`pyms.Experiment.Experiment.dump`.

	:param file_name: Experiment file name.

	:return: The loaded experiment.

	:author: Vladimir Likic, Andrew Isaac, Dominic Davis-Foster (type assertions and pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)
	expr = _pickle_load_path(file_name)

	if not isinstance(expr, Experiment):
		raise TypeError("The loaded file is not an experiment file")

	return expr

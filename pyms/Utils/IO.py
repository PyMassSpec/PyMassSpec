"""
General I/O functions.
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
import gzip
import pathlib
import pickle
from typing import Any, List, Union, cast

# 3rd party
from domdf_python_tools.stringlist import StringList
from domdf_python_tools.typing import PathLike

# this package
from pyms.Utils.Utils import _list_types, _pickle_dump_path, is_number, is_path

__all__ = ["prepare_filepath", "dump_object", "load_object", "file_lines", "save_data"]


def prepare_filepath(
		file_name: PathLike,
		mkdirs: bool = True,
		) -> pathlib.Path:
	"""
	Convert string filename into pathlib.Path object and
	create parent directories if required.

	:param file_name: file_name to process
	:param mkdirs: Whether the parent directory of the file should
		be created if it doesn't exist.

	:return: file_name

	:author: Dominic Davis-Foster
	"""  # noqa: D400

	if not isinstance(file_name, pathlib.Path):
		try:
			file_name = pathlib.Path(file_name)
		except TypeError:
			raise TypeError(f"'file_name' must be a string or a PathLike object, not {type(file_name)}")

	if not file_name.parent.is_dir() and mkdirs:
		file_name.parent.mkdir(parents=True)

	return file_name


def dump_object(obj: Any, file_name: PathLike):
	"""
	Dumps an object to a file through :func:`pickle.dump`.

	:param obj: Object to be dumped
	:param file_name: Name of the file for the object dump

	:authors: Vladimir Likic, Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	_pickle_dump_path(prepare_filepath(file_name), obj)


def load_object(file_name: PathLike) -> object:
	"""
	Loads an object previously dumped with :func:`~.dump_object`.

	:param file_name: Name of the object dump file.

	:return: Object contained in the file.

	:authors: Vladimir Likic, Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name)

	with file_name.open("wb") as fp:
		return pickle.load(fp)


def file_lines(file_name: PathLike, strip: bool = False) -> List[str]:
	"""
	Returns lines from a file, as a list.

	:param file_name: Name of a file
	:param strip: If True, lines are pre-processed. Newline characters are
		removed, leading and trailing whitespaces are removed, and lines
		starting with '#' are discarded

	:return: A list of lines

	:authors: Vladimir Likic, Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	with file_name.open(encoding="UTF-8") as fp:
		lines = fp.readlines()

	if strip:
		# strip leading and trailing whitespaces
		lines_filtered = []
		for line in lines:
			line = line.strip()
			lines_filtered.append(line)

		# discard comments
		lines_to_discard = []
		for line in lines_filtered:
			# remove empty lines and comments
			if len(line) == 0 or line[0] == '#':
				lines_to_discard.append(line)
		for line in lines_to_discard:
			lines_filtered.remove(line)
		lines = lines_filtered

	return lines


def save_data(
		file_name: PathLike,
		data: Union[List[float], List[List[float]]],
		format_str: str = "%.6f",
		prepend: str = '',
		sep: str = ' ',
		compressed: bool = False,
		):
	"""
	Saves a list of numbers or a list of lists of numbers to a file with specific formatting.

	:param file_name: Name of a file
	:param data: A list of numbers, or a list of lists
	:param format_str: A format string for individual entries
	:param prepend: A string, printed before each row
	:param sep: A string, printed after each number
	:param compressed: If :py:obj:`True`, the output will be gzipped.

	:authors: Vladimir Likic, Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name)

	if not isinstance(data, _list_types):
		raise TypeError("'data' must be a list")

	if not isinstance(prepend, str):
		raise TypeError("'prepend' must be a string")

	if not isinstance(sep, str):
		raise TypeError("'sep' must be a string")

	buf = StringList()

	# decide whether data is a vector or matrix
	if is_number(data[0]):
		for item in data:
			if not is_number(item):
				raise TypeError("not all elements of the list are numbers")
		for x_value in data:
			buf.append(prepend + (format_str % x_value))

	else:
		for item in data:
			if not isinstance(item, _list_types):
				raise TypeError("not all elements of the list are lists")

		for x_value in cast(List[List[float]], data):
			line = [prepend]
			for jj, y_value in enumerate(x_value):
				if is_number(y_value):
					line.append(format_str % y_value)
					if jj < (len(x_value) - 1):
						line.append(sep)
				else:
					raise TypeError("'datum' must be a number")
			buf.append(''.join(line))

	if compressed:
		with gzip.open(file_name, "wt") as fp:
			fp.write(str(buf))
	else:
		file_name.write_text(str(buf), encoding="UTF-8")

"""
Functions related to storing and loading a list of Peak objects
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
import pickle
from typing import List, Sequence, Union

# this package
from pyms.Peak.Class import Peak
from pyms.Peak.List.Function import is_peak_list
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Utils import is_path, is_sequence

__all__ = ["store_peaks", "load_peaks"]


def store_peaks(
		peak_list: Sequence[Peak],
		file_name: Union[str, pathlib.Path],
		protocol: int = 1,
		):
	"""
	Store the list of peak objects

	:param peak_list: A list of peak objects
	:param file_name: File name to store peak list
	:param protocol: The :mod:`pickle` protocol to use.

	:authors: Andrew Isaac, Dominic Davis-Foster (type assertions and pathlib support)
	"""

	if not is_peak_list(peak_list):
		raise TypeError("'peak_list' must be a list of Peak objects")

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name)

	fp = file_name.open('wb')
	pickle.dump(peak_list, fp, protocol)
	fp.close()


def load_peaks(file_name: Union[str, pathlib.Path]) -> List[Peak]:
	"""
	Loads the peak_list stored with :func:`~.store_peaks`.

	:param file_name: File name of peak list

	:return: The list of Peak objects

	:authors: Andrew Isaac, Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	fp = file_name.open('rb')
	peak_list = pickle.load(fp)
	fp.close()

	if not is_sequence(peak_list):
		raise IOError("The selected file is not a List")
	if not len(peak_list) > 0 or not isinstance(peak_list[0], Peak):
		raise IOError("The selected file is not a list of Peak objects")

	return peak_list

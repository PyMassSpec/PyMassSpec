"""
Functions for I/O of data in JCAMP-DX format.
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                              #
#                                                                              #
#    Parts based on 'jcamp' by Nathan Hagen									   #
# 	 https://github.com/nzhagen/jcamp										   #
# 	 Licensed under the X11 License											   #
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
from pathlib import Path
from typing import Any, List, MutableMapping, Union

# 3rd party
from domdf_python_tools.paths import PathPlus

# this package
from pyms.GCMS.Class import GCMS_data
from pyms.Spectrum import Scan
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.jcamp import header_info_fields, xydata_tags
from pyms.Utils.Math import is_float
from pyms.Utils.Utils import is_path

__all__ = ["JCAMP_reader"]


def JCAMP_reader(file_name: Union[str, Path]) -> GCMS_data:
	"""
	Generic reader for JCAMP DX files.

	:param file_name: Path of the file to read

	:return: GC-MS data object

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic, David Kainer,
		Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = PathPlus(prepare_filepath(file_name, mkdirs=False))

	print(f" -> Reading JCAMP file {file_name.as_posix()!r}")
	lines_list = file_name.read_lines()
	time_list: List[float] = []
	scan_list: List[Scan] = []

	header_info: MutableMapping[Any, Any] = {}  # Dictionary containing header information

	line_idx = 0
	while line_idx < len(lines_list):
		line = lines_list[line_idx]

		if line.strip():
			if line.startswith("##"):
				# Label
				label, value = line.split('=', 1)
				label = label.lstrip("##").upper()
				value = value.strip()

				if "PAGE" in label:
					if "T=" in value:
						# PAGE contains retention time starting with T=
						# FileConverter Pro style
						time = float(value.lstrip("T="))  # rt for the scan to be submitted
						time_list.append(time)

				elif "RETENTION_TIME" in label:
					# OpenChrom style
					time = float(value)  # rt for the scan to be submitted

					# Check to make sure time is not already in the time list;
					# Can happen when both ##PAGE and ##RETENTION_TIME are specified
					if time_list[-1] != time:
						time_list.append(time)

				elif label in header_info_fields:
					if value.isdigit():
						header_info[label] = int(value)
					elif is_float(value):
						header_info[label] = float(value)
					else:
						header_info[label] = value

				elif label in xydata_tags:
					# Read ahead to find all XY data
					xydata_line_idx = line_idx + 1
					xy_data_lines: List[str] = []
					while xydata_line_idx < len(lines_list):
						xy_data_line = lines_list[xydata_line_idx]
						if xy_data_line.startswith("##"):
							break
						else:
							xy_data_lines.append(xy_data_line)
							xydata_line_idx += 1

					line_idx += len(xy_data_lines)

					mass_list = []
					intensity_list = []
					for xy_data_line in xy_data_lines:
						if not xy_data_line.strip():
							continue

						elements = xy_data_line.strip().rstrip(',').split(',')
						if len(elements) % 2:
							print(elements)
							raise ValueError(f"Expected an even number of values, got {len(elements)}")
						for mass, intensity in zip(elements[::2], elements[1::2]):
							mass_list.append(float(mass.strip()))
							intensity_list.append(float(intensity.strip()))

					scan_list.append(Scan(mass_list, intensity_list))

		line_idx += 1  # End of while loop

	# sanity check
	time_len = len(time_list)
	scan_len = len(scan_list)
	if time_len != scan_len:
		print(time_list)
		print(scan_list)
		raise ValueError(f"Number of time points ({time_len}) does not equal the number of scans ({scan_len})")

	return GCMS_data(time_list, scan_list)

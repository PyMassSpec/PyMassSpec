"""
Functions for I/O of data in JCAMP-DX format
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
import pathlib

# this package
from typing import Union

from pyms.GCMS.Class import GCMS_data
from pyms.Spectrum import Scan
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.jcamp import header_info_fields, xydata_tags
from pyms.Utils.Math import is_float
from pyms.Utils.Utils import is_path


def JCAMP_reader(file_name: Union[str, pathlib.Path] -> GCMS_data):
	"""
	Generic reader for JCAMP DX files

	:param file_name: Path of the file to read
	:type file_name: str or pathlib.Path

	:return: GC-MS data object
	:rtype: :class:`pyms.GCMS.Class.GCMS_data`

	:authors: Qiao Wang, Andrew Isaac, Vladimir Likic, David Kainer, Dominic Davis-Foster (pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	print(f" -> Reading JCAMP file '{file_name}'")
	lines_list = file_name.open('r')
	data = []
	page_idx = 0
	xydata_idx = 0
	time_list = []
	scan_list = []

	header_info = {}  # Dictionary containing header information

	for line in lines_list:

		if len(line.strip()) != 0:
			# prefix = line.find('#')
			# if prefix == 0:
			if line.startswith("##"):
				# key word or information
				fields = line.split('=', 1)
				fields[0] = fields[0].lstrip("##").upper()
				fields[1] = fields[1].strip()

				if "PAGE" in fields[0]:
					if "T=" in fields[1]:
						# PAGE contains retention time starting with T=
						# FileConverter Pro style
						time = float(fields[1].lstrip("T="))  # rt for the scan to be submitted
						time_list.append(time)
					page_idx = page_idx + 1
				elif "RETENTION_TIME" in fields[0]:
					# OpenChrom style
					time = float(fields[1])  # rt for the scan to be submitted

					# Check to make sure time is not already in the time list;
					# Can happen when both ##PAGE and ##RETENTION_TIME are specified
					if time_list[-1] != time:
						time_list.append(time)

				elif fields[0] in xydata_tags:
					xydata_idx = xydata_idx + 1

				elif fields[0] in header_info_fields:
					if fields[1].isdigit():
						header_info[fields[0]] = int(fields[1])
					elif is_float(fields[1]):
						header_info[fields[0]] = float(fields[1])
					else:
						header_info[fields[0]] = fields[1]

			# elif prefix == -1:
			else:
				# Line doesn't start with ##
				# data
				if page_idx > 1 or xydata_idx > 1:
					if len(data) % 2 == 1:
						# TODO: This means the data is not in x, y pairs
						#  Make a better error message
						raise ValueError("data not in pair !")
					mass_list = []
					intensity_list = []
					for i in range(len(data) // 2):
						mass_list.append(data[i * 2])
						intensity_list.append(data[i * 2 + 1])
					if len(mass_list) != len(intensity_list):
						raise ValueError("len(mass_list) is not equal to len(intensity_list)")
					scan_list.append(Scan(mass_list, intensity_list))
					data = []
					data_sub = line.strip().split(',')
					for item in data_sub:
						if not len(item.strip()) == 0:
							data.append(float(item.strip()))
					if page_idx > 1:
						page_idx = 1
					if xydata_idx > 1:
						xydata_idx = 1
				else:
					data_sub = line.strip().split(',')
					for item in data_sub:
						if not len(item.strip()) == 0:
							data.append(float(item.strip()))

	if len(data) % 2 == 1:
		# TODO: This means the data is not in x, y pairs
		#  Make a better error message
		raise ValueError("data not in pair !")

	# get last scan
	mass = []
	intensity = []
	for i in range(len(data) // 2):
		mass.append(data[i * 2])
		intensity.append(data[i * 2 + 1])

	if len(mass) != len(intensity):
		raise ValueError("len(mass) is not equal to len(intensity)")
	scan_list.append(Scan(mass, intensity))

	# sanity check
	time_len = len(time_list)
	scan_len = len(scan_list)
	if time_len != scan_len:
		print(time_list)
		print(scan_list)
		raise ValueError(f"Number of time points ({time_len}) does not equal the number of scans ({scan_len})")

	data = GCMS_data(time_list, scan_list)

	return data

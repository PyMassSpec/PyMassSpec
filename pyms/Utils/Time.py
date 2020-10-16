"""
Time conversion and related functions
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
import math
import re
from typing import Union

# this package
from pyms.IonChromatogram import IonChromatogram

__all__ = ["is_str_num", "time_str_secs", "window_sele_points"]

num_re = re.compile(r'^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$')


def is_str_num(arg: str) -> bool:
	"""
	Returns whether the argument is a string in the format of a number.

	The number can be an integer, or alternatively a floating point number in scientific
	or engineering format.

	:param arg: A string to be evaluate as a number

	:author: Gyro Funch (from Active State Python Cookbook)
	"""

	if num_re.match(str(arg)):
		return True
	return False


def time_str_secs(time_str: str) -> float:
	"""
	Resolves time string of the form ``'<NUMBER>s'`` or ``'<NUMBER>m'`` and returns the time in seconds

	:param time_str: A time string, which must be of the form
		``'<NUMBER>s'`` or ``'<NUMBER>m'`` where ``'<NUMBER>'`` is a valid number

	:return: Time in seconds

	:author: Vladimir Likic
	"""

	if not isinstance(time_str, str):
		raise TypeError("'time_str' must be a string")

	time_number = time_str[:-1]
	time_spec = time_str[-1].lower()

	if not is_str_num(time_number):
		print(f" --> received time string '{time_number}'")
		raise ValueError("improper time string")

	if time_spec not in ["s", "m"]:
		raise ValueError("time string must end with either 's' or 'm'")

	time = float(time_number)

	if time_spec == "m":
		time = time * 60.0

	return time


def window_sele_points(ic: IonChromatogram, window_sele: Union[int, str], half_window: bool = False) -> int:
	"""
	Converts window selection parameter into points based
	on the time step in an ion chromatogram

	:param ic: ion chromatogram object relevant for the conversion

	:param window_sele: The window selection parameter. This can be
		an integer or time string. If an integer, taken as the number
		of points. If a string, must of the form ``'<NUMBER>s'`` or
		``'<NUMBER>m'``, specifying a time in seconds or minutes,
		respectively

	:param half_window: Specifies whether to return half-window

	:return: The number of points in the window

	:author: Vladimir Likic
	"""

	if not isinstance(window_sele, int) and not isinstance(window_sele, str):
		raise TypeError("'window' must be an integer or a string")

	if isinstance(window_sele, int):
		if half_window:
			if window_sele % 2 == 0:
				raise TypeError("window must be an odd number of points")
			else:
				points = int(math.floor(window_sele * 0.5))
		else:
			points = window_sele
	else:
		time = time_str_secs(window_sele)
		time_step = ic.time_step

		if half_window:
			time = time * 0.5

		points = int(math.floor(time / time_step))

	if half_window:
		if points < 1:
			raise ValueError(f"window too small (half window={points:d})")
	else:
		if points < 2:
			raise ValueError(f"window too small (window={points:d})")

	return points

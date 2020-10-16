#############################################################################
#                                                                           #
#    PyMassSpec software for processing of mass-spectrometry data           #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                           #
#                                                                           #
#    This program is free software; you can redistribute it and/or modify   #
#    it under the terms of the GNU General Public License version 2 as      #
#    published by the Free Software Foundation.                             #
#                                                                           #
#    This program is distributed in the hope that it will be useful,        #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#    GNU General Public License for more details.                           #
#                                                                           #
#    You should have received a copy of the GNU General Public License      #
#    along with this program; if not, write to the Free Software            #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#                                                                           #
#############################################################################

__all__ = [
		"test_string",
		"test_int",
		"test_float",
		"test_list_ints",
		"test_list_strs",
		"test_dict",
		"test_numbers",
		"test_tuple",
		"test_lists",
		"test_sequences",
		]

test_string = "abcdefg"
test_int = 1234
test_float = 12.34
test_list_ints = [1, 2, 3, 4]
test_list_strs = ["a", "b", "c", "d"]
test_dict = {"a": 1, "b": 2, "c": 3, "d": 4}

test_numbers = {test_int, test_float}

test_tuple = (*test_numbers, test_string)

test_lists = [test_list_ints, test_list_strs]
test_sequences = [*test_lists, test_tuple]

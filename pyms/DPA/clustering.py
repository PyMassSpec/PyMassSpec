"""
Provides Pycluster.treecluster regardless of which library provides it.
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

try:
	# 3rd party
	from Pycluster import treecluster  # type: ignore[import]
except ModuleNotFoundError:
	try:
		# 3rd party
		from Bio.Cluster import treecluster
	except ModuleNotFoundError:
		raise ModuleNotFoundError(
				"""Neither PyCluster or BioPython is installed.
Please install one of them and try again."""
				) from None

__all__ = ["treecluster"]

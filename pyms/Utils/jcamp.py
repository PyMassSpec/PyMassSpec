"""
Constants required for reading JCAMP files.
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2020 Dominic Davis-Foster                                   #
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

__all__ = ["JcampTagWarning", "header_info_fields", "xydata_tags"]

header_info_fields = [
		"TITLE",
		"JCAMP-DX",
		"SAMPLE_DESCRIPTION",
		"DATE",
		"TIME",
		"SPECTROMETER_SYSTEM",
		"EXPERIMENT_NAME",
		"INLET",
		"IONIZATION_MODE",  # e.g. EI+
		"ELECTRON_ENERGY",  # e.g. 70
		"RESOLUTION",
		"ACCELERATING_VOLTAGE",  # e.g. 8000
		"CALIBRATION_FILE",
		"REFERENCE_FILE",
		"MASS_RANGE",
		"SCAN_LAW",  # e.g. Exponential
		"SCAN_RATE_UNITS",  # Units for SCAN_RATE
		"SCAN_RATE",  # Rate at which scans were acquired
		"SCAN_DELAY_UNITS",  # Units for SCAN_DELAY
		"SCAN_DELAY",  # Time delay in seconds before first scan acquired
		"XUNITS",  # Units for X-Axis e.g. Daltons
		"DATA_FORMAT",  # e.g. Centroid
		"DATA TYPE",  # e.g. MASS SPECTRUM
		"DATA CLASS",  # e.g. NTUPLES
		"ORIGIN",
		"OWNER",
		"CAS REGISTRY NO",
		"$NIST MASS SPEC NO",
		"MOLFORM",
		"MW",
		"$NIST SOURCE.",
		]

xydata_tags = {"XYDATA", "DATA TABLE", "XYPOINTS, PEAK TABLE", "PEAK TABLE", "XYPOINTS"}


class JcampTagWarning(UserWarning):
	"""
	Warning emitted when an unrecognised Jcamp tag is detected.

	:param tag: The name of the tag.
	"""

	#: The name of the tag.
	tag: str

	def __init__(self, tag: str):
		self.tag = str(tag)

	def __repr__(self) -> str:
		return f"JcampTagWarning(Unrecognised tag. tag={self.tag})"

	def __str__(self) -> str:
		return f"Unrecognised tag {self.tag}."

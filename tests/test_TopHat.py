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

# 3rd party
import pytest

# this package
from pyms.IntensityMatrix import IntensityMatrix
from pyms.IonChromatogram import IonChromatogram
from pyms.TopHat import tophat, tophat_im

# this package
from .constants import *


def test_topHat(tic):
	assert isinstance(tic, IonChromatogram)

	# apply noise smoothing and baseline correction
	tic2 = tophat(tic, struct="1.5m")
	assert isinstance(tic2, IonChromatogram)

	tic3 = tophat(tic, struct=None)
	assert isinstance(tic3, IonChromatogram)

	tic4 = tophat(tic, struct=1234)
	assert isinstance(tic4, IonChromatogram)


def test_tophat_im(im):
	# Use TopHat baseline correction on all IC's in the IM
	im_base_corr = tophat_im(im, struct="1.5m")
	assert isinstance(im_base_corr, IntensityMatrix)

	# find the IC for derivatisation product ion before smoothing
	ic = im.get_ic_at_index(73)
	assert isinstance(ic, IonChromatogram)

	# find the IC for derivatisation product ion after smoothing
	ic_base_corr = im_base_corr.get_ic_at_index(73)
	assert isinstance(ic_base_corr, IonChromatogram)


class TestErrors:

	@pytest.mark.parametrize("obj", [test_string, *test_numbers, *test_sequences])
	class TestobjErrors:

		def test_im_errors(self, obj):
			with pytest.raises(TypeError):
				tophat_im(obj, "1m")  # type: ignore

		def test_ic_errors(self, obj):
			with pytest.raises(TypeError):
				tophat(obj, "1m")  # type: ignore

	@pytest.mark.parametrize(
			"struct, expects", [
					(test_float, TypeError),
					(test_string, ValueError),
					] + [(struct, TypeError) for struct in test_sequences]
			)
	class TeststructErrors:

		def test_im_errors(self, im, struct, expects):
			with pytest.raises(expects):
				tophat_im(im, struct)

		def test_ic_errors(self, tic, struct, expects):
			with pytest.raises(expects):
				tophat(tic, struct)


# TODO:
# ic.write("output/ic.dat",minutes=True)
# ic_base_corr.write("output/ic_smooth.dat",minutes=True)
# save smoothed/baseline corrected TIC
# tic.write("output/tic.dat",minutes=True)
# tic1.write("output/tic_smooth.dat",minutes=True)
# tic2.write("output/tic_smooth_bc.dat",minutes=True)

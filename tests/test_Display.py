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

# stdlib
import sys
from typing import Iterator

# 3rd party
import pytest
from domdf_python_tools.paths import PathPlus
from matplotlib import axes, figure  # type: ignore[import]
from matplotlib import pyplot as plt

# this package
from pyms.Display import Display
from pyms.GCMS.Class import GCMS_data
from pyms.IntensityMatrix import IntensityMatrix
from pyms.IonChromatogram import IonChromatogram
from pyms.Spectrum import MassSpectrum

# this package
from .constants import *

if sys.version_info[:2] == (3, 6):
	baseline_dir = str(PathPlus(__file__).parent / "baseline")
	image_hashes = str(PathPlus(__file__).parent / "image_hashes_36.json")
elif sys.version_info[:2] == (3, 7):
	baseline_dir = str(PathPlus(__file__).parent / "baseline")
	image_hashes = str(PathPlus(__file__).parent / "image_hashes_37.json")
elif sys.version_info[:2] == (3, 8):
	baseline_dir = str(PathPlus(__file__).parent / "baseline")
	image_hashes = str(PathPlus(__file__).parent / "image_hashes_38.json")
else:
	baseline_dir = str(PathPlus(__file__).parent / "baseline")
	image_hashes = str(PathPlus(__file__).parent / "image_hashes.json")

check_images = pytest.mark.mpl_image_compare(
		baseline_dir=baseline_dir,
		savefig_kwargs={"dpi": 600},
		hash_library=image_hashes,
		)


def test_Display():
	no_args = Display()
	assert isinstance(no_args.fig, figure.Figure)
	assert isinstance(no_args.ax, axes.Axes)
	plt.close()

	fig = plt.figure()
	fig_arg = Display(fig=fig)
	assert isinstance(fig_arg.fig, figure.Figure)
	assert isinstance(fig_arg.ax, axes.Axes)
	assert fig_arg.fig is fig
	plt.close()

	fig = plt.figure()
	ax = fig.add_subplot(111)
	both_args = Display(fig=fig, ax=ax)
	assert isinstance(both_args.fig, figure.Figure)
	assert isinstance(both_args.ax, axes.Axes)
	assert both_args.fig is fig
	assert both_args.ax is ax
	plt.close()

	for obj in [test_tuple, test_list_strs, test_list_ints, test_string, *test_numbers, test_dict]:
		with pytest.raises(TypeError):
			Display(fig=obj)
		with pytest.raises(TypeError):
			Display(fig=fig, ax=obj)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	with pytest.raises(TypeError):
		Display(ax, fig)
	plt.close()


@pytest.fixture()
def test_plot() -> Iterator[Display]:
	fig = plt.figure()
	ax = fig.add_subplot(111)
	test_plot = Display(fig, ax)
	yield test_plot
	plt.close(fig)


@check_images
def test_plot_ic(im_i: IntensityMatrix, test_plot: Display):
	# Plotting IC with various Line2D options
	test_plot.plot_ic(im_i.get_ic_at_index(5))
	return test_plot.fig


@check_images
def test_plot_ic_label(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_ic(im_i.get_ic_at_index(5), label="IC @ Index 5")
	test_plot.ax.legend()
	return test_plot.fig


@check_images
def test_plot_ic_alpha(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_ic(im_i.get_ic_at_index(5), alpha=0.5)
	return test_plot.fig


@check_images
def test_plot_ic_linewidth(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_ic(im_i.get_ic_at_index(5), linewidth=2)
	return test_plot.fig


@check_images
def test_plot_ic_linestyle(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_ic(im_i.get_ic_at_index(5), linestyle="--")
	return test_plot.fig


@check_images
def test_plot_ic_multiple(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_ic(im_i.get_ic_at_index(5), label="IC @ Index 5")
	test_plot.plot_ic(im_i.get_ic_at_index(10), label="IC @ Index 10")
	test_plot.plot_ic(im_i.get_ic_at_index(20), label="IC @ Index 20")
	test_plot.plot_ic(im_i.get_ic_at_index(40), label="IC @ Index 40")
	test_plot.plot_ic(im_i.get_ic_at_index(80), label="IC @ Index 80")
	test_plot.plot_ic(im_i.get_ic_at_index(160), label="IC @ Index 160")
	test_plot.ax.legend()
	return test_plot.fig


@check_images
def test_plot_ic_title(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_ic(im_i.get_ic_at_index(5))
	test_plot.ax.set_title("Test IC Plot")
	return test_plot.fig


def test_plot_ic_errors(
		im_i: IntensityMatrix,
		test_plot: Display,
		data: GCMS_data,
		ms: MassSpectrum,
		):
	for obj in [*test_sequences, test_string, *test_numbers, test_dict, im_i, data, ms]:
		with pytest.raises(TypeError):
			test_plot.plot_ic(obj)  # type: ignore[arg-type]


# Plotting tic with various Line2D options
@check_images
def test_plot_tic(tic: IonChromatogram, test_plot: Display):
	test_plot.plot_tic(tic)
	return test_plot.fig


@check_images
def test_plot_tic_label(tic: IonChromatogram, test_plot: Display):
	test_plot.plot_tic(tic, label="IC @ Index 5")
	test_plot.ax.legend()
	return test_plot.fig


@check_images
def test_plot_tic_alpha(tic: IonChromatogram, test_plot: Display):
	test_plot.plot_tic(tic, alpha=0.5)
	return test_plot.fig


@check_images
def test_plot_tic_linewidth(tic: IonChromatogram, test_plot: Display):
	test_plot.plot_tic(tic, linewidth=2)
	return test_plot.fig


@check_images
def test_plot_tic_linestyle(tic: IonChromatogram, test_plot: Display):
	test_plot.plot_tic(tic, linestyle="--")
	return test_plot.fig


def test_plot_tic_errors(
		im_i: IntensityMatrix,
		test_plot: Display,
		data: GCMS_data,
		ms: MassSpectrum,
		):
	for obj in [
			*test_sequences,
			*test_numbers,
			test_string,
			test_dict,
			im_i,
			im_i.get_ic_at_index(0),
			data,
			ms,
			]:
		with pytest.raises(TypeError):
			test_plot.plot_tic(obj)  # type: ignore[arg-type]


@check_images
def test_plot_tic_title(tic: IonChromatogram, test_plot: Display):
	test_plot.plot_tic(tic)
	test_plot.ax.set_title("Test TIC Plot")
	return test_plot.fig


# Plotting mass spec with various Line2D options
@check_images
def test_plot_mass_spec(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50))
	return test_plot.fig


@check_images
def test_plot_mass_spec_alpha(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50), alpha=0.5)
	return test_plot.fig


@check_images
def test_plot_mass_spec_width(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50), width=1)
	return test_plot.fig


@check_images
def test_plot_mass_spec_linestyle(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50), linestyle="--")
	return test_plot.fig


def test_plot_mass_spec_errors(
		im_i: IntensityMatrix,
		test_plot: Display,
		data: GCMS_data,
		tic: IonChromatogram,
		):
	for obj in [
			*test_sequences,
			test_string,
			*test_numbers,
			test_dict,
			im_i,
			im_i.get_ic_at_index(0),
			data,
			tic,
			]:
		with pytest.raises(TypeError):
			test_plot.plot_mass_spec(obj)  # type: ignore[arg-type]


@check_images
def test_plot_mass_spec_title(im_i: IntensityMatrix, test_plot: Display):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50))
	test_plot.ax.set_title(f"Mass spec for peak at time {im_i.get_time_at_index(50):5.2f}")
	return test_plot.fig


def test_do_plotting_warning():
	test_plot = Display()

	with pytest.warns(UserWarning) as record:
		test_plot.do_plotting()

	# check that only one warning was raised
	assert len(record) == 1
	# check that the message matches
	args = record[0].message.args  # type: ignore[union-attr]
	expected = """No plots have been created.
Please call a plotting function before calling 'do_plotting()'"""
	assert args[0] == expected

	plt.close()

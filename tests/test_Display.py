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
import os
from pathlib import Path

# 3rd party
import pytest
from matplotlib import axes, figure  # type: ignore
from matplotlib import pyplot as plt  # type: ignore

# this package
from pyms.Display import Display

# this package
from .constants import *

baseline = str(Path(os.path.split(__file__)[0]) / "baseline")


def test_Display():
	no_args = Display()
	assert isinstance(no_args.fig, figure.Figure)
	assert isinstance(no_args.ax, axes.Axes)

	fig = plt.figure()
	fig_arg = Display(fig=fig)
	assert isinstance(fig_arg.fig, figure.Figure)
	assert isinstance(fig_arg.ax, axes.Axes)
	assert fig_arg.fig is fig

	fig = plt.figure()
	ax = fig.add_subplot(111)
	both_args = Display(fig=fig, ax=ax)
	assert isinstance(both_args.fig, figure.Figure)
	assert isinstance(both_args.ax, axes.Axes)
	assert both_args.fig is fig
	assert both_args.ax is ax

	for obj in [test_tuple, test_list_strs, test_list_ints, test_string, *test_numbers, test_dict]:
		with pytest.raises(TypeError):
			Display(fig=obj)
		with pytest.raises(TypeError):
			Display(fig=fig, ax=obj)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	with pytest.raises(TypeError):
		Display(ax, fig)


@pytest.fixture(scope="function")
def test_plot():
	fig = plt.figure()
	ax = fig.add_subplot(111)
	test_plot = Display(fig, ax)
	return test_plot


# Plotting IC with various Line2D options
@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_ic(im_i, test_plot):
	test_plot.plot_ic(im_i.get_ic_at_index(5))
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_ic_label(im_i, test_plot):
	test_plot.plot_ic(im_i.get_ic_at_index(5), label="IC @ Index 5")
	test_plot.ax.legend()
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_ic_alpha(im_i, test_plot):
	test_plot.plot_ic(im_i.get_ic_at_index(5), alpha=0.5)
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_ic_linewidth(im_i, test_plot):
	test_plot.plot_ic(im_i.get_ic_at_index(5), linewidth=2)
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_ic_linestyle(im_i, test_plot):
	test_plot.plot_ic(im_i.get_ic_at_index(5), linestyle="--")
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_ic_multiple(im_i, test_plot):
	test_plot.plot_ic(im_i.get_ic_at_index(5), label="IC @ Index 5")
	test_plot.plot_ic(im_i.get_ic_at_index(10), label="IC @ Index 10")
	test_plot.plot_ic(im_i.get_ic_at_index(20), label="IC @ Index 20")
	test_plot.plot_ic(im_i.get_ic_at_index(40), label="IC @ Index 40")
	test_plot.plot_ic(im_i.get_ic_at_index(80), label="IC @ Index 80")
	test_plot.plot_ic(im_i.get_ic_at_index(160), label="IC @ Index 160")
	test_plot.ax.legend()
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_ic_title(im_i, test_plot):
	test_plot.plot_ic(im_i.get_ic_at_index(5))
	test_plot.ax.set_title("Test IC Plot")
	return test_plot.fig


def test_plot_ic_errors(im_i, test_plot, data, ms):
	for obj in [*test_sequences, test_string, *test_numbers, test_dict, im_i, data, ms]:
		with pytest.raises(TypeError):
			test_plot.plot_ic(obj)


# Plotting tic with various Line2D options
@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_tic(tic, test_plot):
	test_plot.plot_tic(tic)
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_tic_label(tic, test_plot):
	test_plot.plot_tic(tic, label="IC @ Index 5")
	test_plot.ax.legend()
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_tic_alpha(tic, test_plot):
	test_plot.plot_tic(tic, alpha=0.5)
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_tic_linewidth(tic, test_plot):
	test_plot.plot_tic(tic, linewidth=2)
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_tic_linestyle(tic, test_plot):
	test_plot.plot_tic(tic, linestyle="--")
	return test_plot.fig


def test_plot_tic_errors(im_i, test_plot, data, ms):
	for obj in [*test_sequences, *test_numbers, test_string, test_dict, im_i, im_i.get_ic_at_index(0), data, ms]:
		with pytest.raises(TypeError):
			test_plot.plot_tic(obj)


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_tic_title(tic, test_plot):
	test_plot.plot_tic(tic)
	test_plot.ax.set_title("Test TIC Plot")
	return test_plot.fig


# Plotting mass spec with various Line2D options
@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_mass_spec(im_i, test_plot):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50))
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_mass_spec_alpha(im_i, test_plot):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50), alpha=0.5)
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_mass_spec_width(im_i, test_plot):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50), width=1)
	return test_plot.fig


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_mass_spec_linestyle(im_i, test_plot):
	test_plot.plot_mass_spec(im_i.get_ms_at_index(50), linestyle="--")
	return test_plot.fig


def test_plot_mass_spec_errors(im_i, test_plot, data, tic):
	for obj in [*test_sequences, test_string, *test_numbers, test_dict, im_i, im_i.get_ic_at_index(0), data, tic]:
		with pytest.raises(TypeError):
			test_plot.plot_mass_spec(obj)


@pytest.mark.mpl_image_compare(baseline_dir=baseline, savefig_kwargs={'dpi': 1200})
def test_plot_mass_spec_title(im_i, test_plot):
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
	args = record[0].message.args  # type: ignore
	assert args[0] == """No plots have been created.
Please call a plotting function before calling 'do_plotting()'"""

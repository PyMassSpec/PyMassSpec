"""
Class to Display Ion Chromatograms and TIC.
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
import warnings
from typing import Dict, List, Optional, Tuple

# 3rd party
import deprecation  # type: ignore[import]
import matplotlib  # type: ignore[import]
import matplotlib.pyplot as plt  # type: ignore[import]
from matplotlib.axes import Axes  # type: ignore[import]
from matplotlib.container import BarContainer  # type: ignore[import]
from matplotlib.figure import Figure  # type: ignore[import]
from matplotlib.lines import Line2D  # type: ignore[import]

# this package
from pyms import Peak, __version__
from pyms.IonChromatogram import IonChromatogram
from pyms.Peak.List.Function import is_peak_list
from pyms.Spectrum import MassSpectrum, normalize_mass_spec

__all__ = [
		"Display",
		"plot_ic",
		"plot_mass_spec",
		"plot_head2tail",
		"plot_peaks",
		"ClickEventHandler",
		"invert_mass_spec",
		]

default_filetypes = ["png", "pdf", "svg"]

# Ensure that the intersphinx links are correct.
Axes.__module__ = "matplotlib.axes"
Figure.__module__ = "matplotlib.figure"


class Display:
	"""
	Class to display Ion Chromatograms and Total Ion Chromatograms from
	:class:`pyms.IonChromatogram.IonChromatogram` using :mod:`matplotlib.pyplot`.

	:param fig: figure object to use
	:param ax: axes object to use

	If ``fig`` is not given then ``fig`` and ``ax`` default to:

	>>> fig = plt.figure()
	>>> ax = fig.add_subplot(111)


	If only ``fig`` is given then ``ax`` defaults to:

	>>> ax = fig.add_subplot(111)

	:author: Sean O'Callaghan
	:author: Vladimir Likic
	:author: Dominic Davis-Foster
	"""  # noqa: D400

	@deprecation.deprecated(
			"2.2.8",
			"2.4.0",
			current_version=__version__,
			details="Functionality has moved to other functions and classes in this module.",
			)
	def __init__(self, fig: Figure = None, ax: Axes = None):

		if fig is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)

		elif isinstance(fig, matplotlib.figure.Figure) and ax is None:
			ax = fig.add_subplot(111)

		if not isinstance(fig, matplotlib.figure.Figure):
			raise TypeError("'fig' must be a matplotlib.figure.Figure object")

		if not isinstance(ax, matplotlib.axes.Axes):
			raise TypeError("'ax' must be a matplotlib.axes.Axes object")

		self.fig = fig
		self.ax = ax

		# Container to store plots
		self.__tic_ic_plots: List[List[Line2D]] = []

		# Peak list container
		self.__peak_list: List[Peak.Peak] = []

	def do_plotting(self, plot_label: Optional[str] = None) -> None:
		"""
		Plots TIC and IC(s) if they have been created by
		:meth:`~pyms.Display.Display.plot_tic` or
		:meth:`~pyms.Display.Display.plot_ic`.

		Also adds detected peaks if they have been added by
		:meth:`~pyms.Display.Display.plot_peaks`

		:param plot_label: Label for the plot to show e.g. the data origin
		"""  # noqa: D400

		# if no plots have been created advise user
		if len(self.__tic_ic_plots) == 0:
			warnings.warn(
					"""No plots have been created.
Please call a plotting function before calling 'do_plotting()'""",
					UserWarning,
					)
			return

		if plot_label is not None:
			self.ax.set_title(plot_label)

		self.ax.legend()

		self.fig.canvas.draw()

		# If no peak list plot, no mouse click event
		if len(self.__peak_list) != 0:
			self.fig.canvas.mpl_connect("button_press_event", self.onclick)

	# plt.show()

	@staticmethod
	def get_5_largest(intensity_list: List[float]) -> List[int]:
		"""
		Returns the indices of the 5 largest ion intensities.

		:param intensity_list: List of Ion intensities
		"""

		largest = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		# Find out largest value
		for idx, intensity in enumerate(intensity_list):
			if intensity > intensity_list[largest[0]]:
				largest[0] = idx

		# Now find next four largest values
		for j in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
			for idx, intensity in enumerate(intensity_list):
				# if intensity_list[i] > intensity_list[largest[j]] and intensity_list[i] < intensity_list[largest[j-1]]:
				if intensity_list[largest[j]] < intensity < intensity_list[largest[j - 1]]:
					largest[j] = idx

		return largest

	def onclick(self, event) -> None:
		"""
		Finds the 5 highest intensity m/z channels for the selected peak.
		The peak is selected by clicking on it.
		If a button other than the left one is clicked, a new plot of the mass spectrum is displayed.

		:param event: a mouse click by the user
		"""

		intensity_list = []
		mass_list = []

		for peak in self.__peak_list:
			# if event.xdata > 0.9999*peak.rt and event.xdata < 1.0001*peak.rt:
			if 0.9999 * peak.rt < event.xdata < 1.0001 * peak.rt:
				intensity_list = peak.mass_spectrum.mass_spec
				mass_list = peak.mass_spectrum.mass_list

		largest = self.get_5_largest(intensity_list)

		if len(intensity_list) != 0:
			print("mass\t intensity")
			for i in range(10):
				print(mass_list[largest[i]], '\t', intensity_list[largest[i]])
		else:  # if the selected point is not close enough to peak
			print("No Peak at this point")

		# Check if a button other than left was pressed, if so plot mass spectrum
		# Also check that a peak was selected, not just whitespace
		if event.button != 1 and len(intensity_list) != 0:
			# self.plot_mass_spec(event.xdata, mass_list, intensity_list)
			self.plot_mass_spec(MassSpectrum(mass_list, intensity_list))

	def plot_ic(self, ic: IonChromatogram, **kwargs) -> List[Line2D]:
		"""
		Plots an Ion Chromatogram.

		:param ic: Ion Chromatograms m/z channels for plotting

		:Other Parameters: :class:`matplotlib.lines.Line2D` properties.
			Used to specify properties like a line label (for auto legends),
			linewidth, antialiasing, marker face color.

			Example::

			>>> plot_ic(im.get_ic_at_index(5), label='IC @ Index 5', linewidth=2)

			See https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.lines.Line2D.html
			for the list of possible kwargs
		"""

		plot = plot_ic(self.ax, ic, **kwargs)
		self.__tic_ic_plots.append(plot)
		return plot

	def plot_mass_spec(self, mass_spec: MassSpectrum, **kwargs) -> BarContainer:
		"""
		Plots a Mass Spectrum.

		:param mass_spec: The mass spectrum at a given time/index

		:Other Parameters: :class:`matplotlib.lines.Line2D` properties.
			Used to specify properties like a line label (for auto legends),
			linewidth, antialiasing, marker face color.

			Example::

			>>> plot_mass_spec(im.get_ms_at_index(5), linewidth=2)
			>>>	ax.set_title(f"Mass spec for peak at time {im.get_time_at_index(5):5.2f}")

			See https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.lines.Line2D.html
			for the list of possible kwargs
		"""

		plot = plot_mass_spec(self.ax, mass_spec, **kwargs)
		return plot

	def plot_peaks(self, peak_list: List[Peak.Peak], label: str = "Peaks") -> List[Line2D]:
		"""
		Plots the locations of peaks as found by PyMassSpec.

		:param peak_list: List of peaks to plot
		:param label: label for plot legend.
		"""

		plot = plot_peaks(self.ax, peak_list, label)

		# Copy to self.__peak_list for onclick event handling
		self.__peak_list = peak_list

		return plot

	def plot_tic(self, tic: IonChromatogram, minutes: bool = False, **kwargs) -> List[Line2D]:
		"""
		Plots a Total Ion Chromatogram.

		:param tic: Total Ion Chromatogram.
		:param minutes: Whether to show the time in minutes.

		:Other Parameters: :class:`matplotlib.lines.Line2D` properties.
			Used to specify properties like a line label (for auto legends),
			linewidth, antialiasing, marker face color.

			Example::

			>>> plot_tic(data.tic, label='TIC', linewidth=2)

			See https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.lines.Line2D.html
			for the list of possible kwargs
		"""

		if not isinstance(tic, IonChromatogram) or not tic.is_tic():
			raise TypeError("'tic' must be an Ion Chromatogram object representing a total ion chromatogram")

		plot = plot_ic(self.ax, tic, minutes, **kwargs)
		self.__tic_ic_plots.append(plot)
		return plot

	def save_chart(self, filepath: str, filetypes: Optional[List[str]] = None) -> None:
		"""
		Save the chart to the given path with the given filetypes.

		:param filepath: Path and filename to save the chart as. Should not include extension.
		:param filetypes: List of filetypes to use.

		:author: Dominic Davis-Foster
		"""

		# TODO: pathlib and remove extension if given & use that as filetype

		if filetypes is None:
			filetypes = default_filetypes

		# matplotlib.use("Agg")

		for filetype in filetypes:
			# plt.savefig(filepath + ".{}".format(filetype))
			self.fig.savefig(filepath + f".{filetype}")
		plt.close()

	def show_chart(self) -> None:
		"""
		Show the chart on screen.

		:author: Dominic Davis-Foster
		"""

		# matplotlib.use("TkAgg")

		self.fig.show()
		input("Press Enter to close the chart")
		plt.close()


def plot_ic(ax: matplotlib.axes.Axes, ic: IonChromatogram, minutes: bool = False, **kwargs) -> List[Line2D]:
	"""
	Plots an Ion Chromatogram.

	:param ax: The axes to plot the IonChromatogram on
	:param ic: Ion Chromatograms m/z channels for plotting
	:param minutes: Whether the x-axis should be plotted in minutes. Default False (plotted in seconds)

	:Other Parameters: :class:`matplotlib.lines.Line2D` properties.
		Used to specify properties like a line label (for auto legends),
		linewidth, antialiasing, marker face color.

		Example::

		>>> plot_ic(im.get_ic_at_index(5), label='IC @ Index 5', linewidth=2)

		See https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.lines.Line2D.html
		for the list of possible kwargs

	:return: A list of Line2D objects representing the plotted data.
	"""

	if not isinstance(ic, IonChromatogram):
		raise TypeError("'ic' must be an IonChromatogram")

	time_list = ic.time_list
	if minutes:
		time_list = [time / 60 for time in time_list]

	plot = ax.plot(time_list, ic.intensity_array, **kwargs)

	# Set axis ranges
	ax.set_xlim(min(time_list), max(time_list))
	ax.set_ylim(bottom=0)

	return plot


def plot_mass_spec(ax: Axes, mass_spec: MassSpectrum, **kwargs) -> BarContainer:
	"""
	Plots a Mass Spectrum.

	:param ax: The axes to plot the MassSpectrum on
	:param mass_spec: The mass spectrum to plot

	:Other Parameters: :class:`matplotlib.lines.Line2D` properties.
		Used to specify properties like a line label (for auto legends),
		linewidth, antialiasing, marker face color.

		Example::

		>>> plot_mass_spec(im.get_ms_at_index(5), linewidth=2)
		>>>	ax.set_title(f"Mass spec for peak at time {im.get_time_at_index(5):5.2f}")

		See https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.lines.Line2D.html
		for the list of possible kwargs

	:return: Container with all the bars and optionally errorbars.
	:rtype: :class:`matplotlib.container.BarContainer`
	"""

	if not isinstance(mass_spec, MassSpectrum):
		raise TypeError("'mass_spec' must be a MassSpectrum")

	mass_list = mass_spec.mass_list
	intensity_list = mass_spec.mass_spec

	if "width" not in kwargs:
		kwargs["width"] = 0.5

	# to set x axis range find minimum and maximum m/z channels
	min_mz = mass_list[0]
	max_mz = mass_list[-1]

	for idx, mass in enumerate(mass_list):
		if mass_list[idx] > max_mz:
			max_mz = mass_list[idx]

	for idx, mass in enumerate(mass_list):
		if mass_list[idx] < min_mz:
			min_mz = mass_list[idx]

	plot = ax.bar(mass_list, intensity_list, **kwargs)

	# Set axis ranges
	ax.set_xlim(min_mz - 1, max_mz + 1)
	ax.set_ylim(bottom=0)

	return plot


def plot_head2tail(
		ax: Axes,
		top_mass_spec: MassSpectrum,
		bottom_mass_spec: MassSpectrum,
		top_spec_kwargs: Optional[Dict] = None,
		bottom_spec_kwargs: Optional[Dict] = None,
		) -> Tuple[BarContainer, BarContainer]:
	"""
	Plots two mass spectra head to tail.

	:param ax: The axes to plot the MassSpectra on

	:param top_mass_spec: The Mass Spectrum to plot on top
	:param bottom_mass_spec: The Mass Spectrum to plot on the bottom
	:param top_spec_kwargs: A dictionary of keyword arguments for the top mass spectrum.
		Defaults to red with a line width of 0.5
	:no-default top_spec_kwargs:
	:param bottom_spec_kwargs: A dictionary of keyword arguments for the bottom mass spectrum.
		Defaults to blue with a line width of 0.5
	:no-default bottom_spec_kwargs:

	`top_spec_kwargs` and `bottom_spec_kwargs` are used to specify properties like a line label
		(for auto legends), linewidth, antialiasing, marker face color.

		See https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.lines.Line2D.html
		for the list of possible kwargs

	:return: A tuple of container with all the bars and optionally errorbars for the top and bottom spectra.
	:rtype: tuple of :class:`matplotlib.container.BarContainer`
	"""

	if not isinstance(top_mass_spec, MassSpectrum):
		raise TypeError("'top_mass_spec' must be a MassSpectrum")

	if not isinstance(bottom_mass_spec, MassSpectrum):
		raise TypeError("'bottom_mass_spec' must be a MassSpectrum")

	if top_spec_kwargs is None:
		top_spec_kwargs = dict(color="red", width=0.5)
	elif not isinstance(top_spec_kwargs, dict):
		raise TypeError("'top_spec_kwargs' must be a dictionary of keyword arguments for the top mass spectrum.")

	if bottom_spec_kwargs is None:
		bottom_spec_kwargs = dict(color="blue", width=0.5)
	elif not isinstance(bottom_spec_kwargs, dict):
		raise TypeError(
				"'bottom_spec_kwargs' must be a dictionary of keyword arguments for the bottom mass spectrum."
				)

	# Plot a line at y=0 with same width and colour as Spines
	ax.axhline(y=0, color=ax.spines["bottom"].get_edgecolor(), linewidth=ax.spines["bottom"].get_linewidth())

	# Normalize the mass spectra
	top_mass_spec = normalize_mass_spec(top_mass_spec)
	bottom_mass_spec = normalize_mass_spec(bottom_mass_spec)

	# Invert bottom mass spec
	invert_mass_spec(bottom_mass_spec, inplace=True)

	top_plot = plot_mass_spec(ax, top_mass_spec, **top_spec_kwargs)
	bottom_plot = plot_mass_spec(ax, bottom_mass_spec, **bottom_spec_kwargs)

	# Set ylim to 1.1 times max/min values
	ax.set_ylim(
			bottom=min(bottom_mass_spec.intensity_list) * 1.1,
			top=max(top_mass_spec.intensity_list) * 1.1,
			)

	# ax.spines['bottom'].set_position('zero')

	return top_plot, bottom_plot


def plot_peaks(ax: Axes, peak_list: List[Peak.Peak], label: str = "Peaks", style: str = 'o') -> List[Line2D]:
	"""
	Plots the locations of peaks as found by PyMassSpec.

	:param ax: The axes to plot the peaks on
	:param peak_list: List of peaks to plot
	:param label: label for plot legend.
	:param style: The marker style. See `https://matplotlib.org/3.1.1/api/markers_api.html` for a complete list

	:return: A list of Line2D objects representing the plotted data.
	"""

	if not is_peak_list(peak_list):
		raise TypeError("'peak_list' must be a list of Peak objects")

	time_list = []
	height_list = []

	if "line" in style.lower():
		lines = []
		for peak in peak_list:
			lines.append(ax.axvline(x=peak.rt, color="lightgrey", alpha=0.8, linewidth=0.3))

		return lines

	else:
		for peak in peak_list:
			time_list.append(peak.rt)
			height_list.append(sum(peak.mass_spectrum.intensity_list))
			# height_list.append(peak.height)
			# print(peak.height - sum(peak.mass_spectrum.intensity_list))
			# print(sum(peak.mass_spectrum.intensity_list))

		return ax.plot(time_list, height_list, style, label=label)


# TODO: Change order of arguments and use plt.gca() a la pyplot


class ClickEventHandler:
	"""
	Class to enable clicking of chromatogram to view the intensities top n most intense
	ions at that peak, and viewing of the mass spectrum with a right click
	"""  # noqa: D400

	peak_list: List[Peak.Peak]
	fig: Figure
	ax: Axes
	ms_fig: Optional[Figure]
	ms_ax: Optional[Axes]
	n_intensities: int

	def __init__(
			self,
			peak_list: List[Peak.Peak],
			fig: Optional[Figure] = None,
			ax: Optional[Axes] = None,
			tolerance: float = 0.005,
			n_intensities: int = 5,
			):
		if fig is None:
			self.fig = plt.gcf()
		else:
			self.fig = fig

		if ax is None:
			self.ax = plt.gca()
		else:
			self.ax = ax

		self.peak_list = peak_list

		self.ms_fig: Optional[Figure] = None
		self.ms_ax: Optional[Axes] = None

		self._min = 1 - tolerance
		self._max = 1 + tolerance
		self.n_intensities = n_intensities

		# If no peak list plot, no mouse click event
		if len(self.peak_list) != 0:
			self.cid = self.fig.canvas.mpl_connect("button_press_event", self.onclick)
		else:
			self.cid = None

	def onclick(self, event) -> None:
		"""
		Finds the n highest intensity m/z channels for the selected peak.
		The peak is selected by clicking on it.
		If a button other than the left one is clicked, a new plot of the mass spectrum is displayed.

		:param event: a mouse click by the user
		"""

		for peak in self.peak_list:
			# if event.xdata > 0.9999*peak.rt and event.xdata < 1.0001*peak.rt:
			if self._min * peak.rt < event.xdata < self._max * peak.rt:
				intensity_list = peak.mass_spectrum.mass_spec
				mass_list = peak.mass_spectrum.mass_list

				largest = self.get_n_largest(intensity_list)

				print(f"RT: {peak.rt}")
				print("Mass\t Intensity")
				for i in range(self.n_intensities):
					print(f"{mass_list[largest[i]]}\t {intensity_list[largest[i]]}")

				# Check if right mouse button pressed, if so plot mass spectrum
				# Also check that a peak was selected, not just whitespace
				if event.button == 3 and len(intensity_list) != 0:
					# from pyms.Display import plot_mass_spec

					if self.ms_fig is None or self.ms_ax is None:
						self.ms_fig, self.ms_ax = plt.subplots(1, 1)
					else:
						self.ms_ax.clear()

					plot_mass_spec(self.ms_ax, peak.mass_spectrum)
					self.ms_ax.set_title(f"Mass Spectrum at RT {peak.rt}")
					self.ms_fig.show()
				# TODO: Add multiple MS to same plot window and add option to close one of them
				# TODO: Allow more interaction with MS, e.g. adjusting mass range?
				return

		# if the selected point is not close enough to peak
		print("No Peak at this point")

	def get_n_largest(self, intensity_list: List[float]) -> List[int]:
		"""
		Computes the indices of the largest n ion intensities for writing to console.

		:param intensity_list: List of Ion intensities

		:return: Indices of largest ``n`` ion intensities
		"""

		largest = [0] * self.n_intensities

		# Find out largest value
		for idx, intensity in enumerate(intensity_list):
			if intensity > intensity_list[largest[0]]:
				largest[0] = idx

		# Now find next four largest values
		for j in list(range(1, self.n_intensities)):
			for idx, intensity in enumerate(intensity_list):
				# if intensity_list[i] > intensity_list[largest[j]] and intensity_list[i] < intensity_list[largest[j-1]]:
				if intensity_list[largest[j]] < intensity < intensity_list[largest[j - 1]]:
					largest[j] = idx

		return largest


def invert_mass_spec(mass_spec: MassSpectrum, inplace: bool = False) -> MassSpectrum:
	"""
	Invert the mass spectrum for display in a head2tail plot.

	:param mass_spec: The Mass Spectrum to normalize
	:param inplace: Whether the inversion should be applied to the
		:class:`~pyms.Spectrum.MassSpectrum` object given, or to a copy (default behaviour).

	:return: The normalized mass spectrum
	"""

	inverted_intensity_list = [-x for x in mass_spec.intensity_list]

	if inplace:
		mass_spec.intensity_list = inverted_intensity_list
		return mass_spec
	else:
		return MassSpectrum(mass_spec.mass_list, inverted_intensity_list)

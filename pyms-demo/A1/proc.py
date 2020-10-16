"""
proc.py

Plot detected peaks using matplotlib
"""

import sys
sys.path.append("../..")

import pathlib

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from pyms.BillerBiemann import BillerBiemann, rel_threshold, num_ions_threshold
from pyms.Display import plot_ic, plot_peaks
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat


data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

# Read raw data
andi_file = data_directory / "MM-10.0_1_no_processing.cdf"
data = ANDI_reader(andi_file)

# Build Intensity Matrix
im = build_intensity_matrix_i(data)

# Perform pre-filtering and peak detection.

n_scan, n_mz = im.size

for ii in range(n_mz):
	ic = im.get_ic_at_index(ii)
	ic_smooth = savitzky_golay(ic)
	ic_bc = tophat(ic_smooth, struct="1.5m")
	im.set_ic_at_index(ii, ic_bc)

# Detect Peaks
peak_list = BillerBiemann(im, points=9, scans=2)

print("Number of peaks found: ", len(peak_list))

# Filter the peak list, first by removing all intensities in a peak less than a
# given relative threshold, then by removing all peaks that have less than a
# given number of ions above a given value

# Parameters
# percentage ratio of ion intensity to max ion intensity
percent = 2
# minimum number of ions, n
n = 3
# greater than or equal to threshold, t
cutoff = 10000

# trim by relative intensity
pl = rel_threshold(peak_list, percent)

# trim by threshold
new_peak_list = num_ions_threshold(pl, n, cutoff)

print("Number of filtered peaks: ", len(new_peak_list))

# TIC from raw data
tic = data.tic

# Get Ion Chromatograms for all m/z channels
n_mz = len(im.mass_list)

# Create a subplot
fig, ax = plt.subplots(1, 1)

# Plot the peaks
plot_peaks(ax, new_peak_list, style="lines")
# Note: No idea why, but the dots for the peaks consistently appear 2e7 below the apex of the peak.
# As an alternative, the positions of the peaks can be shown with thin grey lines as in this example.
# The peak positions seem to appear OK in the other examples.
# See pyms-demo/scripts/Displaying_Detected_Peaks.py for a better example

# Plot the TIC
plot_ic(ax, tic, label="TIC")

# Plot the ICs
# for m in range(n_mz):
# 	plot_ic(ax, im.get_ic_at_index(m))

# Set the title
ax.set_title('TIC and PyMS Detected Peaks')

# Add the legend
plt.legend()

# Show the plot
plt.show()

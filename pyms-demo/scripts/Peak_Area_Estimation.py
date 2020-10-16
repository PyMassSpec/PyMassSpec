#!/usr/bin/env python
# coding: utf-8

# # Peak Area Estimation
#
# The |Peak| object does not contain any information about the width or area of
# the peak when it is first created. This information can be added after the
# instantiation of a Peak object. The area of the peak can be set with the
# attribute |pyms.Peak.Class.Peak.area|.
#
# The total peak area can by obtained by the |peak_sum_area()| function in
# |pyms.Peak.Function|. The function determines the total area as the sum of the
# ion intensities for all masses that apex at the given peak. To calculate the
# peak area of a single mass, the intensities are added from the apex of the
# mass peak outwards.
#
# Edge values are added until the following conditions are met:
# * the added intensity adds less than 0.5\% to the accumulated area; or
# * the added intensity starts increasing (i.e. when the ion is common to
#   co-eluting compounds).
#
# To avoid noise effects, the edge value is taken at the midpoint of three consecutive edge values.
#
# First, build the Peak list as before

# In[ ]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat
from pyms.BillerBiemann import BillerBiemann

jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
im = build_intensity_matrix(data)

n_scan, n_mz = im.size

for ii in range(n_mz):
    ic = im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    im.set_ic_at_index(ii, ic_bc)

peak_list = BillerBiemann(im, points=9, scans=2)

from pyms.Noise.Analysis import window_analyzer
tic = data.tic
noise_level = window_analyzer(tic)

from pyms.BillerBiemann import num_ions_threshold
filtered_peak_list = num_ions_threshold(peak_list, n=3, cutoff=noise_level)
print(filtered_peak_list[:10])


# Given a list of peaks, areas can be determined and added as follows:

# In[ ]:


from pyms.Peak.Function import peak_sum_area
for peak in peak_list:
    area = peak_sum_area(im, peak)
    peak.area = area

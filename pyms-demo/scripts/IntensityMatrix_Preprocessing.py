#!/usr/bin/env python
# coding: utf-8

# # Pre-processing the IntensityMatrix
#
# Noise smoothing and baseline correction can be applied to each |IonChromatogram| in an |IntensityMatrix|.
#
# First, setup the paths to the datafiles and the output directory, then import the required functions.

# In[1]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat


# Read the raw data files and build the |IntensityMatrix|:

# In[2]:


jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
im = build_intensity_matrix(data)


# Perform Savitzky-Golay smoothing and Tophat baseline correction
#

# In[4]:


n_scan, n_mz = im.size

for ii in range(n_mz):
    # print("Working on IC#", ii+1)
    ic = im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    im.set_ic_at_index(ii, ic_bc)


# Alternatively, the filtering may be performed on the |IntensityMatrix| without
# using a ``for`` loop, as outlined in previous examples. However filtering by
# |IonChromatogram| in a ``for`` loop as described here is much faster.
#

#!/usr/bin/env python

# ## Example: Peak Detection
#
# First, setup the paths to the datafiles and the output directory, then import JCAMP_reader and build_intensity_matrix.

# In[1]:

import pathlib

data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix

# Read the raw data file and build the |IntensityMatrix|.

# In[2]:

jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
im = build_intensity_matrix(data)

# Preprocess the data (Savitzky-Golay smoothing and Tophat baseline detection).

# In[3]:

from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat

n_scan, n_mz = im.size

for ii in range(n_mz):
	ic = im.get_ic_at_index(ii)
	ic_smooth = savitzky_golay(ic)
	ic_bc = tophat(ic_smooth, struct="1.5m")
	im.set_ic_at_index(ii, ic_bc)

# Now the Biller and Biemann based technique can be applied to detect peaks.

# In[4]:

from pyms.BillerBiemann import BillerBiemann

peak_list = BillerBiemann(im)
print(peak_list[:10])

# In[5]:

len(peak_list)

# Note that this is nearly as many peaks as there are scans in the data
# (9865 scans). This is due to noise and the simplicity of the technique.
#
# The number of detected peaks can be constrained by the selection of better
# parameters. Parameters can be determined by counting the number of points
# across a peak, and examining where peaks are found. For example, the peak
# list can be found with the parameters of a window of 9 points and by
# combining 2 neighbouring scans if they apex next to each other:
#

# In[6]:

peak_list = BillerBiemann(im, points=9, scans=2)
print(peak_list[:10])

# In[7]:

len(peak_list)

# The number of detected peaks has been reduced, but there are still many more
# than would be expected from the sample. Functions to filter the peak list are
# covered in the next example.
#
# ## Example: Peak List Filtering
#
# There are two functions to filter the list of Peak objects.
#
# The first, |rel_threshold()| modifies the mass spectrum stored in each peak so
# any intensity that is less than a given percentage of the maximum intensity for the peak is removed.
#
# The second, |num_ions_threshold()|, removes any peak that has less than a given
# number of ions above a given threshold.
#
# Once the peak list has been constructed, the filters can be applied as follows:

# In[8]:

from pyms.BillerBiemann import num_ions_threshold, rel_threshold

pl = rel_threshold(peak_list, percent=2)
print(pl[:10])

# In[9]:

new_peak_list = num_ions_threshold(pl, n=3, cutoff=10000)
print(new_peak_list[:10])

# In[10]:

len(new_peak_list)

# The number of detected peaks is now more realistic of what would be expected in
# the test sample.

#!/usr/bin/env python

# ## Example: Creating an Experiment
#
# Before aligning peaks from multiple experiments, the peak objects need to be
# created and encapsulated into |Experiment| objects. During this process it is
# often useful to pre-process the peaks in some way, for example to null certain
# m/z channels and/or to select a certain retention time range.
#
# The procedure starts the same as in the previous examples, namely:
#
# 1. read a file,
# 2. bin the data into fixed mass values,
# 3. smooth the data,
# 4. remove the baseline,
# 5. deconvolute peaks,
# 6. filter the peaks,
# 7. set the mass range,
# 8. remove uninformative ions, and
# 9. estimate peak areas.
#
#
# First, setup the paths to the datafiles and the output directory, then import ANDI_reader and build_intensity_matrix_i.

# In[1]:

import pathlib

data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i

# Read the raw data file and build the |IntensityMatrix|.

# In[2]:

andi_file = data_directory / "a0806_077.cdf"
data = ANDI_reader(andi_file)
im = build_intensity_matrix_i(data)

# Preprocess the data (Savitzky-Golay smoothing and Tophat baseline detection)

# In[3]:

from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat

n_scan, n_mz = im.size

for ii in range(n_mz):
	ic = im.get_ic_at_index(ii)
	ic1 = savitzky_golay(ic)
	ic_smooth = savitzky_golay(ic1)  # Why the second pass here?
	ic_bc = tophat(ic_smooth, struct="1.5m")
	im.set_ic_at_index(ii, ic_bc)

# Now the Biller and Biemann based technique can be applied to detect peaks.

# In[4]:

from pyms.BillerBiemann import BillerBiemann

pl = BillerBiemann(im, points=9, scans=2)
len(pl)

# Trim the peak list by relative intensity

# In[5]:

from pyms.BillerBiemann import num_ions_threshold, rel_threshold

apl = rel_threshold(pl, percent=2)
len(apl)

# Trim the peak list by noise threshold

# In[6]:

peak_list = num_ions_threshold(apl, n=3, cutoff=3000)
len(peak_list)

# Set the mass range, remove unwanted ions and estimate the peak area

# In[7]:

from pyms.Peak.Function import peak_sum_area

for peak in peak_list:
	peak.crop_mass(51, 540)

	peak.null_mass(73)
	peak.null_mass(147)

	area = peak_sum_area(im, peak)
	peak.area = area

# Create an |Experiment|.

# In[8]:

from pyms.Experiment import Experiment

expr = Experiment("a0806_077", peak_list)

# Set the time range for all Experiments

# In[9]:

expr.sele_rt_range(["6.5m", "21m"])

# Save the experiment to disk.

# In[10]:

expr.dump(output_directory / "experiments" / "a0806_077.expr")

#

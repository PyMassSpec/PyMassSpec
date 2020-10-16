#!/usr/bin/env python
# coding: utf-8

# # Noise smoothing
#
# The purpose of noise smoothing is to remove high-frequency noise from
# data, and thereby increase the contribution of the signal relative to
# the contribution of the noise.
#
# First, setup the paths to the datafiles and the output directory, then import JCAMP_reader.

# In[1]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader


# Read the raw data files and extract the TIC.

# In[2]:


jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
tic = data.tic


# ## Window averaging
#
# A simple approach to noise smoothing is moving average window smoothing.
# In this approach the window of a fixed size (:math:`2N+1` points) is moved
# across the ion chromatogram, and the intensity value at each point is
# replaced with the mean intensity calculated over the window size.
# The example below illustrates smoothing of TIC by window averaging.
#
# To apply mean window smoothing with a 5-point window:

# In[3]:


from pyms.Noise.Window import window_smooth
tic1 = window_smooth(tic, window=5)


# To apply median window smoothing with a 5-point window:

# In[4]:


tic2 = window_smooth(tic, window=5, use_median=True)


# To apply the mean windows smoothing, but specifying the window as
# a time string (in this example, 7 seconds):

# In[5]:


tic3 = window_smooth(tic, window='7s')


# Write the original TIC and the smoothed TICs to disk:

# In[6]:


tic.write(output_directory / "noise_smoothing_tic.dat",minutes=True)
tic1.write(output_directory / "noise_smoothing_tic1.dat",minutes=True)
tic2.write(output_directory / "noise_smoothing_tic2.dat",minutes=True)


# ## Window Averaging on Intensity Matrix
#
# In the previous section, window averaging was applied to an
# Ion Chromatogram object (in that case a TIC). Where filtering
# is to be performed on all Ion Chromatograms, the
# |window_smooth_im()| function may be used instead.
#
# The use of this function is identical to the Ion Chromatogram
# |window_smooth()| function, except that an Intensity Matrix
# is passed to it.
#
# For example, to perform window smoothing on an |IntensityMatrix|
# object with a 5 point window and mean window smoothing:

# In[7]:


from pyms.IntensityMatrix import build_intensity_matrix
from pyms.Noise.Window import window_smooth_im
im = build_intensity_matrix(data)
im_smooth1 = window_smooth_im(im, window = 5, use_median = False)


# Write the IC for mass 73 to disk for both the original and smoothed |IntensityMatrix|:

# In[8]:


ic = im.get_ic_at_index(73)
ic_smooth1 = im_smooth1.get_ic_at_index(73)

ic.write(output_directory/"noise_smoothing_ic.dat",minutes=True)
ic_smooth1.write(output_directory/"noise_smoothing_ic_smooth1.dat",minutes=True)


# ## Savitzky--Golay noise filter
#
# A more sophisticated noise filter is the Savitzky-Golay filter.
# Given the data loaded as above, this filter can be applied as
# follows:

# In[9]:


from pyms.Noise.SavitzkyGolay import savitzky_golay
tic4 = savitzky_golay(tic)


# Write the smoothed TIC to disk:

# In[10]:


tic4.write(output_directory / "noise_smoothing_tic4.dat",minutes=True)


# In this example the default parameters were used.
#
# ### Savitzky-Golay Noise filtering of Intensity Matrix Object
#
# The |savitzky_golay()| function described above acts on a single
# |IonChromatogram|. Where it is desired to perform Savitzky Golay
# filtering on the whole |IntensityMatrix| the function
# |savitzky_golay_im()| may be used as follows:

# In[11]:


from pyms.Noise.SavitzkyGolay import savitzky_golay_im
im_smooth2 = savitzky_golay_im(im)


# Write the IC for mass 73 in the smoothed |IntensityMatrix| to disk:

# In[12]:


ic_smooth2 = im_smooth2.get_ic_at_index(73)
ic_smooth2.write(output_directory/"noise_smoothing_ic_smooth2.dat",minutes=True)


#

#!/usr/bin/env python
# coding: utf-8

# # Baseline Correction
#
# Baseline distortion originating from instrument imperfections and experimental
# setup is often observed in mass spectrometry data, and off-line baseline
# correction is often an important step in data pre-processing. There are many
# approaches for baseline correction. One advanced approach is based on the
# top-hat transform developed in mathematical morphology, and used extensively in
# digital image processing for tasks such as image enhancement. Top-hat baseline
# correction was previously applied in proteomics based mass spectrometry.
# PyMS currently implements only the top-hat baseline corrector, using the SciPy
# package ``ndimage``.
#
# Application of the top-hat baseline corrector requires the size of the
# structural element to be specified. The structural element needs to be larger
# than the features one wants to retain in the spectrum after the top-hat
# transform. In the example below, the top-hat baseline corrector is applied to
# the TIC of the data set ``gc01_0812_066.cdf``, with the structural element of
# 1.5 minutes:
#
# The purpose of noise smoothing is to remove high-frequency noise from data,
# and thereby increase the contribution of the signal relative to the contribution
# of the noise.
#
# First, setup the paths to the datafiles and the output directory, then import ANDI_reader, savitzky_golay and tophat.

# In[1]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat


# Read the raw data files and extract the TIC.

# In[2]:


andi_file = data_directory / "gc01_0812_066.cdf"
data = ANDI_reader(andi_file)
tic = data.tic


# Perform Savitzky-Golay smoothing
#

# In[3]:


tic1 = savitzky_golay(tic)


# Perform Tophat baseline correction

# In[4]:


tic2 = tophat(tic1, struct="1.5m")


# Save the output to disk

# In[5]:


tic.write(output_directory / "baseline_correction_tic.dat",minutes=True)
tic1.write(output_directory / "baseline_correction_tic_smooth.dat",minutes=True)
tic2.write(output_directory / "baseline_correction_tic_smooth_bc.dat",minutes=True)


# ## Tophat Baseline correction on an Intensity Matrix object
#
# The |tophat()| function acts on a single |IonChromatogram|. To perform baseline
# correction on an |IntensityMatrix| object (i.e. on all `Ion Chromatograms`)
# the |tophat_im()| function may be used.
#
# Using the same value for `struct` as above, |tophat_im()| is used as follows:

# In[6]:


from pyms.TopHat import tophat_im
from pyms.IntensityMatrix import build_intensity_matrix
im = build_intensity_matrix(data)
im_base_corr = tophat_im(im, struct="1.5m")


# Write the IC for mass 73 to disk for both the original and smoothed |IntensityMatrix|:

# In[7]:


ic = im.get_ic_at_index(73)
ic_base_corr = im_base_corr.get_ic_at_index(73)

ic.write(output_directory/"baseline_correction_ic.dat",minutes=True)
ic_base_corr.write(output_directory/"baseline_correction_ic_base_corr.dat",minutes=True)


#

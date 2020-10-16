#!/usr/bin/env python
# coding: utf-8

# # Example: IonChromatogram Objects
#
# First, setup the paths to the datafiles and the output directory, then import JCAMP_reader and build_intensity_matrix.

# In[1]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix


# Read the raw data files and create the IntensityMatrix.

# In[2]:


jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
im = build_intensity_matrix(data)


# An |IonChromatogram| object is a
# one dimensional vector containing mass intensities as a function of
# retention time. This can can be either |m/z| channel intensities
# (for example, the ion chromatogram at 73 |m/z|), or cumulative
# intensities over all measured |m/z| (TIC).
#
# An |IonChromatogram| object for the
# TIC can be obtained as follows:

# In[3]:


print(data.tic)


# The |IonChromatogram| at index 0 can be obtained with:

# In[4]:


im.get_ic_at_index(0)


# The |IonChromatogram| for the closest mass to 73 can be obtained with:

# In[5]:


im.get_ic_at_mass(73)


# An ion chromatogram object has a method |is_tic()| which returns ``True`` if
# the ion chromatogram is a TIC, ``False`` otherwise.

# In[6]:


data.tic.is_tic()


# In[7]:


im.get_ic_at_mass(73).is_tic()

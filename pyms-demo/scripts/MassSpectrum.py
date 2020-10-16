#!/usr/bin/env python
# coding: utf-8

# # Example: MassSpectrum Objects
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


# A |MassSpectrum| object contains two attributes, |mass_list| and
# |intensity_list|, a list of mass values and corresponding intensities,
# respectively. A |MassSpectrum| is returned by the |IntensityMatrix| method
# |get_ms_at_index(index)|.
#
# For example, the properties of the first |MassSpectrum| object can be obtained
# as follows:

# In[7]:


ms = im.get_ms_at_index(0)

print(ms)


# In[4]:


len(ms)


# In[5]:


len(ms.mass_list)


# In[6]:


len(ms.intensity_list)


# The length of all attributes should be the same.

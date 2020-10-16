#!/usr/bin/env python
# coding: utf-8

# # Example: IntensityMatrix Resizing
#
# Once an IntensityMatrix has been constructed from the raw GC-MS data, the
# entries of the matrix can be modified. These modifications can operate on the
# entire matrix, or individual masses or scans.
#
# First, setup the paths to the datafiles and the output directory, then import JCAMP_reader and build_intensity_matrix.

# In[26]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix


# Read the raw data files and create the IntensityMatrix.

# In[27]:


jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
im = build_intensity_matrix(data)


# ## Retention time range
#
# A basic operation on the GC-MS data is to select a specific time range for
# processing. In PyMassSpec, any data outside the chosen time range is discarded.
# The |trim()| method operates on the raw data, so any subsequent processing only
# refers to the trimmed data.
#
# The data can be trimmed to specific scans:

# In[28]:


data.trim(1000, 2000)
data.info()


# or specific retention times (in ``seconds`` or ``minutes``):

# In[29]:


data.trim("700s", "15m")
data.info()


# ## Mass Spectrum range and entries
#
# An |IntensityMatrix| object has a set mass range and interval that is derived
# from the data at the time of building the intensity matrix. The range of mass
# values can be cropped. This is done, primarily, to ensure that the range of
# masses used are consistent when comparing samples.
#
# The mass range of the intensity matrix can be "cropped" to a new (smaller)
# range as follows:

# In[30]:


im.crop_mass(60, 400)
print(im.min_mass)


# In[31]:


print(im.max_mass)


# It is also possible to set all intensities for a given mass to zero. This is
# useful for ignoring masses associated with sample preparation. The mass can be
# "nulled" with:

# In[32]:


im.null_mass(73)
sum(im.get_ic_at_mass(73).intensity_array)


# As expected, the sum of the intensity array is `0`

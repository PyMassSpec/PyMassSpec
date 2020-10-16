#!/usr/bin/env python
# coding: utf-8

# # Example: Comparing two GC-MS data sets
#
# Occasionally it is useful to compare two data sets. For example,
# one may want to check the consistency between the data set
# exported in netCDF format from the manufacturer's software, and
# the JCAMP format exported from a third party software.
#
#
#
# First, setup the paths to the datafiles and the output directory, then import JCAMP_reader and AMDI_reader.

# In[10]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.GCMS.IO.ANDI import ANDI_reader


# Then the raw data is read as before.

# In[11]:


andi_file = data_directory / "gc01_0812_066.cdf"
data1 = ANDI_reader(andi_file)
print(data1)


# In[12]:


jcamp_file = data_directory / "gc01_0812_066.jdx"
data2 = JCAMP_reader(jcamp_file)
print(data2)


# To compare the two data sets, use the function |diff()|

# In[13]:


from pyms.GCMS.Function import diff

diff(data1, data2)


# If the data cannot be compared, for example because of
# different number of scans, or inconsistent number of m/z values
# in between two scans, |diff()|
# will report the difference. For example:

# In[14]:


data2.trim(begin=1000, end=2000)


# In[15]:


diff(data1, data2)

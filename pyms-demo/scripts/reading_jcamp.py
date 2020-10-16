#!/usr/bin/env python
# coding: utf-8

# # Example: Reading JCAMP GC-MS data
#
# The PyMS package |pyms.GCMS.IO.JCAMP| provides capabilities to read the raw
# GC-MS data stored in the JCAMP-DX format.
#
#
# First, setup the paths to the datafile and the output directory,
# then import JCAMP_reader.

# In[1]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader


# Read the raw JCAMP-dx data.

# In[2]:


jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
print(data)


# ### A GCMS_data Object
#
# The object ``data`` (from the two previous examples) stores the raw data as a
# |pyms.GCMS.Class.GCMS_data| object.
# Within the |GCMS_data|
# object, raw data are stored as a list of
# |pyms.Spectrum.Scan| objects and a list of
# retention times. There are several methods available to access data and
# attributes of the |GCMS_data|
# and |Scan| objects.
#
# The |GCMS_data| object's methods relate to the raw data.
# The main properties relate to the masses, retention times and scans. For example, the
# minimum and maximum mass from all of the raw data can be returned by the
# following:

# In[3]:


print(data.min_mass)


# In[4]:


print(data.max_mass)


# A list of the first 10 retention times can be returned with:

# In[5]:


print(data.time_list[:10])


# The index of a specific retention time (in seconds) can be returned with:
#

# In[6]:


data.get_index_at_time(400.0)


# Note that this returns the index of the retention time in the
# data closest to the given retention time of 400.0 seconds.
#
# The |GCMS_data.tic| attribute
# returns a total ion chromatogram (TIC) of the data
# as an |IonChromatogram| object:

# In[7]:


print(data.tic)


# The |IonChromatogram| object is explained in a later example.
#
# ### A Scan Object
#
# A |pyms.Spectrum.Scan| object contains a list of masses and a corresponding list of intensity values
# from a single mass-spectrum scan in the raw data. Typically only non-zero (or non-threshold) intensities and
# corresponding masses are stored in the raw data.
#
# A list of the first 10 |pyms.Spectrum.Scan| objects can be returned with:

# In[8]:


scans = data.scan_list
print(scans[:10])


# A list of the first 10 masses in a scan (e.g. the 1st scan) is returned with:

# In[9]:


print(scans[0].mass_list[:10])


# A list of the first 10 corresponding intensities in a scan is returned with:

# In[10]:


print(scans[0].intensity_list[0])


# The minimum and maximum mass in an individual scan (e.g. the 1st scan) are
# returned with:

# In[11]:


print(scans[0].min_mass)


# In[12]:


print(scans[0].max_mass)


# ### Exporting data and obtaining information about a data set
#
# Often it is of interest to find out some basic information about the
# data set, e.g. the number of scans, the retention time range, and
# m/z range and so on. The |GCMS_data|
# class provides a method |info()|
# that can be used for this purpose.
#
# The entire raw data of a |GCMS_data| object can be exported to a file
# with the method |write()|:

# In[13]:


data.write(output_directory / "data")


# This method takes the filename ("output/data", in this example)
# and writes two CSV files. One has extension ".I.csv" and
# contains the intensities ("output/data.I.csv" in this example),
# and the other has the extension ".mz" and contains the
# corresponding table of m/z value ("output/data.mz.csv" in
# this example). In general, these are not two-dimensional matrices,
# because different scans may have different number of m/z
# values recorded.

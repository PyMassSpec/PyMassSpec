#!/usr/bin/env python
# coding: utf-8

# ## Example: Building an Intensity Matrix
#
# First, setup the paths to the datafiles and the output directory, then import JCAMP_reader.

# In[21]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.GCMS.IO.JCAMP import JCAMP_reader


# Read the raw data files.

# In[22]:


jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
print(data)


# Then the data can be converted to an |IntensityMatrix| using the function
# |build_intensity_matrix| from |pyms.IntensityMatrix|.
#
# The default operation of |build_intensity_matrix| is to use a bin interval of
#  one and treat the masses as floating point numbers. The default intensity
#  matrix can be built as follows:

# In[23]:


from pyms.IntensityMatrix import build_intensity_matrix

im = build_intensity_matrix(data)

print(im)


# The size as the number of scans and the number of bins can be returned with:

# In[24]:


print(im.size)


# There are 9865 scans and 551 bins in this example.
#
# The raw masses have been binned into new mass units based on the minimum mass
# in the raw data and the bin size. A list of the first ten new masses can be
# obtained as follows:

# In[25]:


print(im.mass_list[:10])


# The attributes |im.min_mass| and |im.max_mass| return the minimum and maximum mass:

# In[26]:


print(im.min_mass)


# In[27]:


print(im.max_mass)


# It is also possible to search for a particular mass, by finding the index of
# the binned mass closest to the desired mass. For example, the index of the
# closest binned mass to a mass of 73.3 |m/z| can be found by using the
# methods |im.get_index_of_mass()|:

# In[28]:


index = im.get_index_of_mass(73.3)

print(index)


# The value of the closest mass can be returned by the method
# |im.get_mass_at_index()|:

# In[29]:


im.get_mass_at_index(index)


# A mass of 73.0 is returned in this example.
#
# ## Build intensity matrix parameters
#
# The bin interval can be set to values other than one, and binning boundaries
# can also be adjusted. In the example below, to fit the 0.5 bin interval, the
# upper and lower boundaries are set to ± 0.25.

# In[30]:


im = build_intensity_matrix(data, 0.5, 0.25, 0.25)

print(im)


# The size of the intensity matrix will reflect the change in the number of bins:

# In[31]:


print(im.size)


# In[32]:


print(im.mass_list[:10])


# In this example there are 9865 scans (as before), but 1101 bins.
#
# The index and binned mass of the mass closest to 73.3 should also reflect the
# different binning.

# In[33]:


index = im.get_index_of_mass(73.3)

print(index)


# In[34]:


im.get_mass_at_index(index)


# ## Build integer mass intensity matrix
#
# It is also possible to build an intensity matrix with integer masses and a bin
# interval of one using |build_intensity_matrix_i()|. The default range for the
# binning is -0.3 and +0.7 mass units. The function is imported from
# |pyms.IntensityMatrix|:

# In[35]:


from pyms.IntensityMatrix import build_intensity_matrix_i

im = build_intensity_matrix_i(data)
print(im)


# In[36]:


print(im.size)


# In[37]:


print(im.mass_list[:10])


# The masses are now integers.

# In[38]:


index = im.get_index_of_mass(73.3)
print(index)


# In[39]:


im.get_mass_at_index(index)


# The lower and upper bounds can be adjusted with |build_intensity_matrix_i(data, lower, upper)|.

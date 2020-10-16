"""proc.py
"""
# This file has been replaced by jupyter/IntensityMatrix_Resizing.ipynb

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i

# read the raw data
andi_file = "data/gc01_0812_066.cdf"

data = ANDI_reader(andi_file)

# original mass range
print(data.min_mass, data.max_mass)

# Build an integer mass intensity matrix
im = build_intensity_matrix_i(data)

# Crop the mass range to 60 to 400 m/z
im.crop_mass(60, 400)

# new mass range
print(im.min_mass, im.max_mass)

# Set intensities of masses related to solvents to null (zero)
im.null_mass(73)
im.null_mass(147)

# test if all intensities zeroed
print(sum(im.get_ic_at_mass(73).intensity_array))
print(sum(im.get_ic_at_mass(147).intensity_array))

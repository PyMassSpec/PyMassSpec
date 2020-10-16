"""proc.py
"""
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix
from pyms.Utils.IO import save_data

# read the raw data as a GCMS_data object
jcamp_file = "data/gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)

# IntensityMatrix
# must build intensity matrix before accessing any intensity matrix methods.

# default, float masses with interval (bin interval) of one from min mass
print("default intensity matrix, bin interval = 1, boundary +/- 0.5")
im = build_intensity_matrix(data)

#
# Saving data
#

# save the intensity matrix values to a file
mat = im.intensity_array
print("saving intensity matrix intensity values...")
save_data("output/im.dat", mat)

# Export the entire IntensityMatrix as CSV. This will create
# data.im.csv, data.mz.csv, and data.rt.csv where
# these are the intensity matrix, retention time
# vector, and m/z vector in the CSV format
print("exporting intensity matrix data...")
im.export_ascii("output/data")

# Export the entire IntensityMatrix as LECO CSV. This is
# useful for import into AnalyzerPro
print("exporting intensity matrix data to LECO CSV format...")
im.export_leco_csv("output/data_leco.csv")

#
# Import saved data
#

from pyms.IntensityMatrix import import_leco_csv

# import LECO CSV file
print("importing intensity matrix data from LECO CSV format...")
iim = import_leco_csv("output/data_leco.csv")

# Check size to original
print("Output dimensions:", im.size, " Input dimensions:", iim.size)

"""proc.py
"""

# TODO: mzML demo; need example mzML file

import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

from pyms.GCMS.IO.MZML import mzML_reader


# read the raw data
mzml_file = data_directory / ".mzML"
data = mzML_reader(mzml_file)
print(data)

# raw data operations
print("minimum mass found in all data: ", data.min_mass)
print("maximum mass found in all data: ", data.max_mass)

# time
time = data.time_list
print(time)
print("number of retention times: ", len(time))
print("retention time of 1st scan: ", time[0], "sec")
print("index of 400sec in time_list: ", data.get_index_at_time(400.0))

# TIC
tic = data.tic
print(tic)
print("number of scans in TIC: ", len(tic))
print("start time of TIC: ", tic.get_time_at_index(0), "sec")

# raw scans
scans = data.scan_list
print(scans)
print(scans[0].mass_list)
print("1st mass value for 1st scan: ", scans[0].mass_list[0])
print("1st intensity value for 1st scan: ", scans[0].intensity_list[0])

print("minimum mass found in 1st scan: ", scans[0].min_mass)
print("maximum mass found in 1st scan: ", scans[0].max_mass)

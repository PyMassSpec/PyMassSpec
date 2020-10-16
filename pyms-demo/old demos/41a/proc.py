"""proc.py
"""
# This file has been replaced by jupyter/NoiseSmoothing.ipynb

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.Noise.Window import window_smooth

# read the raw data
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# get the TIC
tic = data.tic

# apply window smoothing: mean and median, in both cases
# the window is 5 points
tic1 = window_smooth(tic, window=5)
tic2 = window_smooth(tic, window=5, median=True)

# an example of how to specify window as a time string
# (7 seconds in this case)
tic3 = window_smooth(tic, window='7s')

# write the origina IC and smoothed ICs to a file
tic.write("output/tic.dat",minutes=True)
tic1.write("output/tic1.dat",minutes=True)
tic2.write("output/tic2.dat",minutes=True)

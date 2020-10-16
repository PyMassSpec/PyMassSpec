"""proc.py
"""
# This file has been replaced by jupyter/BaselineCorrection.ipynb

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat

# read the raw data
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# get the TIC
tic = data.tic

# apply noise smoothing and baseline correction
tic1 = savitzky_golay(tic)
tic2 = tophat(tic1, struct="1.5m")

# save smoothed/baseline corrected TIC
tic.write("output/tic.dat", minutes=True)
tic1.write("output/tic_smooth.dat", minutes=True)
tic2.write("output/tic_smooth_bc.dat", minutes=True)

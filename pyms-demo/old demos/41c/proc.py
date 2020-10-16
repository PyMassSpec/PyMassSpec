"""proc.py
"""
# This file has been replaced by jupyter/NoiseSmoothing.ipynb

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.Noise.SavitzkyGolay import savitzky_golay

# read the raw data
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# get the TIC
tic = data.tic

tic1 = savitzky_golay(tic)

tic.write("output/tic.dat", minutes=True)
tic1.write("output/tic1.dat",minutes=True)

"""proc.py
"""
# This file has been replaced by jupyter/NoiseSmoothing.ipynb

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix
from pyms.Noise.SavitzkyGolay import savitzky_golay_im

# read the raw data
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# build an intensity matrix object from the data
im = build_intensity_matrix(data)

# Use Savitzky-Golay filtering to smooth all IC's in the IM
print("Smoothing ...")
im_smooth = savitzky_golay_im(im)
print("Done")

# find the IC for derivatisation product ion before smoothing
ic = im.get_ic_at_index(73)

# find the IC for derivatisation product ion after smoothing
ic_smooth = im_smooth.get_ic_at_index(73)

ic.write("output/ic.dat",minutes=True)
ic_smooth.write("output/ic_smooth.dat",minutes=True)

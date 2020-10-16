"""proc.py
"""

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix
from pyms.TopHat import tophat_im

# This file has been replaced by jupyter/BaselineCorrection.ipynb

# read the raw data
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# build an intensity matrix object from the data
im = build_intensity_matrix(data)

# Use TopHat baseline correction on all IC's in the IM
print("Smoothing ...")
im_base_corr = tophat_im(im, struct="1.5m")
print("Done")

# find the IC for derivatisation product ion before smoothing
ic = im.get_ic_at_index(73)

# find the IC for derivatisation product ion after smoothing
ic_base_corr = im_base_corr.get_ic_at_index(73)

ic.write("output/ic.dat",minutes=True)
ic_base_corr.write("output/ic_smooth.dat",minutes=True)

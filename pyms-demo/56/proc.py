"""proc.py
"""

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat

from pyms.BillerBiemann import BillerBiemann, rel_threshold, num_ions_threshold

from pyms.Peak.Function import peak_top_ion_areas

# read the raw data as a GCMS_data object
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

im = build_intensity_matrix_i(data)

n_scan, n_mz = im.size

print("Intensity matrix size (scans, masses):", (n_scan, n_mz))

# noise filter and baseline correct
for ii in range(n_mz):
    ic = im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    im.set_ic_at_index(ii, ic_bc)

# Use Biller and Biemann technique to find apexing ions at a scan.
peak_list = BillerBiemann(im, points=9, scans=2)

# percentage ratio of ion intensity to max ion intensity
r = 2
# minimum number of ions, n
n = 3
# greater than or equal to threshold, t
t = 10000

# trim by relative intensity
pl = rel_threshold(peak_list, r)

# trim by threshold
new_peak_list = num_ions_threshold(pl, n, t)

print("Number of filtered peaks: ", len(new_peak_list))

# find and set areas
print("Top 5 most abundant ions for each peak ")

for peak in new_peak_list:
    rt = peak.rt
    # Only test interesting sub-set from 29.5 to 32.5 minutes
    if rt >= 29.5*60.0 and rt <= 32.5*60.0:
        # determine and set ion areas, use default num of ions =5
        areas_dict = peak_top_ion_areas(im, peak)
        peak.ion_areas = areas_dict

        area_dict = peak.ion_areas
        # print the top 5 ions for each peak
        print(area_dict.keys())

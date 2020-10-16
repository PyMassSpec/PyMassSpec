import pathlib

from pyms.BillerBiemann import BillerBiemann, num_ions_threshold, rel_threshold
from pyms.Display import Display
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import peak_sum_area
from pyms.TopHat import tophat

data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

# read raw data
andi_file = data_directory / "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

data.trim(4101, 4350)

# Build Intensity Matrix
real_im = build_intensity_matrix_i(data)

n_scan, n_mz = real_im.size

 # perform necessary pre filtering
for ii in range(n_mz):
    ic = real_im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    real_im.set_ic_at_index(ii, ic_bc)


 # Detect Peaks
peak_list = BillerBiemann(real_im, points=3, scans=2)

print("Number of peaks found in real data: ", len(peak_list))

######### Filter peaks###############
# Filter the peak list,
# first by removing all intensities in a peak less than a given relative
# threshold,
# then by removing all peaks that have less than a given number of ions above
# a given value

# Parameters
# percentage ratio of ion intensity to max ion intensity
r = 1

# minimum number of ions, n
n = 3
# greater than or equal to threshold, t
t = 10000

# trim by relative intensity
pl = rel_threshold(peak_list, r)

# trim by threshold
real_peak_list = num_ions_threshold(pl, n, t)
print("Number of filtered peaks in real data: ", len(real_peak_list))

# Set the peak areas
for peak in real_peak_list:
    area = peak_sum_area(real_im, peak)
    peak.area = area


# real_peak_list is PyMassSpec' best guess at the true peak list

################## Run Simulator ######################
# Simulator takes a peak list, time_list and mass_list
# and returns an IntensityMatrix object.
# The mass_list and time_list are the same for the real
# data and the simulated data.

time_list = real_im.time_list
mass_list = real_im.mass_list

sim_im = gcms_sim(time_list, mass_list, real_peak_list)
# sim_im is an IntensityMatrix object

# Now add noise to the simulated intensity matrix object
scale = 1000
cutoff = 10000
prop = 0.0003
add_gaussv_noise(sim_im, scale, cutoff, prop)


### Now display the ics from the simulated data
ics = []
for i in range(n_mz):
    ics.append(sim_im.get_ic_at_index(i))

display = Display()
for ic in ics:
    display.plot_ic(ic)
display.do_plotting('ICs, and PyMassSpec Detected Peaks of Simulated Data')
display.show_chart()

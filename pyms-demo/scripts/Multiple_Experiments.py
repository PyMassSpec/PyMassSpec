#!/usr/bin/env python
# coding: utf-8

# ## Example: Creating Multiple Experiments
#
# In example three GC-MS experiments are prepared for peak alignment. The
# experiments are named ``a0806_077``, ``a0806_078``, ``a0806_079``, and
# represent separate GC-MS sample runs from the same biological sample.
#
# The procedure is the same as for the previous example, but is repeated three times.
#
# First, setup the paths to the datafiles and the output directory, then import the required functions.

# In[11]:


import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

from pyms.BillerBiemann import BillerBiemann, num_ions_threshold, rel_threshold
from pyms.Experiment import Experiment
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas
from pyms.TopHat import tophat


# Define the data files to process

# In[12]:


expr_codes = ["a0806_077", "a0806_078", "a0806_079"]
# expr_codes = ["a0806_140", "a0806_141", "a0806_142"]


# Loop over the experiments and perform the processing.

# In[13]:


for expr_code in expr_codes:

	print(f" -> Processing experiment '{expr_code}'")

	andi_file = data_directory / f"{expr_code}.cdf"

	data = ANDI_reader(andi_file)

	im = build_intensity_matrix_i(data)

	n_scan, n_mz = im.size

	# Preprocess the data (Savitzky-Golay smoothing and Tophat baseline detection)

	for ii in range(n_mz):
		ic = im.get_ic_at_index(ii)
		ic1 = savitzky_golay(ic)
		ic_smooth = savitzky_golay(ic1)  # Why the second pass here?
		ic_bc = tophat(ic_smooth, struct="1.5m")
		im.set_ic_at_index(ii, ic_bc)

	# Peak detection
	pl = BillerBiemann(im, points=9, scans=2)

	# Trim the peak list by relative intensity
	apl = rel_threshold(pl, percent=2)

	# Trim the peak list by noise threshold
	peak_list = num_ions_threshold(apl, n=3, cutoff=3000)

	print("\t -> Number of Peaks found:", len(peak_list))

	print("\t -> Executing peak post-processing and quantification...")

	# Set the mass range, remove unwanted ions and estimate the peak area
	# For peak alignment, all experiments must have the same mass range

	for peak in peak_list:
		peak.crop_mass(51, 540)

		peak.null_mass(73)
		peak.null_mass(147)

		area = peak_sum_area(im, peak)
		peak.area = area
		area_dict = peak_top_ion_areas(im, peak)
		peak.ion_areas = area_dict

	# Create an Experiment
	expr = Experiment(expr_code, peak_list)

	# Use the same retention time range for all experiments
	lo_rt_limit = "6.5m"
	hi_rt_limit = "21m"

	print(f"\t -> Selecting retention time range between '{lo_rt_limit}' and '{hi_rt_limit}'")

	expr.sele_rt_range([lo_rt_limit, hi_rt_limit])

	# Save the experiment to disk.
	output_file = output_directory / "experiments" / f"{expr_code}.expr"
	print(f"\t -> Saving the result as '{output_file}'")

	expr.dump(output_file)


# The previous set of data all belong to the same experimental condition. That is,
# they represent one group and any comparison between the data is a within group
# comparison. For the original experiment, another set of GC-MS data was collected
# for a different experimental condition. This group must also be stored as a set
# of experiments, and can be used for between group comparison.
#
# The second set of data files are named ``a0806_140``, ``a0806_141``, and ``a0806_142``, and are
# processed and stored as above.
#
# In the example notebook, you can uncomment the line in code cell 2 and run the
# notebook again to process the second set of data files.

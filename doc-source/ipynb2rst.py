"""
Convert ipynb notebook to rst
"""

# stdlib
import io
import pathlib

# 3rd party
# Install pandoc
from pandoc_downloader import download_pandoc

download_pandoc()

# 3rd party
# Import the RST exproter
from nbconvert import RSTExporter

# Instantiate it
rst_exporter = RSTExporter()
# from nbconvert.preprocessors import ExecutePreprocessor
# ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
# rst_exporter.register_preprocessor(ep, True)

replacements = {
		"pyms.GCMS.IO.JCAMP":
				":mod:`pyms.GCMS.IO.JCAMP`",
		"pyms.GCMS.IO.ANDI":
				":mod:`pyms.GCMS.IO.ANDI`",
		"pyms.GCMS.Class.GCMS_data":
				":class:`pyms.GCMS.Class.GCMS_data`",
		"GCMS_data":
				":class:`~pyms.GCMS.Class.GCMS_data`",
		"pyms.Spectrum.Scan":
				":class:`pyms.Spectrum.Scan`",
		"Scan":
				":class:`~pyms.Spectrum.Scan`",
		"GCMS_data.tic":
				":attr:`GCMS_data.tic <pyms.GCMS.Class.GCMS_data.tic>`",
		"IonChromatogram":
				":class:`~pyms.IonChromatogram.IonChromatogram`",
		"info()":
				":py:meth:`info() <pyms.GCMS.Class.GCMS_data.info()>`",
		"write()":
				":py:meth:`write() <pyms.GCMS.Class.GCMS_data.write()>`",
		"diff()":
				":py:meth:`diff() <pyms.GCMS.Function.diff()>`",
		"build_intensity_matrix()":
				":meth:`build_intensity_matrix() <pyms.IntensityMatrix.build_intensity_matrix>`",
		"build_intensity_matrix_i()":
				":meth:`build_intensity_matrix_i() <pyms.IntensityMatrix.build_intensity_matrix_i>`",
		"pyms.IntensityMatrix":
				":mod:`pyms.IntensityMatrix <pyms.IntensityMatrix>`",
		"IntensityMatrix":
				":class:`~pyms.IntensityMatrix.IntensityMatrix`",
		"im.min_mass":
				":attr:`im.min_mass <pyms.IntensityMatrix.IntensityMatrix.min_mass>`",
		"im.max_mass":
				":attr:`im.max_mass <pyms.IntensityMatrix.IntensityMatrix.max_mass>`",
		"im.get_index_of_mass()":
				":meth:`im.get_index_of_mass() <pyms.IntensityMatrix.IntensityMatrix.get_index_of_mass>`",
		"m/z":
				":math:`m/z`",
		"im.get_mass_at_index()":
				":meth:`im.get_mass_at_index() <pyms.IntensityMatrix.IntensityMatrix.get_mass_at_index>`",
		"build_intensity_matrix_i(data, lower, upper)":
				":meth:`build_intensity_matrix_i(data, lower, upper) <pyms.IntensityMatrix.build_intensity_matrix_i>`",
		"MassSpectrum":
				":class:`~pyms.Spectrum.MassSpectrum`",
		"mass_list":
				":attr:`~pyms.Spectrum.MassSpectrum.mass_list`",
		"intensity_list":
				":attr:`~pyms.Spectrum.MassSpectrum.intensity_list`",
		"get_ms_at_index(index)":
				":meth:`get_ms_at_index(index) <pyms.IntensityMatrix.IntensityMatrix.get_ms_at_index()>`",
		"is_tic()":
				":meth:`is_tic() <pyms.IonChromatogram.IonChromatogram.is_tic>`",
		"trim()":
				":meth:`trim() <pyms.GCMS.Class.GCMS_data.trim>`",
		"savitzky_golay()":
				":meth:`savitzky_golay() <pyms.Noise.SavitzkyGolay.savitzky_golay>`",
		"savitzky_golay_im()":
				":meth:`savitzky_golay_im() <pyms.Noise.SavitzkyGolay.savitzky_golay_im>`",
		"window_smooth()":
				":meth:`window_smooth() <pyms.Noise.Window.window_smooth>`",
		"window_smooth_im()":
				":meth:`window_smooth_im() <pyms.Noise.Window.window_smooth_im>`",
		"tophat()":
				":meth:`tophat() <pyms.TopHat.tophat>`",
		"tophat_im()":
				":meth:`tophat_im() <pyms.TopHat.tophat_im>`",
		"Peak":
				":class:`~pyms.Peak.Class.Peak`",
		"pyms.Peak.Class.Peak.rt":
				":attr:`pyms.Peak.Class.Peak.rt`",
		"pyms.Peak.Class.Peak.mass_spectrum":
				":attr:`pyms.Peak.Class.Peak.mass_spectrum`",
		"pyms.Peak.Class.Peak.UID":
				":attr:`pyms.Peak.Class.Peak.UID`",
		"crop_mass()":
				":meth:`crop_mass() <pyms.Peak.Class.Peak.crop_mass>`",
		"null_mass()":
				":meth:`null_mass() <pyms.Peak.Class.Peak.null_mass>`",
		"rel_threshold()":
				":meth:`rel_threshold() <pyms.BillerBiemann.rel_threshold>`",
		"num_ions_threshold()":
				":meth:`num_ions_threshold() <pyms.BillerBiemann.num_ions_threshold>`",
		"window_analyzer()":
				":meth:`window_analyzer() <pyms.Noise.Analysis.window_analyzer>`",
		"pyms.Peak.Class.Peak.area":
				":attr:`~pyms.Peak.Class.Peak.area`",
		"set_ion_areas()":
				":meth:`set_ion_areas() <pyms.Peak.Class.Peak.set_ion_areas>`",
		"peak_sum_area()":
				":meth:`peak_sum_area() <pyms.Peak.Function.peak_sum_area>`",
		"pyms.Peak.Function":
				":mod:`pyms.Peak.Function`",
		"Experiment":
				":class:`~pyms.Experiment.Experiment`",
		"Alignment":
				":class:`~pyms.DPA.Alignment.Alignment` ",
		"exprl2alignment()":
				":meth:`exprl2alignment() <pyms.DPA.Function.exprl2alignment>`.",
		"pyms.DPA.PairwiseAlignment.PairwiseAlignment":
				":class:`pyms.DPA.PairwiseAlignment.PairwiseAlignment`",
		"PairwiseAlignment":
				":class:`~pyms.DPA.PairwiseAlignment.PairwiseAlignment`",
		"align_with_tree()":
				":meth:`align_with_tree() <pyms.DPA.Alignment.align_with_tree>`",
		"Display":
				":class:`pyms.Display.Display`",
		"plot_ic()":
				":py:meth:`plot_ic() <pyms.Display.plot_ic>`",
		"plot_peaks()":
				":py:meth:`plot_peaks() <pyms.Display.plot_peaks>`",
		"ClickEventHandler":
				":class:`pyms.Display.ClickEventHandler`",
		"ClickEventHandler(peak_list=new_peak_list)":
				":class:`pyms.Display.ClickEventHandler`",
		}

string_replacements = {
		"developed in\nmathematical morphology": "developed in mathematical morphology [1]_",
		"proteomics based mass spectrometry": "proteomics based mass spectrometry [2]_",
		}

notebooks = [
		"reading_jcamp",
		"reading_andi",
		"comparing_datasets",
		"IntensityMatrix",
		"IntensityMatrix_Resizing",
		"IntensityMatrix_Preprocessing",
		"MassSpectrum",
		"IonChromatogram",
		"NoiseSmoothing",
		"BaselineCorrection",
		"Peak",
		"Peak_Detection",
		"Peak_Filtering_Noise_Analysis",
		"Peak_Area_Estimation",
		"Experiment",
		"Multiple_Experiments",
		"DPA",
		"Displaying_TIC",
		"Displaying_Multiple_IC",
		"Displaying_Mass_Spec",
		"Displaying_Detected_Peaks",
		"Display_User_Interaction",
		]

demo_rst_dir = pathlib.Path("./demo_rst").resolve()
if not demo_rst_dir.is_dir():
	demo_rst_dir.mkdir()

images_dir = pathlib.Path("./graphics").resolve()
if not images_dir.is_dir():
	images_dir.mkdir()

for notebook in notebooks:
	# Convert the notebook to RST format
	(body, resources) = rst_exporter.from_file(f"../pyms-demo/jupyter/{notebook}.ipynb")
	for original, replacement in replacements.items():
		body = body.replace(f"\\|{original}\\|", replacement)

		# Sometimes trailing slash doesn't get escaped
		body = body.replace(f"\\|{original}|", replacement)

		original = original.replace("_", "\\_")
		body = body.replace(f"\\|{original}\\|", replacement)

		# Sometimes trailing slash doesn't get escaped
		body = body.replace(f"\\|{original}|", replacement)

	for original, replacement in string_replacements.items():
		body = body.replace(original, replacement)

	i = 1
	while True:
		# print(i)
		# replace nbconvert syntax with nbsphinx syntax
		old_body = body
		# body = old_body.replace(".. code:: ipython2", f".. nbinput:: ipython3", 1) #
		body = old_body.replace(
				".. code:: ipython2", f".. nbinput:: ipython3\n    :execution-count: {i}", 1
				)  # \n    :execution-count: {i}
		if body == old_body:
			break
		# body = body.replace(".. parsed-literal::", f".. nboutput::", 1)
		body = body.replace(".. parsed-literal::", f".. parsed-literal::", 1)
		i += 1

	outputs = resources["outputs"]
	if outputs:
		# Embedded images present
		# Put images in the images_dir

		# Write images to file
		for name, data in outputs.items():
			# print(name)
			# print(data)
			with io.open(images_dir / f"{notebook}_{name}", 'wb') as f:
				f.write(data)

		# Replace `.. image:: output` with `.. image:: {notebook}_output`
		body = body.replace(".. image:: output", f".. image:: graphics/{notebook}_output")

	# Write rst to file
	with open(f"{demo_rst_dir/notebook}.rst", "w") as fp:
		fp.write(body)

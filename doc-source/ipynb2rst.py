"""
Convert ipynb notebook to rst

Based on https://github.com/bebraw/pypandoc
MIT Licensed
"""

# stdlib
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple
from urllib.request import urlopen

# 3rd party
from domdf_python_tools.paths import PathPlus


def get_pandoc_urls(version: str = "latest") -> Tuple[Dict[str, List[str]], str]:
	"""Get the urls of pandoc's binaries
	Uses sys.platform keys, but removes the 2 from linux2
	Adding a new platform means implementing unpacking in "DownloadPandocCommand"
	and adding the URL here

	:param version: pandoc version.
		Valid values are either a valid pandoc version e.g. "1.19.1", or "latest".

	:return: A dictionary with mapping platforms to the url of the relevant binaries,
		and the actual pandoc version (e.g. "latest" will be resolved to the actual version).
	"""

	# url to pandoc download page
	url = f"https://github.com/jgm/pandoc/releases/{('tag/' if version != 'latest' else '')}{version}"

	# read the HTML content
	response = urlopen(url)
	content = response.read()

	# regex for the binaries
	regex = re.compile(r"/jgm/pandoc/releases/download/.*\.(?:msi|deb|pkg)")

	# a list of urls to the binaries
	pandoc_urls_list = regex.findall(content.decode("utf-8"))

	# actual pandoc version
	version = pandoc_urls_list[0].split('/')[5]

	# dict that lookup the platform from binary extension
	ext2platform = {"msi": "win32", "deb": "linux", "pkg": "darwin"}

	# parse pandoc_urls from list to dict
	pandoc_urls = {ext2platform[url_frag[-3:]]: ("https://github.com" + url_frag) for url_frag in pandoc_urls_list}

	return pandoc_urls, version


def _make_executable(path):
	mode = os.stat(path).st_mode
	mode |= (mode & 0o444) >> 2  # copy R bits to X
	print("* Making %s executeable..." % (path))
	os.chmod(path, mode)


def download_pandoc(url: Optional[str] = None, targetfolder: str = "~/bin", version: str = "latest"):
	"""
	Download and unpack pandoc.

	Downloads prebuild binaries for pandoc from ``url`` and unpacks it into ``targetfolder``.

	:param url: URL for the to be downloaded pandoc binary distribution for
		the platform under which this python runs.

	:param str targetfolder: directory, where the binaries should be installed
		to. If no `targetfolder` is give, uses a platform specific user
		location: `~/bin` on Linux, `~/Applications/pandoc` on Mac OS X, and
		`~\\AppData\\Local\\Pandoc` on Windows.
	"""

	pandoc_urls, _ = get_pandoc_urls(version)

	if url is None:
		url = pandoc_urls["linux"]

	filename = url.split('/')[-1]
	if os.path.isfile(filename):
		print(f"* Using already downloaded file {filename}")
	else:
		print(f"* Downloading pandoc from {url} ...")

		# https://stackoverflow.com/questions/30627937/tracebaclk-attributeerroraddinfourl-instance-has-no-attribute-exit
		response = urlopen(url)
		with open(filename, "wb") as out_file:
			shutil.copyfileobj(response, out_file)

	targetfolder = PathPlus(targetfolder).expanduser()

	# Make sure target folder exists...
	if not targetfolder.is_dir():
		targetfolder.mkdir(parents=True)

	print(f"* Unpacking {filename} to tempfolder...")

	tempfolder = PathPlus(tempfile.mkdtemp())
	cur_wd = os.getcwd()
	filename = os.path.abspath(filename)
	try:
		os.chdir(tempfolder)
		cmd = ["ar", 'x', filename]
		subprocess.check_call(cmd)

		dir_listing = set(os.listdir(tempfolder))
		if "data.tar.gz" in dir_listing:
			cmd = ["tar", "xzf", "data.tar.gz"]
		elif "data.tar.xz" in dir_listing:
			cmd = ["tar", "xJf", "data.tar.xz"]
		elif "data.tar.bz" in dir_listing:
			cmd = ["tar", "xjf", "data.tar.bz"]
		else:
			raise FileNotFoundError(f"'data' archive not found. Files in the download are:\n{dir_listing}")

		subprocess.check_call(cmd)
		src = tempfolder / "usr" / "bin" / "pandoc"
		dst = targetfolder / "pandoc"

		if not dst.parent.is_dir():
			dst.parent.mkdir(parents=True)

		print(f"* Copying 'pandoc' to {str(targetfolder)!r} ...")
		shutil.copyfile(src, dst)
		_make_executable(dst)

		src = tempfolder / "usr" / "share" / "doc" / "pandoc" / "copyright"
		dst = targetfolder / "copyright.pandoc"
		print(f"* Copying copyright to {str(targetfolder)!r} ...")
		shutil.copyfile(src, dst)
	finally:
		os.chdir(cur_wd)
		shutil.rmtree(tempfolder)


# download_pandoc()

# Import the RST exporter and instantiate it

# 3rd party
from nbconvert import RSTExporter

rst_exporter = RSTExporter()

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
				":meth:`exprl2alignment() <pyms.DPA.Alignment.exprl2alignment>`.",
		"pyms.DPA.PairwiseAlignment.PairwiseAlignment":
				":class:`pyms.DPA.PairwiseAlignment.PairwiseAlignment`",
		"PairwiseAlignment":
				":class:`~pyms.DPA.PairwiseAlignment.PairwiseAlignment`",
		"align_with_tree()":
				":meth:`align_with_tree() <pyms.DPA.PairwiseAlignment.align_with_tree>`",
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

demo_rst_dir = PathPlus("./demo_rst").resolve()
if not demo_rst_dir.is_dir():
	demo_rst_dir.mkdir()

images_dir = PathPlus("./graphics").resolve()
if not images_dir.is_dir():
	images_dir.mkdir()

for notebook in notebooks:
	# Convert the notebook to RST format
	(body, resources) = rst_exporter.from_file(f"../pyms-demo/jupyter/{notebook}.ipynb")
	for original, replacement in replacements.items():
		body = body.replace(f"\\|{original}\\|", replacement)

		# Sometimes trailing slashes doesn't get escaped
		body = body.replace(f"\\|{original}|", replacement)

		original = original.replace('_', "\\_")
		body = body.replace(f"\\|{original}\\|", replacement)

		# Sometimes trailing slashes doesn't get escaped
		body = body.replace(f"\\|{original}|", replacement)

	for original, replacement in string_replacements.items():
		body = body.replace(original, replacement)

	i = 1
	while True:
		old_body = body

		# replace nbconvert syntax with nbsphinx syntax
		body = re.sub(r".. code:: ipython\d?", f".. nbinput:: ipython3\n    :execution-count: {i}", body, 1)
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
			(images_dir / f"{notebook}_{name}").write_bytes(data)

		# Replace `.. image:: output` with `.. image:: {notebook}_output`
		body = body.replace(".. image:: output", f".. image:: graphics/{notebook}_output")

	# Write rst to file
	PathPlus(f"{demo_rst_dir/notebook}.rst").write_clean(body)

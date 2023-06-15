# 3rd party
from docutils import nodes
from sphinx.application import Sphinx
from sphinx.errors import NoUri


def handle_missing_xref(app: Sphinx, env, node: nodes.Node, contnode: nodes.Node) -> None:
	if node.get("reftarget", '') in {
			"pyms.Experiment.Experiment.get_expr_code",
			"pyms.Experiment.Experiment.get_peak_list",
			"pyms.Experiment.Experiment.store",
			"pyms.Experiment.store_expr",
			"pyms.GCMS.Class.GCMS_data.get_scan_list",
			"pyms.GCMS.Class.GCMS_data.get_tic",
			"pyms.Gapfill.Class.MissingPeak.get_common_ion",
			"pyms.Gapfill.Class.MissingPeak.get_common_ion_area",
			"pyms.Gapfill.Class.MissingPeak.get_exact_rt",
			"pyms.Gapfill.Class.MissingPeak.get_qual_ion1",
			"pyms.Gapfill.Class.MissingPeak.get_qual_ion2",
			"pyms.Gapfill.Class.MissingPeak.get_rt",
			"pyms.Gapfill.Class.MissingPeak.set_common_ion_area",
			"pyms.Gapfill.Class.MissingPeak.set_exact_rt",
			"pyms.Gapfill.Class.Sample.get_missing_peaks",
			"pyms.Gapfill.Class.Sample.get_mp_rt_area_dict",
			"pyms.Gapfill.Class.Sample.get_name",
			"pyms.Gapfill.Function.transposed",
			"pyms.IonChromatogram.IonChromatogram.get_mass",
			"pyms.IonChromatogram.IonChromatogram.get_time_step",
			"pyms.IonChromatogram.IonChromatogram.set_intensity_array",
			"pyms.Mixins.MaxMinMassMixin.get_max_mass",
			"pyms.Mixins.MaxMinMassMixin.get_min_mass",
			"pyms.Mixins.MassListMixin.get_mass_list",
			"pyms.Mixins.TimeListMixin.get_time_list",
			"pyms.Mixins.IntensityArrayMixin.get_intensity_array",
			"pyms.Mixins.IntensityArrayMixin.get_matrix_list",
			"pyms.Peak.Class.Peak.get_area",
			"pyms.Peak.Class.Peak.get_ic_mass",
			"pyms.Peak.Class.Peak.get_ion_areas",
			"pyms.Peak.Class.Peak.get_mass_spectrum",
			"pyms.Peak.Class.Peak.get_pt_bounds",
			"pyms.Peak.Class.Peak.get_rt",
			"pyms.Peak.Class.Peak.get_UID",
			"pyms.Peak.Class.Peak.set_area",
			"pyms.Peak.Class.Peak.set_ic_mass",
			"pyms.Peak.Class.Peak.set_ion_areas",
			"pyms.Peak.Class.Peak.set_mass_spectrum",
			"pyms.Peak.Class.Peak.set_pt_bounds",
			"pyms.Peak.Class.Peak.pt_bounds",
			"pyms.Utils.Utils.is_positive_int",
			"pyms.Utils.Utils.is_list_of_dec_nums",
			"pyms.Gapfill.Function.file2matrix",
			"pyms.IonChromatogram.IonChromatogram`",
			}:
		raise NoUri


def setup(app: Sphinx):
	app.connect("missing-reference", handle_missing_xref, priority=950)

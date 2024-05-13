# stdlib
import json

expected_order = [
		"tests.test_Display.test_plot_mass_spec_width",
		"tests.test_Display.test_plot_ic",
		"tests.test_Display.test_plot_tic",
		"tests.test_Display.test_plot_ic_title",
		"tests.test_Display.test_plot_ic_linestyle",
		"tests.test_Display.test_plot_tic_linewidth",
		"tests.test_Display.test_plot_mass_spec",
		"tests.test_Display.test_plot_ic_multiple",
		"tests.test_Display.test_plot_mass_spec_alpha",
		"tests.test_Display.test_plot_tic_alpha",
		"tests.test_Display.test_plot_tic_title",
		"tests.test_Display.test_plot_tic_linestyle",
		"tests.test_Display.test_plot_mass_spec_linestyle",
		"tests.test_Display.test_plot_ic_label",
		"tests.test_Display.test_plot_tic_label",
		"tests.test_Display.test_plot_ic_alpha",
		"tests.test_Display.test_plot_ic_linewidth",
		"tests.test_Display.test_plot_mass_spec_title",
		]

for filename in [
		"tests/image_hashes_36.json",
		"tests/image_hashes_37.json",
		"tests/image_hashes_38.json",
		"tests/image_hashes.json",
		]:
	with open(filename, encoding="UTF-8") as fp:
		data = list(json.load(fp).items())

	data.sort(key=lambda x: expected_order.index(x[0]))
	with open(filename, 'w', encoding="UTF-8") as fp:
		json.dump(dict(data), fp, indent=2)
		fp.write('\n')

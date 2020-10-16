# stdlib
import tempfile

# 3rd party
from domdf_python_tools.paths import PathPlus
from pytest_regressions.file_regression import FileRegressionFixture

# this package
from pyms.Gapfill.Function import file2dataframe


def test_file2dataframe(file_regression: FileRegressionFixture):
	area_file = PathPlus(__file__).parent / "area.csv"

	with tempfile.TemporaryDirectory() as tmpdir:
		tmpdir_p = PathPlus(tmpdir)
		file2dataframe(area_file).to_csv(tmpdir_p / "area.csv", index=False, na_rep="NA")

		file_regression.check(PathPlus(tmpdir_p / "area.csv").read_text(), encoding="UTF-8", extension=".csv")

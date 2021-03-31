# 3rd party
from coincidence.regressions import AdvancedFileRegressionFixture
from domdf_python_tools.paths import PathPlus

# this package
from pyms.Gapfill.Function import file2dataframe


class MaxPrecisionFloatFormat(str):
	__slots__ = ()

	def __new__(cls, max_precision: int):
		return super().__new__(cls, f"%.{max_precision}f")

	def __mod__(self, other):
		modded_string = super().__mod__(other).rstrip('0')

		if modded_string.endswith('.'):
			modded_string += '0'

		if modded_string == "-0.0":
			return "0.0"

		return modded_string


def test_file2dataframe(tmp_pathplus: PathPlus, advanced_file_regression: AdvancedFileRegressionFixture):
	area_file = PathPlus(__file__).parent / "area.csv"

	file2dataframe(area_file).to_csv(
			tmp_pathplus / "area.csv",
			index=False,
			na_rep="NA",
			float_format=MaxPrecisionFloatFormat(3),
			)

	advanced_file_regression.check_file(tmp_pathplus / "area.csv")

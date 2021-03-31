# stdlib
import math
import random
import statistics
from decimal import Decimal
from fractions import Fraction
from typing import List

# 3rd party
import pytest

# this package
from pyms.Utils import Math


class TestMean:

	def test_torture_pep(self):
		# "Torture Test" from PEP-450.
		assert statistics.mean([1e100, 1, 3, -1e100]) == 1
		assert Math.mean([1e100, 1, 3, -1e100]) == 1
		assert statistics.mean([1e100, 1, 3, -1e100]) == Math.mean([1e100, 1, 3, -1e100])

	def test_ints(self):
		# Test mean with ints.
		data = [0, 1, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 7, 7, 8, 9]
		random.shuffle(data)
		assert statistics.mean(data) == 4.8125
		assert Math.mean(data) == 4.8125
		assert statistics.mean(data) == Math.mean(data)

	def test_floats(self):
		# Test mean with floats.
		data = [17.25, 19.75, 20.0, 21.5, 21.75, 23.25, 25.125, 27.5]
		random.shuffle(data)
		assert statistics.mean(data) == 22.015625
		assert Math.mean(data) == 22.015625
		assert statistics.mean(data) == Math.mean(data)

	def test_decimals(self):
		# Test mean with Decimals.
		data = [Decimal("1.634"), Decimal("2.517"), Decimal("3.912"), Decimal("4.072"), Decimal("5.813")]
		random.shuffle(data)
		assert statistics.mean(data) == Decimal("3.5896")
		assert Math.mean(data) == Decimal("3.5896")
		assert statistics.mean(data) == Math.mean(data)

	def test_fractions(self):
		# Test mean with Fractions.
		data = [
				Fraction(1, 2),
				Fraction(2, 3),
				Fraction(3, 4),
				Fraction(4, 5),
				Fraction(5, 6),
				Fraction(6, 7),
				Fraction(7, 8),
				]
		random.shuffle(data)
		assert statistics.mean(data) == Fraction(1479, 1960)
		assert Math.mean(data) == Fraction(1479, 1960)
		assert statistics.mean(data) == Math.mean(data)

	def test_inf(self):
		# Test mean with infinities.
		raw = [1, 3, 5, 7, 9]  # Use only ints, to avoid TypeError later.
		for kind in (float, Decimal):
			for sign in (1, -1):
				inf = kind("inf") * sign
				data = raw + [inf]
				assert math.isinf(statistics.mean(data))
				assert statistics.mean(data) == inf
				assert math.isinf(Math.mean(data))
				assert Math.mean(data) == inf
				assert statistics.mean(data) == Math.mean(data)

	def test_mismatched_infs(self):
		# Test mean with infinities of opposite sign.
		data = [2, 4, 6, float("inf"), 1, 3, 5, float("-inf")]
		assert math.isnan(statistics.mean(data))
		assert math.isnan(Math.mean(data))

	def test_nan(self):
		# Test mean with NANs.
		raw = [1, 3, 5, 7, 9]  # Use only ints, to avoid TypeError later.
		for kind in (float, Decimal):
			inf = kind("nan")
			data = raw + [inf]
			assert math.isnan(statistics.mean(data))
			assert math.isnan(Math.mean(data))

	def test_big_data(self):
		# Test adding a large constant to every data point.
		c = 1e9
		data = [3.4, 4.5, 4.9, 6.7, 6.8, 7.2, 8.0, 8.1, 9.4]
		expected = statistics.mean(data) + c
		assert expected != c
		assert statistics.mean([x + c for x in data]) == expected
		assert Math.mean([x + c for x in data]) == expected
		assert statistics.mean([x + c for x in data]) == Math.mean([x + c for x in data])

	def test_doubled_data(self):
		# Mean of [a,b,c...z] should be same as for [a,a,b,b,c,c...z,z].
		data = [random.uniform(-3, 5) for _ in range(1000)]
		expected = statistics.mean(data)
		assert statistics.mean(data * 2) == expected
		assert Math.mean(data * 2) == expected
		assert statistics.mean(data * 2) == Math.mean(data * 2)

	def test_regression_20561(self):
		# Regression test for issue 20561.
		# See http://bugs.python.org/issue20561
		d = Decimal("1e4")
		assert statistics.mean([d]) == d
		assert Math.mean([d]) == d
		assert Math.mean([d]) == statistics.mean([d])

	def test_regression_25177(self):
		# Regression test for issue 25177.
		# Ensure very big and very small floats don't overflow.
		# See http://bugs.python.org/issue25177.

		data = [8.988465674311579e307, 8.98846567431158e307]
		assert statistics.mean(data) == 8.98846567431158e307
		assert Math.mean(data) == 8.98846567431158e307
		assert statistics.mean(data) == Math.mean([8.988465674311579e307, 8.98846567431158e307])

		big = 8.98846567431158e307
		tiny = 5e-324
		for n in (2, 3, 5, 200):
			assert statistics.mean([big] * n) == big
			assert statistics.mean([tiny] * n) == tiny
			assert Math.mean([big] * n) == big
			assert Math.mean([tiny] * n) == tiny
			assert statistics.mean([big] * n) == Math.mean([big] * n)
			assert statistics.mean([tiny] * n) == Math.mean([tiny] * n)


class TestMedian:

	def test_even_ints(self):
		# Test median with an even number of int data points.
		data = [1, 2, 3, 4, 5, 6]
		assert len(data) % 2 == 0
		assert statistics.median(data) == 3.5
		assert Math.median(data) == 3.5
		assert statistics.median(data) == Math.median(data)

	def test_odd_ints(self):
		# Test median with an odd number of int data points.
		data = [1, 2, 3, 4, 5, 6, 9]
		assert len(data) % 2 == 1
		assert statistics.median(data) == 4
		assert Math.median(data) == 4
		assert statistics.median(data) == Math.median(data)

	def test_odd_fractions(self):
		# Test median works with an odd number of Fractions.
		data = [Fraction(1, 7), Fraction(2, 7), Fraction(3, 7), Fraction(4, 7), Fraction(5, 7)]
		assert len(data) % 2 == 1
		random.shuffle(data)
		assert statistics.median(data) == Fraction(3, 7)
		assert Math.median(data) == Fraction(3, 7)
		assert statistics.median(data) == Math.median(data)

	def test_even_fractions(self):
		# Test median works with an even number of Fractions.
		data = [
				Fraction(1, 7),
				Fraction(2, 7),
				Fraction(3, 7),
				Fraction(4, 7),
				Fraction(5, 7),
				Fraction(6, 7),
				]
		assert len(data) % 2 == 0
		random.shuffle(data)
		assert statistics.median(data) == Fraction(1, 2)
		assert Math.median(data) == Fraction(1, 2)
		assert statistics.median(data) == Math.median(data)

	def test_odd_decimals(self):
		# Test median works with an odd number of Decimals.
		data = [Decimal("2.5"), Decimal("3.1"), Decimal("4.2"), Decimal("5.7"), Decimal("5.8")]
		assert len(data) % 2 == 1
		random.shuffle(data)
		assert statistics.median(data) == Decimal("4.2")
		assert Math.median(data) == Decimal("4.2")
		assert statistics.median(data) == Math.median(data)

	def test_even_decimals(self):
		# Test median works with an even number of Decimals.
		data = [
				Decimal("1.2"),
				Decimal("2.5"),
				Decimal("3.1"),
				Decimal("4.2"),
				Decimal("5.7"),
				Decimal("5.8"),
				]
		assert len(data) % 2 == 0
		random.shuffle(data)
		assert statistics.median(data) == Decimal("3.65")
		assert Math.median(data) == Decimal("3.65")
		assert statistics.median(data) == Math.median(data)


class TestStdev:
	# Tests for sample variance.
	def setUp(self):
		self.func = statistics.variance

	def test_single_value(self):
		for x in (81, 203.74, 3.9e14, Fraction(5, 21), Decimal("35.719")):
			with pytest.raises(statistics.StatisticsError):
				statistics.stdev([x])  # type: ignore

	def test_ints(self):
		# Test sample variance with int data.
		data: List[int] = [4, 7, 13, 16]
		exact = math.sqrt(30)
		assert statistics.stdev(data) == exact
		assert Math.std(data) == exact
		assert statistics.stdev(data) == Math.std(data)

	def test_fractions(self):
		# Test sample variance with Fraction data.
		data = [Fraction(1, 4), Fraction(1, 4), Fraction(3, 4), Fraction(7, 4)]
		expected = 0.7071067811865476
		assert statistics.stdev(data) == expected
		assert Math.std(data) == expected
		assert statistics.stdev(data) == Math.std(data)

	def test_decimals(self):
		# Test sample variance with Decimal data.
		data = [Decimal(2), Decimal(2), Decimal(7), Decimal(9)]
		exact = (4 * Decimal("9.5") / Decimal(3)).sqrt()
		assert statistics.stdev(data) == exact
		assert Math.std(data) == exact
		assert statistics.stdev(data) == Math.std(data)
		assert isinstance(statistics.stdev(data), Decimal)

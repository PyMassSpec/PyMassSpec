# TODO:
"""proc.py
"""

# this package
from pyms.GCMS.Function import diff
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.IO.JCAMP import JCAMP_reader

# read the raw data
andi_file = "data/gc01_0812_066.cdf"
jcamp_file = "data/gc01_0812_066.jdx"

data1 = ANDI_reader(andi_file)
data2 = JCAMP_reader(jcamp_file)

# trim data2 between scans 1000 and 2000
# data2.trim(begin=1000,end=2000)

diff(data1, data2)

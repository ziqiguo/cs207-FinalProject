from chemkin import *
from parser import *
import numpy as np

# Input variables
x = [2.0, 1.0, 0.5, 1.0, 1.0, 0.0, 0.0, 0.25]
T = 1500
fname = 'rxnset_long.xml'

testcase1 = chemkin()
testcase1.parse(fname)
print(testcase1.reaction_rates(x, T))

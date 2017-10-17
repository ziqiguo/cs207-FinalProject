from chemkin import *
from parser import *
import numpy as np

# Input variables
x = [2., 1., .5, 1., 1., 1.]
T = 1500
fname = 'rxns.xml'

testcase1 = chemkin()
testcase1.parse(fname)
print(testcase1.reaction_rates(x, T))

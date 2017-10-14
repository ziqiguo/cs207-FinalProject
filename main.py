
from chemkin import *
from parser import *

# Input variables
x = [1,1,1,1,1,1]
T = 1500
v1, v2, k = parseXML('rxns.xml', T)
print(chemkin.reaction(k, x, v1, v2))
from chemkin import *
from parser import *
import numpy as np

# Input variables
x = [1,1,1,1,1,1]
T = 1500
reactions_dict = parseXML('rxns.xml')
c = chemkin()

# Calculate reaction rates
k = c.reaction_rates(reactions_dict['rates'], T)

# Get stoichiometric coefficient
v1 = []
v2 = []
for key, value in reactions_dict['reactants'].items():
	v1.append(value)
for key, value in reactions_dict['products'].items():
	v2.append(value)

v1 = np.array(v1).T
v2 = np.array(v2).T

# Calculate progress rates
print(c.reaction(k, x, v1, v2))

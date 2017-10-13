import xml.etree.ElementTree as ET
tree = ET.parse('rxns.xml')
root = tree.getroot()

Arrhenius = []
import numpy as np
for child in root.find('reactionData').findall('reaction'):
    
    # Parse reactants and products
    species_dict = {'reactants':{}, 'products':{}}
    reactants = child.find('reactants').text.split()
    products = child.find('products').text.split()
    for pair in reactants:
        key, value = pair.split(':')
        species_dict['reactants'][key] = float(value)
    for pair in products:
        key, value = pair.split(':')
        species_dict['products'][key] = float(value)
    
    # Parse reaction coefficients
    coeffs = child.find('rateCoeff')
    coeffs_dict = {}
    if coeffs.find('Arrhenius'):
        coeffs_dict['A'] = float(coeffs.find('Arrhenius').find('A').text)
        coeffs_dict['E'] = float(coeffs.find('Arrhenius').find('E').text)
        
    if coeffs.find('modifiedArrhenius'):
        coeffs_dict['A'] = float(coeffs.find('modifiedArrhenius').find('A').text)
        coeffs_dict['E'] = float(coeffs.find('modifiedArrhenius').find('E').text)
        coeffs_dict['b'] = float(coeffs.find('modifiedArrhenius').find('b').text)
    if coeffs.find('Constant'):
        coeffs_dict['k'] = float(coeffs.find('Constant').find('k').text)
        

    print(species_dict)
    print(coeffs_dict)
        
        
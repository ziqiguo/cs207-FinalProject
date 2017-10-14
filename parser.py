import xml.etree.ElementTree as ET
from chemkin import *
import numpy as np

def parseXML(file):
    """
    Parses the XML file containing reactions information. Returns reaction coefficients and reaction rates.
    
    INPUTS
    ========
    file: string, required
          path of the XML file
          
    RETURNS
    ========
    v1: nested list of floats with shape m*n, where m is number of species and n is number of reactions in the system
        coefficients of reactants following the species array's order
    v2: nested list of floats with shape m*n, where m is number of species and n is number of reactions in the system
        coefficients of reactants following the species array's order
    reaction_rate: list of floats with shape n, where n is the number of reactions in the system
    
    EXAMPLES
    ========
    >> parseXML('rxns.xml')
    ([[1.0, 0.0, 0.0],
     [0.0, 1.0, 0.0],
     [0.0, 0.0, 1.0],
     [0.0, 1.0, 1.0],
     [0.0, 0.0, 0.0],
     [1.0, 0.0, 0.0]],
     [[0.0, 1.0, 1.0],
     [1.0, 0.0, 0.0],
     [1.0, 1.0, 0.0],
     [0.0, 0.0, 0.0],
     [0.0, 0.0, 1.0],
     [0.0, 0.0, 0.0]],
     [114837571.22536749, 2310555.9199959813, 1000.0])
    """
    
    tree = ET.parse(file)
    root = tree.getroot()

    # initialize variables
    T = 1500
    species_dict = {'reactants':{}, 'products':{}}
    reaction_rate = []
    reactant_coeffs = []
    product_coeffs = []
    for species in root.find('phase').find('speciesArray').text.split():
        species_dict['reactants'][species] = []
        species_dict['products'][species] = []
    v1 = []
    v2 = []

    for child in root.find('reactionData').findall('reaction'):

        # Parse reactants and products
        species = list(species_dict['reactants'].keys())
        for reactant in child.find('reactants').text.split():
            key, value = reactant.split(':')
            species_dict['reactants'][key].append(float(value))
            species.remove(key)
        for reactant in species:
            species_dict['reactants'][reactant].append(0.0)

        species = list(species_dict['products'].keys())
        for product in child.find('products').text.split():
            key, value = product.split(':')
            species_dict['products'][key].append(float(value))
            species.remove(key)
        for product in species:
            species_dict['products'][product].append(0.0)

        # Parse reaction coefficients
        coeffs = child.find('rateCoeff')
        if coeffs.find('Arrhenius'):
            coeffs_dict['A'] = float(coeffs.find('Arrhenius').find('A').text)
            coeffs_dict['E'] = float(coeffs.find('Arrhenius').find('E').text)
            reaction_rate.append(chemkin.arrhenius(coeffs_dict['A'], coeffs_dict['E'], T))
        if coeffs.find('modifiedArrhenius'):
            coeffs_dict['A'] = float(coeffs.find('modifiedArrhenius').find('A').text)
            coeffs_dict['E'] = float(coeffs.find('modifiedArrhenius').find('E').text)
            coeffs_dict['b'] = float(coeffs.find('modifiedArrhenius').find('b').text)
            reaction_rate.append(chemkin.modified(coeffs_dict['A'], coeffs_dict['b'], coeffs_dict['E'], T))
        if coeffs.find('Constant'):
            coeffs_dict['k'] = float(coeffs.find('Constant').find('k').text)
            reaction_rate.append(chemkin.constant(coeffs_dict['k']))

    for key, value in species_dict['reactants'].items():
        v1.append(value)
    for key, value in species_dict['products'].items():
        v2.append(value)
        
    return v1, v2, reaction_rate
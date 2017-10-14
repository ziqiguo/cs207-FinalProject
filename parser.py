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
    reactions_dict = {'reactants': {}, 'products': {}, 'rates': []}
    coeffs_dict = []

    for species in root.find('phase').find('speciesArray').text.split():
        reactions_dict['reactants'][species] = []
        reactions_dict['products'][species] = []

    for child in root.find('reactionData').findall('reaction'):

        # Parse reactants and products
        species = list(reactions_dict['products'].keys())
        for reactant in child.find('reactants').text.split():
            key, value = reactant.split(':')
            reactions_dict['reactants'][key].append(float(value))
            species.remove(key)
        for reactant in species:
            reactions_dict['reactants'][reactant].append(0.0)

        species = list(reactions_dict['products'].keys())
        for product in child.find('products').text.split():
            key, value = product.split(':')
            reactions_dict['products'][key].append(float(value))
            species.remove(key)
        for product in species:
            reactions_dict['products'][product].append(0.0)

        # Parse reaction coefficients
        coeffs = child.find('rateCoeff')
        if coeffs.find('Arrhenius'):
            new_dict = {}
            new_dict['type'] = 'Arrhenius'
            new_dict['A'] = float(coeffs.find('Arrhenius').find('A').text)
            new_dict['E'] = float(coeffs.find('Arrhenius').find('E').text)
            reactions_dict['rates'].append(new_dict)
            
        if coeffs.find('modifiedArrhenius'):
            new_dict = {}
            new_dict['type'] = 'modifiedArrhenius'
            new_dict['A'] = float(coeffs.find('modifiedArrhenius').find('A').text)
            new_dict['E'] = float(coeffs.find('modifiedArrhenius').find('E').text)
            new_dict['b'] = float(coeffs.find('modifiedArrhenius').find('b').text)
            reactions_dict['rates'].append(new_dict)

        if coeffs.find('Constant'):
            new_dict = {}
            new_dict['type'] = 'Constant'
            new_dict['k'] = float(coeffs.find('Constant').find('k').text)
            reactions_dict['rates'].append(new_dict)

    return reactions_dict
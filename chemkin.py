import numpy as np
from parser import *

class chemkin:
    def __init__(self):
        self.r = 8.314
        self.v1 = None
        self.v2 = None
        self.rates = None
        self.k = None

    def parse(self, file):
        """Calls parser and stores instance attributes

        INPUTS
        =======
        file: string, path of XML file to be parsed

        RETURNS
        =======
        None. Instance attributes added:
        rates: dictionary of reaction rate coefficients information
               need to call reaction_rate() to calculate list k
        v1: coefficients of reactants
        v2: coefficients of products

        EXAMPLES
        ========
        >>> c = chemkin()
        >>> c.parse('rxns.xml')
        >>> c.v1
        array([[ 1.,  0.,  0.,  0.,  0.,  1.],
               [ 0.,  1.,  0.,  1.,  0.,  0.],
               [ 0.,  0.,  1.,  1.,  0.,  0.]])
        """
        reactions_dict = parseXML('rxns.xml')
        self.rates = reactions_dict['rates']

        self.v1 = []
        self.v2 = []
        for key, value in reactions_dict['reactants'].items():
            self.v1.append(value)
        for key, value in reactions_dict['products'].items():
            self.v2.append(value)

        self.v1 = np.array(self.v1).T
        self.v2 = np.array(self.v2).T



    def k_constant(self, k):
        """Returns the constant reaction rate coefficient.

        INPUTS
        =======
        k: int or float, the constant reaction rate coefficient

        RETURNS
        ========
        k: float, the constant reaction rate coefficient

        EXAMPLES
        =========
        >>> chemkin().k_constant(10.0)
        10.0
        """
        try:
            k = float(k)
        except:
            raise TypeError("Error: unable to convert k to float!")
        return k

    def k_arrhenius(self, a, e, t):
        """Returns the Arrhenius reaction rate coefficient.

        INPUTS
        =======
        a: int or float, Arrhenius prefactor, A, is strictly positive
        e: int or float, the activation energy for the reaction
        t: int or float, the temperature T, must be positive (assuming a Kelvin scale)

        RETURNS
        ========
        k: float, the Arrhenius reaction rate coefficient

        EXAMPLES
        =========
        >>> chemkin().k_arrhenius(10,10,10)
        8.8667297841210573
        """
        try:
            a = float(a)
            e = float(e)
            t = float(t)
        except:
            raise ValueError("Error: unable to convert all parameters to float!")
        if a <= 0:
            raise ValueError("The Arrhenius prefactor A must be positive")
        if t <= 0:
            raise ValueError("The temperature T must be positive (assume a Kelvin scale)")
        return a*np.exp(-e/(self.r*t))

    def k_modified(self, a, b, e, t):
        """Returns the modified Arrhenius reaction rate coefficient.

        INPUTS
        =======
        a: int or float, Arrhenius prefactor, A, is strictly positive
        b: int or float, fitted rate constant
        e: int or float, the activation energy for the reaction
        t: int or float, the temperature T, must be positive (assuming a Kelvin scale)

        RETURNS
        ========
        k: float, the modified Arrhenius reaction rate coefficient

        EXAMPLES
        =========
        >>> chemkin().k_modified(10**7,0.5,10**3,10**2)
        30035490.889639609
        """
        if isinstance(b, complex):
            raise ValueError("The modified Arrhenius parameter b must be real")
        try:
            a = float(a)
            b = float(b)
            e = float(e)
            t = float(t)
        except:
            raise ValueError("Error: unable to convert all parameters to float!")
        if a <= 0:
            raise ValueError("The Arrhenius prefactor A must be positive")
        if t <= 0:
            raise ValueError("The temperature T must be positive (assume a Kelvin scale)")
        return a*(t**b)*np.exp(-e/(self.r*t))

    def k_system(self, T):
        """Calculates a list of k, one for each reaction.

        INPUTS
        =======
        rates: dictionary containing reaction rate type and all parameters pertaining to set of reactions
        T: int or float, environment temperature

        RETURNS
        =======
        None. Class attributes are added:
        k: list of floats, has length m where m is the number of reactions

        EXAMPLES
        >>> c = chemkin()
        >>> c.parse('rxns.xml')
        >>> c.k_system(1500)
        >>> c.k
        [114837571.22536749, 2310555.9199959813, 1000.0]
        """
        if self.rates is None:
            raise ValueError("Rates not initialized. Please call Chemkin().parse() to parse XML first.")

        self.k = []
        for reaction in self.rates:
            if reaction['type'] == 'Arrhenius':
                self.k.append(self.k_arrhenius(reaction['A'], reaction['E'], T))
            elif reaction['type'] == 'modifiedArrhenius':
                self.k.append(self.k_modified(reaction['A'], reaction['b'], reaction['E'], T))
            elif reaction['type'] == 'Constant':
                self.k.append(self.k_constant(reaction['k']))

    def progress_reaction(self, k, x, v):
        """Returns the progress rate of a single reaction.

        INPUTS
        =======
        k: int or float, the reaction rate coefficient
        x: list of concentration of each species
        v: list of Stoichiometric coefficients of reactants

        RETURNS
        ========
        progress rate: float, the progress rate of a reaction

        EXAMPLES
        =========
        >>> chemkin().progress_reaction(10,[1,2,3],[2,1,0])
        20
        """
        if len(x) != len(v):
            raise ValueError("Dimensions of concentration and coefficients do not match!")
        if k==0:
            raise ValueError("Reaction rate of reaction is 0. No reaction.")
        if x==[0]*len(x):
            raise ValueError("All reactant concentrations for reaction are 0.")
        if list(v)==[0]*len(v):
            raise ValueError("All stoich coefficients for reaction are 0.")

        progress = 1
        for i in range(len(x)):
            progress = progress * (x[i]**v[i])
        return k*progress

    def progress_system(self, x):
        """Returns the progress rate of a system of reactions.

        INPUTS
        =======
        x: list of concentration of each species

        RETURNS
        ========
        progress rate: list of progress rate of each reaction
        
        EXAMPLES
        ========
        >>> c = chemkin()
        >>> c.parse('rxns.xml')
        >>> c.k_system(1500)
        >>> c.progress_system([2., 1., .5, 1., 1., 1.])
        [229675142.45073497, 2310555.9199959813, 500.0]
        """
        if self.k is None:
            raise ValueError("Reaction rate coefficients not initialized. Please call Chemkin().k_system() to calculate list k first.")

        if len(self.v1) != len(self.k):
            raise ValueError("Number of k does not much number of reactions!")

        for lst in self.v1:
            if any(isinstance(d, complex) for d in lst) == True:
                raise ValueError('Complex value in reactant coefficient detected!')

        progress = []
        for i in range(len(self.v1)):
            if len(self.v1[i]) != len(x):
                raise ValueError("Not enough coefficient values! Check the dimension of coefficient matrix.")
            progress.append(self.progress_reaction(self.k[i], x, self.v1[i]))
        return progress

    def reaction_rates(self, x, T):
        """Returns the reaction rate of a system of reactions for each specie.

        INPUTS
        =======
        x: float, list of n, where n is the number of species
           concentrations of each species
        T: float, required
           environment temperature, used to calculate reaction rate coefficients

        RETURNS
        ========
        reaction rate: list of rate of consumption or formation of specie
        
        EXAMPLES
        ========
        >>> c = chemkin()
        >>> c.parse('rxns.xml')
        >>> c.reaction_rates([2., 1., .5, 1., 1., 1.], 1500)
        [-227364086.53073898, 227364586.53073898, 231985198.37073097, -2311055.9199959813, 500.0, -229675142.45073497]
        """
        self.k_system(T)

        if np.array(self.v1).shape != np.array(self.v2).shape:
            raise ValueError("Dimensions of coefficients of reactants and products do not match.")
        for lst in self.v2:
            if any(isinstance(d, complex) for d in lst) == True:
                raise ValueError('Complex value in product coefficient detected!')

        progress = self.progress_system(x)
        reaction_rates = []

        for i in range(len(x)):
            reaction_rates.append(0)
            for j in range(len(progress)):
                reaction_rates[i] = reaction_rates[i] + ((self.v2[j][i]-self.v1[j][i])*progress[j])
        return reaction_rates

    def __str__(self):
        return str(vars(self))

    def __repr__(self):
        class_name = type(self).__name__
        return class_name + '()'

import numpy as np

class chemkin:
    def __init__(self):
        self.r = 8.314

    def constant(self, k):
        """Returns the constant reaction rate coefficient.

        INPUTS
        =======
        k: int or float, the constant reaction rate coefficient

        RETURNS
        ========
        k: float, the constant reaction rate coefficient

        EXAMPLES
        =========
        >>> constant(10.0)
        10.0
        """
        try:
            k = float(k)
        except ValueError:
            return "Error: unable to convert k to float!"
        return k

    def arrhenius(self, a, e, t):
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
        >>> arrhenius(10,10,10)
        8.8667297841210573
        """
        try:
            a = float(a)
            e = float(e)
            t = float(t)
        except TypeError:
            return "Error: unable to convert all parameters to float!"
        if a <= 0:
            raise ValueError("Arrhenius prefactor a is non-positive!")
        if t <= 0:
            raise ValueError("Temperature t is non-positive!")
        return a*np.exp(-e/(self.r*t))

    def modified(self, a, b, e, t):
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
        >>> modified(10**7,0.5,10**3,10**2)
        30035490.889639609
        """
        try:
            a = float(a)
            b = float(b)
            e = float(e)
            t = float(t)
        except TypeError:
            return "Error: unable to convert all parameters to float!"
        if a <= 0:
            raise ValueError("Arrhenius prefactor a is non-positive!")
        if isinstance(b, complex):
            raise ValueError("b is complex number!")
        if t <= 0:
            raise ValueError("Temperature t is non-positive!")
        return a*(t**b)*np.exp(-e/(self.r*t))

    def reaction_rates(self, rates, T):
        """Calculates reaction rates for a list of reactions.
        
        INPUTS
        =======
        rates: dictionary containing reaction rate type and all needed parameters
        T: int or float, environment temperature

        RETURNS
        =======
        k: list of floats, has length m where m is the number of reactions

        EXAMPLES
        >> reaction_rates([{'type': 'Arrhenius', 'A': 35200000000.0, 'E': 71400.0}, {'type': 'modifiedArrhenius', 'A': 0.0506, 'E': 26300.0, 'b': 2.7}, {'type': 'Constant', 'k': 1000.0}], 1500)
        [114837571.22536749, 2310555.9199959813, 1000.0]
        """

        k = []
        for reaction in rates:
            if reaction['type'] == 'Arrhenius':
                k.append(self.arrhenius(reaction['A'], reaction['E'], T))
            elif reaction['type'] == 'modifiedArrhenius':
                k.append(self.modified(reaction['A'], reaction['b'], reaction['E'], T))
            elif reaction['type'] == 'Constant':
                k.append(self.constant(reaction['k']))
        return k

    def progress_u(self, k, x, v):
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
        >>> progress(10,[1,2,3],[2,1,0])
        20
        """
        if len(x) != len(v):
            raise ValueError("Length of x and v does not match!")
        n = len(x)
        p = 1
        for i in range(n):
            p = p * (x[i]**v[i])
        return k*p

    def progress(self, k, x, v1):
        """Returns the progress rate of a system of reactions.

        INPUTS
        =======
        k: list of reaction rate coefficients for each reaction
        x: list of concentration of each species
        v1: matrix (list of list) of Stoichiometric coefficients of reactants of each reaction

        RETURNS
        ========
        progress rate: list of progress rate of each reaction

        EXAMPLES
        =========
        >>> progress_m([10,10],[1,2,1],[[1,2,0],[2,0,2]])
        [40, 10]
        """
        m = len(v1)
        if m != len(k):
            raise ValueError("Number of k does not much number of reactions!")
        n = len(x)
        w = []
        for i in range(m):
            if len(v1[i]) != n:
                raise ValueError("Error in dimension of v values!")
            w.append(self.progress_u(k[i],x,v1[i]))
        return w

    def reaction(self, k, x, v1, v2):
        """Returns the reaction rate of a system of reactions for each specie.

        INPUTS
        =======
        k: list of reaction rate coefficients for each reaction
        x: list of concentration of each species
        v1: matrix (list of list) of Stoichiometric coefficients of reactants of each reaction
        v2: matrix (list of list) of Stoichiometric coefficients of products of each reaction

        RETURNS
        ========
        reaction rate: list of rate of consumption or formation of specie

        EXAMPLES
        =========
        >>> reaction([10,10],[1,2,1],[[1,2,0],[0,0,2]],[[0,0,1],[1,2,0]])
        [-30, -60, 20]
        """
        w = self.progress(k,x,v1)
        f = []
        m = len(w)
        if m != len(k):
            raise ValueError("Number of k does not much number of reactions!")
        n = len(x)
        for i in range(n):
            f.append(0)
            for j in range(m):
                if len(v1[j]) != n or len(v2[j]) != n:
                    raise ValueError("Error in dimension of v values!")
                f[i] = f[i] + ((v2[j][i]-v1[j][i])*w[j])
        return f


    def __repr__(self):
        class_name = type(self).__name__
        return "This is {0}".format(class_name)

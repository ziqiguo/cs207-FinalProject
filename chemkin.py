def progress(k,x,v):
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

def progress_m(k,x,v1):
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
        w.append(progress(k[i],x,v1[i]))
    return w

def reaction(k,x,v1,v2):
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
    w = progress_m(k,x,v1)
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
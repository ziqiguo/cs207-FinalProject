import numpy as np

def constant(k):
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

def arrhenius(a,e,t):
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
    return a*np.exp(-e/(8.314*t))
    
def modified(a,b,e,t):
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
    return a*(t**b)*np.exp(-e/(8.314*t))
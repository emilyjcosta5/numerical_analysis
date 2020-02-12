from scipy import misc
from math import exp

def fixed_point(p_0, f, E, N_0):
    '''
    Finds a fixed point in a function, which is the value
    of the function that does not change when the function
    is applied.

    Parameters
    ----------
    p_0: float
        Initial approximation.
    f: function
        Function to find fixed-points in.
    E: float
        Tolerance.
    N_0: int
        Maximum number of iterations.

    Returns
    -------
    p: float
        Approximate solution.
    '''
    p = 0
    for n in N_0:
        p = f(p_0)
        if(abs(p-p_0)<E):
            return p
        p_0 = p
    print('No solution found.')
    return None

def newtons_method(f, p_0, E, N):
    p = p_0
    n = 1
    for n in range(0,N):
        p = p_0 - f(p_0)/misc.derivative(f, p_0)
        if abs(p-p_0)<E:
            break
        p_0 = p
    return p_0

def f(x):
    return x**3 + x**2 + 2*x 

def g(x):
    return exp(x) - x - 1

def g_1(x);
    return (1/2) * (x**3 + 1)

def g_2(x):
    return (2/x) - (1/(x^2))

if __name__=='__main__':
    # Standardize parameters
    TOL = 0.01
    p_0 = 0.5
    fp = 10
    E = 0.0001
    N = 1000

    # 1. a.
    print(fixed_point(p_0, g_1, TOL, N_0))
    print(fixed_point(p_0, g_2, TOL, N_0))

    # 3. 1. a. 
    print(newtons_method(f, fp, E, N))
    # b.
    print(newtons_method(g, fp, E, N))

    # 2. a.

    # b.

    # 3. a. 

    # b.

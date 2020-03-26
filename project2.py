# Author: Emily Costa
# Created on: March 26, 2020
# Various composite methods for computing the integral of a function.

import numpy as np

def trapezoidal(f, a, b, N, I):
    '''
    Implements the trapezoidal rule to approx. the integral of f(x),
    with the N subintervals from a to b.

    Parameters
    ----------
    f : function
        Single variable function to be integrated.
    a , b : numbers
        Interval of integration [a,b]
    N : integer
        Number of subintervals of [a,b]
    I : float
        Exact value of the solution.

    Returns
    -------
    approx : float
        Final approx. of the integral of f(x).
    err : float
        Absolute error of the approx.

    '''
    x = np.linspace(a,b,N+1)
    y = f(x)
    y_right = y[1:] 
    y_left = y[:-1]
    dx = (b - a)/N
    approx = (dx/2) * np.sum(y_right + y_left)
    err = np.abs(I - approx)
    return approx, err

def simpsons(f, a, b, N, I):
    '''
    Implements Simpson's rule to approx. the integral of f(x),
    with the N subintervals from a to b.

    Parameters
    ----------
    f : function
        Single variable function to be integrated.
    a , b : numbers
        Interval of integration [a,b]
    N : integer
        Number of subintervals of [a,b]
    I : float
        Exact value of the solution.

    Returns
    -------
    approx : float
        Final approx. of the integral of f(x).
    err : float
        Absolute error of the approx.
    '''
    dx = (b-a)/N
    x = np.linspace(a,b,N+1)
    y = f(x)
    approx = dx/3 * np.sum(y[0:-1:2] + 4*y[1::2] + y[2::2])
    err = np.abs(I - approx)
    return approx, err

def midpoint(f, a, b, N, I):
    '''
    Implements Simpson's rule to approx. the integral of f(x),
    with the N subintervals from a to b.

    Parameters
    ----------
    f : function
        Single variable function to be integrated.
    a , b : numbers
        Interval of integration [a,b]
    N : integer
        Number of subintervals of [a,b]
    I : float
        Exact value of the solution.

    Returns
    -------
    approx : float
        Final approx. of the integral of f(x).
    err : float
        Absolute error of the approx.
    '''

    dx = np.float(b-a)/N
    x = np.linspace(a+dx/2, b-dx/2, N)
    approx = dx*np.sum(f(x))
    err = np.abs(I - approx)
    return approx, err

if __name__ == "__main__":
    # Exact value for calculating abs. err.
    I = np.exp(4)*(2*np.sin(6)-3*np.cos(6))/13 + 3/13
    a = 0
    b = 2
    f0 = lambda x : np.exp(2*x)*np.sin(3*x)

    print("(Approx., Abs. err.)")

    # Test composite trapezoidal.
    print("Composite trapezoidal.")
    print("N=4: ", trapezoidal(f0, a, b, 4, I))
    print("N=8: ", trapezoidal(f0, a, b, 8, I))
    print("N=16: ", trapezoidal(f0, a, b, 16, I))

    # Test composite Simpson.
    print("Composite Simpsons's.")
    print("N=4: ", simpsons(f0, a, b, 4, I))
    print("N=8: ", simpsons(f0, a, b, 8, I))
    print("N=16: ", simpsons(f0, a, b, 16, I))

    # Test composite midpoint.
    print("Composite midpoint.")
    print("N=4: ", midpoint(f0, a, b, 4, I))
    print("N=8: ", midpoint(f0, a, b, 8, I))
    print("N=16: ", midpoint(f0, a, b, 16, I))

# Output.
'''
(Approx., Abs. err.)
Composite trapezoidal.
N=4:  (-11.753892882857468, 2.460084247005053)
N=8:  (-13.57597939179939, 0.637997738063131)
N=16:  (-14.053231044872895, 0.16074608498962561)
Composite Simpsons's.
N=4:  (-13.476843513846811, 0.7371336160157096)
N=8:  (-14.183341561446696, 0.03063556841582482)
N=16:  (-14.212314929230729, 0.001662200631791677)
Composite midpoint.
N=4:  (-15.398065900741308, 1.1840887708787875)
N=8:  (-14.530482697946404, 0.3165055680838833)
N=16:  (-14.294200161426023, 0.08022303156350219)
'''

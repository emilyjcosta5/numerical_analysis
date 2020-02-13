# Written by: Emily Costa

from scipy import misc
from math import exp
import pandas as pd

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
        Epsilon, threshold for error tolerance analysis.
    N_0: int
        Maximum number of iterations.

    Returns
    -------
    p: float
        Approximate solution.
    '''
    p = 0
    for n in range(0, N_0):
        print("Iteration ", n, " is ", p_0)
        p = f(p_0)
        if(abs(p-p_0)<E):
            return p
        p_0 = p
    print('No solution found.')
    return None

def newtons_method(f, p_0, E, N):
    '''
    Finds a zero of a function determined by Newton's 
    Method.

    Parameters
    ----------
    f: function
        Function to find zeros in.
    p_0: float
        Initial approximation to zero.
    E: float
        Epsilon, threshold for error tolerance analysis.
    N_0: int
        Maximum number of iterations.

    Returns
    -------
    p_0: float
        A zero of the function determined by method.
    '''
    p = p_0
    for n in range(0,N):
        print("Iteration ", n, " is ", p)
        p = p_0 - f(p_0)/misc.derivative(f, p_0)
        if abs(p-p_0)<E:
            break
        p_0 = p
    return p_0

def bisection_method(f, a, b, E, N_0):
    '''
    Finds a zero of a function determined by Newton's 
    Method.

    Parameters
    ----------
    f: function
        Function to find zeros in.
    a, b: floats
        Range for finding the zero of the function.
    E: float
        Epsilon, threshold for error tolerance analysis.
    N_0: int
        Maximum number of iterations.

    Returns
    -------
    p_0: float
        A zero of the function determined by method.
    '''
    fa = f(a)
    p = 0
    for n in range(0, N_0):
        print("Iteration ", n, " is ", p)
        p = a + (b - a) / 2
        fp = f(p)
        if fp == 0 or b-a < E:
            return p
        elif fa*fp > 0:
            a = p 
            fa = fp
        else:
            b = p
    return None

#def secant_method

def f(x):
    return x**3 + x**2 + 2*x 

def g(x):
    return exp(x) - x - 1

def g_1(x):
    return (1/2) * (x**3 + 1)

def g_2(x):
    return (2/x) - (1/(x**2))

if __name__=='__main__':
    # Standardize parameters
    TOL = 0.01
    p_0 = 0.5
    fp = 10
    E = 0.0001
    N = 1000

    # 1. a.
    print("1.a.i.: ")
    print("Solution: ", fixed_point(p_0, g_1, TOL, N))
    print("1.a.ii.: ")
    print("Solution: ", fixed_point(p_0, g_2, TOL, N))


    print("Newton's Method")
    # 3. 1. a. 
    print("a.")
    print("Final solution is ", newtons_method(f, fp, E, N))
    # b.
    print("b.")
    print("Final solution is ", newtons_method(g, fp, E, N))

    print("Bisection Method")
    # 2. a.
    print("a.")
    print("Final solution is ", bisection_method(f, -2, 4, E, N))
    # b.
    print("b.")
    print("Final solution is ", bisection_method(g, -2, 4, E, N))

    # 3. a. 

    # b.

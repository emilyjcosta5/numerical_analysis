'''
Project 1 for Numerical Analysis, MAD 3401, at Florida
International University. This code is posted for the
sake of turning in to the professor-- please, no 
plagarism.

This project includes function that determine roots
or fixed points of a function using various numerical
methods.

Citations
---------
All code is based off the psuedocode in Chapter two of:
Numerical Analysis,
Tenth Edition
Richard L. Burden, J. Douglas Faires,
Annette M. Burden
'''

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
        if p_0 == 0:
            return 0
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
    for n in range(1,N):
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
    for n in range(1, N_0):
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

def secant_method(p_0, p_1, f, E, N_0):
    '''
    Finds a zero of a function determined by the Secant Method.

    Parameters
    ----------
    p_0, p_1: floats
        Initial points on x-axis to estimate the next points.
    f: function
        Function to find zeros in.
    E: float
        Epsilon, threshold for error tolerance analysis.
    N_0: int
        Maximum number of iterations.

    Returns
    -------
    p_0: float
        A zero of the function determined by method.
    '''
    q_0 = f(p_0)
    q_1 = f(p_1)
    p = 0
    for n in range(2, N_0):
        p = p_1 - q_1*(p_1-p_0)/(q_1-q_0)
        print("Iteration ", n, " is ", p)
        if(abs(p-p_1)<E):
            return p
        p_0 = p_1
        q_0 = q_1
        p_1 = p
        q_1 = f(p)
    return None

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
    a = -2
    b = 4

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
    print("Final solution is ", bisection_method(f, a, b, E, N))
    # b.
    print("b.")
    print("Final solution is ", bisection_method(g, a, b, E, N))

    print("Secant Method")
    # 3. a. 
    print("a.")
    print("Final solution is ", secant_method(a, b, f, E, N))
    # b.
    print("b. ")
    print("Final solution is ", secant_method(a, b, g, E, N))

# Output:
'''
(base) emily@emily-XPS-13-9370:~/Documents/numerical_analysis/numerical_analysis$ python3 project.py 
1.a.i.: 
Iteration  0  is  0.5
Iteration  1  is  0.5625
Iteration  2  is  0.5889892578125
Iteration  3  is  0.602162644566306
Solution:  0.6091720424515518
1.a.ii.: 
Iteration  0  is  0.5
Iteration  1  is  0.0
Solution:  0
Newton's Method
a.
Iteration  1  is  10
Iteration  2  is  6.53250773993808
Iteration  3  is  4.210944932149095
Iteration  4  is  2.6506588264700066
Iteration  5  is  1.5971679620149126
Iteration  6  is  0.8880289733277671
Iteration  7  is  0.43087128789945545
Iteration  8  is  0.17573159318934778
Iteration  9  is  0.06314172821274934
Iteration  10  is  0.021550931294736342
Iteration  11  is  0.007237781930462588
Iteration  12  is  0.0024185121053225227
Iteration  13  is  0.0008068242708295856
Iteration  14  is  0.0002690138891456164
Iteration  15  is  8.967934236907125e-05
Final solution is  8.967934236907125e-05
b.
Iteration  1  is  10
Iteration  2  is  9.149473962415037
Iteration  3  is  8.299396796425302
Iteration  4  is  7.450266702958602
Iteration  5  is  6.6031092799224815
Iteration  6  is  5.759991967649747
Iteration  7  is  4.924960789433266
Iteration  8  is  4.1055964099469815
Iteration  9  is  3.315191462300527
Iteration  10  is  2.5747745547278016
Iteration  11  is  1.912640000781003
Iteration  12  is  1.3580647582554017
Iteration  13  is  0.9293214642274574
Iteration  14  is  0.6240055072287617
Iteration  15  is  0.4208992789607361
Iteration  16  is  0.29127521013746094
Iteration  17  is  0.20943864869335233
Iteration  18  is  0.15699600172119574
Iteration  19  is  0.12234050981936521
Iteration  20  is  0.09857523563434947
Iteration  21  is  0.08166258336403984
Iteration  22  is  0.069209623276474
Iteration  23  is  0.059760718396981
Iteration  24  is  0.05240215972486112
Iteration  25  is  0.04654167855524939
Iteration  26  is  0.04178340776085847
Iteration  27  is  0.03785515996143648
Iteration  28  is  0.034564920679800444
Iteration  29  is  0.031774067205105797
Iteration  30  is  0.02938042522716828
Iteration  31  is  0.02730727038406888
Iteration  32  is  0.025496016901909817
Iteration  33  is  0.023901249714862835
Iteration  34  is  0.022487280474183125
Iteration  35  is  0.02122571567455202
Iteration  36  is  0.020093710259541292
Iteration  37  is  0.019072693931136
Iteration  38  is  0.018147428905866087
Iteration  39  is  0.01730530367344223
Iteration  40  is  0.016535797205819415
Iteration  41  is  0.015830067905796795
Iteration  42  is  0.015180634966624045
Iteration  43  is  0.014581128975846892
Iteration  44  is  0.014026094957177488
Iteration  45  is  0.013510835517905297
Iteration  46  is  0.013031284954704146
Iteration  47  is  0.012583907464833795
Iteration  48  is  0.012165614279870765
Iteration  49  is  0.011773695767349848
Iteration  50  is  0.011405765457586235
Iteration  51  is  0.01105971363611031
Iteration  52  is  0.010733668658262082
Iteration  53  is  0.010425964535612326
Iteration  54  is  0.010135113645512143
Iteration  55  is  0.009859783648260602
Iteration  56  is  0.009598777877854967
Iteration  57  is  0.009351018614400327
Iteration  58  is  0.009115532758325898
Iteration  59  is  0.00889143951534729
Iteration  60  is  0.008677939771940072
Iteration  61  is  0.008474306897852342
Iteration  62  is  0.008279878757886203
Iteration  63  is  0.008094050752242886
Iteration  64  is  0.007916269734782959
Iteration  65  is  0.007746028683198061
Iteration  66  is  0.0075828620152512466
Iteration  67  is  0.007426341461883248
Iteration  68  is  0.007276072421729898
Iteration  69  is  0.007131690733035256
Iteration  70  is  0.006992859808439095
Iteration  71  is  0.006859268086111774
Iteration  72  is  0.006730626757341256
Iteration  73  is  0.006606667736354206
Iteration  74  is  0.006487141842858474
Iteration  75  is  0.0063718171718444416
Iteration  76  is  0.006260477628577821
Iteration  77  is  0.006152921609650714
Iteration  78  is  0.00604896081343738
Iteration  79  is  0.0059484191654023055
Final solution is  0.0059484191654023055
Bisection Method
a.
Iteration  1  is  0
Iteration  2  is  1.0
Iteration  3  is  -0.5
Iteration  4  is  0.25
Iteration  5  is  -0.125
Iteration  6  is  0.0625
Iteration  7  is  -0.03125
Iteration  8  is  0.015625
Iteration  9  is  -0.0078125
Iteration  10  is  0.00390625
Iteration  11  is  -0.001953125
Iteration  12  is  0.0009765625
Iteration  13  is  -0.00048828125
Iteration  14  is  0.000244140625
Iteration  15  is  -0.0001220703125
Iteration  16  is  6.103515625e-05
Iteration  17  is  -3.0517578125e-05
Final solution is  1.52587890625e-05
b.
Iteration  1  is  0
Iteration  2  is  1.0
Iteration  3  is  2.5
Iteration  4  is  3.25
Iteration  5  is  3.625
Iteration  6  is  3.8125
Iteration  7  is  3.90625
Iteration  8  is  3.953125
Iteration  9  is  3.9765625
Iteration  10  is  3.98828125
Iteration  11  is  3.994140625
Iteration  12  is  3.9970703125
Iteration  13  is  3.99853515625
Iteration  14  is  3.999267578125
Iteration  15  is  3.9996337890625
Iteration  16  is  3.99981689453125
Iteration  17  is  3.999908447265625
Final solution is  3.9999542236328125
Secant Method
a.
Iteration  2  is  -1.5
Iteration  3  is  -1.2537313432835822
Iteration  4  is  -0.666449093507692
Iteration  5  is  -0.26228572710584236
Iteration  6  is  0.007081901129849244
Iteration  7  is  -0.0007635823097413615
Iteration  8  is  -2.712260791112972e-06
Iteration  9  is  1.0351199691201376e-09
Final solution is  1.0351199691201376e-09
b. 
Iteration  2  is  -2.1405616189355294
Iteration  3  is  -2.300382680926525
Iteration  4  is  -0.7290134575938705
Iteration  5  is  -0.4496802382830732
Iteration  6  is  -0.2523606355721279
Iteration  7  is  -0.15291569782171738
Iteration  8  is  -0.0921926930724927
Iteration  9  is  -0.05641244754494613
Iteration  10  is  -0.0345890199329647
Iteration  11  is  -0.021288662150817428
Iteration  12  is  -0.013120061128228668
Iteration  13  is  -0.008095408761529569
Iteration  14  is  -0.004998003607203054
Iteration  15  is  -0.0030869876795535994
Iteration  16  is  -0.0019071091460149592
Iteration  17  is  -0.0011783730336519032
Iteration  18  is  -0.0007281650821634237
Iteration  19  is  -0.00044998910076748667
Iteration  20  is  -0.00027809261165860546
Iteration  21  is  -0.00017186460383372074
Iteration  22  is  -0.0001062158417794725
Final solution is  -0.0001062158417794725
'''
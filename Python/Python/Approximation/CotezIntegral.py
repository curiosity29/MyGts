from sympy import *
import numpy as np
import math

cotezCoefs = [[1/2,     1/2],
             [1/6,      4/6,        1/6],
             [1/8,      3/8,        3/8,        1/8],
             [7/90,     32/90,      12/90,      32/90,      7/90],
             [19/288,   75/288,     50/288,     50/288,     75/288,19/288]]

#h^n - error term
convergeSpeeds = [2, 4, 4, 6, 6]

multiplier = [12, 90, 80/3, 945/8]

def CotezCoef(n):
    return cotezCoefs[n-1]
# M is sup(|f^(k)|) where f^(k) is the k-derivative of f, n is the degree of polynomial used
# k is the converge speed
# n < 5 pls ...
def Integral(f, M, rangeX, epsilon, n = 2):
    cotezCoef = CotezCoef(n)
    speed = convergeSpeeds[n]
    mult = multiplier[n-1]

    a,b = rangeX
    length = b-a

    #M/multiplier * length * h^k < eps
    h = math.pow(epsilon*mult/M/length, 1/speed)

    steps = (int) (length/h) + 1
    dx = length/steps/n

    integral = 0
    xn = a
    x = symbols('x')
    for step in range(steps):
        for t in range(n+1):
            integral += f.subs(x, xn)*cotezCoef[t]
            xn += dx
        xn -= dx
    
    return integral*dx*n


# test 
#

def test1():
    x = symbols('x')
    f = sympify("x**3 -3*x + 7*x**2 - 5 + 10*sin(100*x) + 2**x")
    print(f)
    integral = Integral(f, 1000, (-3,3), 0.001, 2)
    return integral.evalf()



#print(test1())



def Ln12(t, epsilon):
    if not 1<=t<=2:
        raise ArithmeticError

    x = symbols('x')
    f = 1/x
    return Integral(f, 2**5, (1,t), epsilon, 4)


#ln2 = Ln12(2, 10**-12)

def Ln(x, epsilon):
    n = 0
    while x>2:
        x/=2
        n+=1
    while x<1:
        x*=2
        n+=1

    return n*ln2 + Ln12(x, epsilon)
#print(Ln(123456, 10**-5))

x = symbols('x')
f = 1/(x**2+1)
integral = 4*Integral(f, 10000, (0,1), 10**-7, 4)

print(integral)
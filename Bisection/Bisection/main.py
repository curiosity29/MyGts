from sympy import *
from math import *
import sys
from Bisection import *

expr = "x^2 - 15"
a    = 3
b    = 5
eps  = 10**-5

uu = bisection_oop(a, b, eps, expr);
sol = uu.Solve();

print(f"Nghiệm của phương trình {expr} = 0 trên khoảng [{a}, {b}] với sai số {eps} là x = {sol}");
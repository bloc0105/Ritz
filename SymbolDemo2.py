import sympy as sym
import numpy as num

a,b,c = sym.symbols('a b c')

q =[a * 2,(b + c)**2,c]

z = sum(q)

print(sym.latex(z))

'''
Created on Feb 17, 2020

@author: Andy Block
'''

import sympy as sym

x,y = sym.symbols('x y')
f, g = sym.symbols('f g', cls=sym.Function)

f_x = f(x)

di = f(x).diff(x,x,y)
diff_eq = sym.Eq(f(x).diff(x,x) - 2 * f(x).diff(x) + f(x),sym.sin(x))

dorp = sym.dsolve(diff_eq, f(x))

print(sym.latex(di))
print (dorp)
print(sym.latex(dorp))
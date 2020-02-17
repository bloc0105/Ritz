'''
Created on Feb 17, 2020

@author: Andrew Block
'''


import sympy as sym

x,y,z,mu =sym.symbols('x y z mu')

expr = sym.exp(x * y * z)
expr = sym.Derivative(expr,x)

zopp = sym.Integral(sym.Integral(x + 7 * y,x),y)

mu = sym.Matrix([[1,1,1,1],[1,1,2,3]])

linear = sym.linsolve(mu,(x,y,z))

system = [sym.exp(x) - sym.sin(y),1/y - 3]

vars = [x,y]

serp = sym.nonlinsolve(system,vars)
print(sym.latex(mu))
# print (linear)
# print(sym.latex(linear))
# print (serp)

# print (sym.latex(serp))

# zeep = sym.Limit(sym.sin(x)/x,x,sym.oo)

# print (sym.latex(expr))
# print (sym.latex(zopp))
# print (sym.latex(zeep))

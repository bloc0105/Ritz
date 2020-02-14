'''
Created on Dec 17, 2019

@author: andy
'''
import sympy as sym
from sympy.core.numbers import oo #infinity

x,t,z,nu, k_1 = sym.symbols('x t z nu varrho_1')
sym.init_printing(use_unicode = True)

func = sym.sin(x) * sym.exp(x)

deri = sym.diff(func,x)

integ = sym.integrate(deri,x)

#print(deri)
#print(sym.latex(deri))
#print(integ)

zoop = sym.integrate(sym.sin(x**2),(x, -oo, oo))
#print(zoop)
#dooba = zoop.evalf
#print(dooba)
zerp = sym.limit(sym.sin(x)/x,x,0)

#print(zerp)

zapp = sym.solve(x**2 - 2,x)

print(zapp)

print(k_1)
print(sym.latex(k_1))

y = sym.Function('y')

q = sym.dsolve(sym.Eq(y(t).diff(t,t) - y(t),sym.exp(t)),y(t))
g = sym.Eq(x,nu**3)

#print(q)
#print(g)
#print(sym.latex(g))

zoob = sym.Matrix([[1,2],[2,2]]).eigenvals()
#print(zoob)


#print(func)

if __name__ == '__main__':
    pass
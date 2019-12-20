'''
Created on Dec 17, 2019

@author: andy
'''
import sympy as sym
import numpy as num


x,y,z,t = sym.symbols('x y z t')

expr = sym.cos(x) + 1

zoop = expr.subs(x,sym.exp(y))
#print(zoop)
yayyyy = zoop.evalf(subs={y:5}) # note that since this evaluation contains trig, the result is in radians
#print(yayyyy)

a = num.linspace(0,2 * num.pi, 9)
print (a)

simpleExpr = sym.sin(t)

f = sym.lambdify(t,simpleExpr,"numpy")

print(f(a))



if __name__ == '__main__':
    pass
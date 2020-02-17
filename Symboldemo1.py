'''
Created on Dec 17, 2019

@author: andy
'''
import sympy as sym

x,y,z= sym.symbols('x y z')

expr = expr = (x*y**2 - 2*x*y*z + x*z**2 + y**2 - 2*y*z + z**2)/(x**2 - 1)

expr = sym.cancel(expr)


#expr = sym.expand(expr)

print(expr)
print(sym.latex(expr))



if __name__ == '__main__':
    pass
'''
Created on Dec 17, 2019

@author: andy
'''
import sympy as sym

x,y = sym.symbols('x y')

expr = x + 2 * y

expr = expr * x

#expr = sym.expand(expr)

print(sym.latex(expr))



if __name__ == '__main__':
    pass
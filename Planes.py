'''
Created on Mar 22, 2020

@author: Andy Block

'''

import numpy as num
import matplotlib.pyplot as plott
import sympy as sym

x0 = 2
x1 = 3

y0 = 2
y1 = 3

phi0 = 1
phix = 2.6
phiy = 4


m = (phix - phi0)/(x1 - x0)
n = (phiy - phi0)/(y1 - y0)
b = phi0 - m * x0 - n * y0

print(b)



zerp = num.linspace(0,10,11)

x_Grid,Y_Grid = num.meshgrid(zerp,zerp)



vals = num.zeros([len(zerp),len(zerp)])

for counterx in range(len(x_Grid)):
    for countery in range(len(x_Grid[0])):
        vals[countery][counterx] = m * x_Grid[countery][counterx] + n * Y_Grid[countery][counterx] + b 



print(sym.latex(sym.Matrix(x_Grid)))
print(sym.latex(sym.Matrix(Y_Grid)))
print(sym.latex(sym.Matrix(vals)))


plott.contourf(x_Grid,Y_Grid,vals,150)
plott.show()



# u = m * x + n * y + b
import numpy as num
import sympy as sym
from matplotlib import pyplot as plotting

'''
Created on Dec 19, 2019

@author: andy
'''

q= num.linspace(1,2,8)

g = num.linspace(5,6,8)

X_Grid1,y_Grid1 = num.meshgrid(q,q)
X_Grid2,y_Grid2 = num.meshgrid(g,g)

print(sym.latex(sym.Matrix(X_Grid1)))
p = num.zeros([len(g),len(g)])
d = num.zeros([len(q),len(q)])

p = p + 5

d = d + 3

plotting.contourf(X_Grid1,y_Grid1,p)
plotting.contourf(X_Grid2,y_Grid2,d)

# print(p)

plotting.colorbar()
plotting.show()
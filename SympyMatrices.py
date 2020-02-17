'''
Created on Feb 17, 2020

@author: Andy Block
'''
import sympy as sym
A,B = sym.symbols('A B')

A = sym.Matrix([[4,4],[2,2]])
B = sym.Matrix([[1,2,3],[1,3,3]])

out_matrix = A * B

derp = sym.Eq(A * B, out_matrix)

# print (out_matrix)

print(sym.latex(A))
print(sym.latex(B))
print(sym.latex(out_matrix))
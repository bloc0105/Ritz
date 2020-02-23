'''
Created on Feb 20, 2020

@author: Andy Block
'''
import sympy as sym


global_matrix_size = 4
global_symbol_matrix = [[sym.symbols('k_'+ str(x+1) +str(y+1)) for x in range(global_matrix_size)] for y in range(global_matrix_size)]

print (sym.latex(sym.Matrix(global_symbol_matrix)))

delta_1, delta_2 = sym.symbols('delta_1 delta_2')

F1,F2 = sym.symbols('F_1 F_2')


zerp = sym.Matrix([[k_11, k_12],[-k_21,-k_22]])

zorp = sym.Matrix([[delta_1],[delta_2]])
zapp = [F1,F2]

doopy = zerp * zorp

solvethese = []
for equation_count in range(len(doopy)):
    solvethese.append(sym.Eq(doopy[equation_count],zapp[equation_count]))

# print(sym.latex(doopy))
# print(sym.latex(solvethese))

dapp = sym.linear_eq_to_matrix(solvethese,[delta_1,delta_2])
dopp = sym.linear_eq_to_matrix(solvethese,[F1,F2])
# print (sym.latex(dapp))
# print (sym.latex(dopp))
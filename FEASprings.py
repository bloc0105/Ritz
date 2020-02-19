'''
Created on Feb 13, 2020

@author: Andy Block
'''

import numpy as np
import sympy as sym
# Step 1 - Discretize the Structure

global_mapping = [[1,2],[2,3],[2,3],[3,4]]



k_aa,k_bb, k_ab, k_ba = sym.symbols('k_aa k_bb k_ab k_ba')
d1,d2,d3,d4 = sym.symbols('delta_1 delta_2 delta_3 delta_4')




k_11,k_12,k_13,k_14 = sym.symbols('k_11 k_12 k_13 k_14')
k_21,k_22,k_23,k_24 = sym.symbols('k_21 k_22 k_23 k_24')
k_31,k_32,k_33,k_34 = sym.symbols('k_31 k_32 k_33 k_34')
k_41,k_42,k_43,k_44 = sym.symbols('k_41 k_42 k_43 k_44')

global_position_matrix = sym.Matrix([])
global_position_matrix = global_position_matrix.row_insert(0, sym.Matrix([[k_41,k_42,k_43,k_44]]))
global_position_matrix = global_position_matrix.row_insert(0, sym.Matrix([[k_31,k_32,k_33,k_34]]))
global_position_matrix = global_position_matrix.row_insert(0, sym.Matrix([[k_21,k_22,k_23,k_24]]))
global_position_matrix = global_position_matrix.row_insert(0, sym.Matrix([[k_11,k_12,k_13,k_14]]))
final_matrix = [[0 for x in range(4)] for y in range(4)]
elemental_stiffness_matrix = [[1,-1],[-1,1]]

for element_set in global_mapping:

        stiffness_matrix = [[global_position_matrix.row(element_set[y] - 1).col(element_set[x] - 1).tolist()[0][0] * elemental_stiffness_matrix[y][x] 
                             for x in range(2)] for y in range(2)]
        print (sym.latex(stiffness_matrix))
        for rowcounter in range(2):
            for columncounter in range(2):
                
                row_pos = element_set[rowcounter] - 1
                col_pos = element_set[columncounter] - 1
                
                globally_mapped_element = stiffness_matrix[rowcounter][columncounter]
                final_matrix[row_pos][col_pos] = final_matrix[row_pos][col_pos] + globally_mapped_element

# Step 2 - Find all the element properties

# print (sym.latex(global_position_matrix))
# print (global_position_matrix.col(1).row(2).tolist()[0][0])
final_matrix = sym.Matrix(final_matrix)
print (sym.latex(final_matrix))



elemental_whole_matirix  = []

indiv_stiffness_matrices =[]


k = []
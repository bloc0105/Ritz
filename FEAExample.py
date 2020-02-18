'''
Created on Feb 13, 2020

@author: Andy Block
'''

import numpy as np
import sympy as sym
# Step 1 - Discretize the Structure

global_mapping = [];
global_mapping.append(1)
global_mapping.append(2)
global_mapping.append(3)

k_a,k_b = sym.symbols('k_a k_b')

k_11,k_12,k_13,k_14 = sym.symbols('k_11 k_12 k_13 k_14')
k_21,k_22,k_23,k_24 = sym.symbols('k_21 k_22 k_23 k_24')
k_31,k_32,k_33,k_34 = sym.symbols('k_31 k_32 k_33 k_34')
k_41,k_42,k_43,k_44 = sym.symbols('k_41 k_42 k_43 k_44')

# Step 2 - Find all the element properties
area = [32.3 * 10 ** -4,38.7 * 10 ** -4,25.8 * 10 ** -4]

intersection_angle = [90 * np.pi/180,0,135 * np.pi/180]
element_length = [2.54,2.54,2.54 * np.sqrt(2)]
modulus = [6.9 * 10**10, 20.7 * 10**10, 20.7 * 10**10]

elemental_stiffness_matrix = sym.Matrix([[1,-1],[-1,1]])

indiv_stiffness_matrices =[]


k = []

for counter in range(len(intersection_angle)):
    
    rotation_matrix = sym.Matrix([[sym.cos(intersection_angle[counter]),0],[sym.sin(intersection_angle[counter]),0],[0,sym.cos(intersection_angle[counter])],[0,sym.sin(intersection_angle[counter])]])
    k.append( (modulus[counter] * area[counter]/element_length[counter]) *  rotation_matrix * elemental_stiffness_matrix*np.transpose(rotation_matrix))
#     print(np.round(k[counter],0))
    pass


# sym.sin()



# a = sym.Matrix([[1,2],[3,4]])
# b = np.transpose(a)


# c = a * b
print(sym.latex(k))
# print (a)
# print (b)

# print(c)
# print (global_mapping)
'''
Created on Feb 13, 2020

@author: Andy Block
'''

import numpy as np
import sympy as sym
# Step 1 - Discretize the Structure

global_mapping = [[1,2],[2,3],[1,3]]
u1,v1,u2,v2,u3,v3 = sym.symbols('u1 v1 u2 v2 u3 v3')
symbol_mapping = [[u1,v1],[u2,v2],[u3,v3]]

# Step 2 - Find all the element properties
area = [32.3 * 10 ** -4,38.7 * 10 ** -4,25.8 * 10 ** -4]

intersection_angle = [90 * np.pi/180,0,135 * np.pi/180]
element_length = [2.54,2.54,2.54 * np.sqrt(2)]
modulus = [6.9 * 10**10, 20.7 * 10**10, 20.7 * 10**10]

elemental_stiffness_matrix = np.matrix([[1,-1],[-1,1]])

indiv_stiffness_matrices =[]

global_force_vector = [u1,v1,u2,v2,u3,v3]

global_stiffness_matrix = np.zeros((6,6))
k = []

uv = []

#Step 3 - Assemble the Element Matrix
for counter in range(len(global_mapping)):
    
    rotation_matrix = np.matrix([[np.cos(intersection_angle[counter]),0],[np.sin(intersection_angle[counter]),0],[0,np.cos(intersection_angle[counter])],[0,np.sin(intersection_angle[counter])]])
    k.append( (modulus[counter] * area[counter]/element_length[counter]) *  rotation_matrix * elemental_stiffness_matrix*np.transpose(rotation_matrix))
#     print(np.round(k[counter],0))
    zoop = []
    for counter2 in global_mapping[counter]:
        zoop.append(symbol_mapping[counter2 - 1][0])
        zoop.append(symbol_mapping[counter2 - 1][1])
    uv.append(zoop)
    pass
k = np.round(k)

print (global_stiffness_matrix)

for main_counter in range(len(k)):
    k_matrix = k[main_counter]
    uv_vector = uv[main_counter]
    
    for row in k_matrix:
        for column_counter in range(len(row)):
            pass

    


# np.sin()


# print(k)
# print(uv)

# a = np.matrix([[1,2],[3,4]])
# b = np.transpose(a)


# c = a * b

# print (a)
# print (b)

# print(c)
# print (global_mapping)
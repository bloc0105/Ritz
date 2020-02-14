'''
Created on Feb 13, 2020

@author: Andy Block
'''

import numpy as np
# Step 1 - Discretize the Structure

global_mapping = [];
global_mapping.append(1)
global_mapping.append(2)
global_mapping.append(3)

# Step 2 - Find all the element properties
area = [32.3 * 10 ** -4,38.7 * 10 ** -4,25.8 * 10 ** -4]

intersection_angle = [90 * np.pi/180,0,135 * np.pi/180]
element_length = [2.54,2.54,2.54 * np.sqrt(2)]
modulus = [6.9 * 10**10, 20.7 * 10**10, 20.7 * 10**10]

elemental_stiffness_matrix = np.matrix([[1,-1],[-1,1]])

indiv_stiffness_matrices =[]


k = []

for counter in range(len(intersection_angle)):
    
    rotation_matrix = np.matrix([[np.cos(intersection_angle[counter]),0],[np.sin(intersection_angle[counter]),0],[0,np.cos(intersection_angle[counter])],[0,np.sin(intersection_angle[counter])]])
    k.append( (modulus[counter] * area[counter]/element_length[counter]) *  rotation_matrix * elemental_stiffness_matrix*np.transpose(rotation_matrix))
#     print(np.round(k[counter],0))
    pass


# np.sin()



# a = np.matrix([[1,2],[3,4]])
# b = np.transpose(a)


# c = a * b

# print (a)
# print (b)

# print(c)
# print (global_mapping)
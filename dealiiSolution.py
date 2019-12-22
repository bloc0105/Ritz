'''
Created on Dec 22, 2019

@author: andy
'''

from matplotlib import pyplot as plott
import numpy as np

exes = []
whies = []
zees = []

unique_exes = []
unique_whies = []

with open("solution.txt", "r+") as f:
    for q in f.readlines():
        g = q.split()
        if len(g) == 3:
            g0 = float(g[0])
            g1 = float(g[1])
            g2 = float(g[2])
            if g1 >= 0 and g0 >=0:
                exes.append(g0)
                whies.append(g1)
                zees.append(g2)
    

for x in exes:
    if x not in unique_exes:
        unique_exes.append(x)

for y in whies:
    if x not in unique_whies:
        unique_whies.append(y)        

[X_Grid,Y_Grid] = np.meshgrid(unique_exes, unique_whies)

values = np.copy(X_Grid)
for xcounter in range(len(unique_exes)):
    for ycounter in range(len(unique_whies)):
        for counter in range(len(exes)):
            if X_Grid[ycounter][xcounter] == exes[counter] and Y_Grid[ycounter][xcounter] == whies[counter]:
                values[ycounter][xcounter] = zees[counter]

# x = np.array([1,2,3])
# y = np.array([1,2,3])
# xx, yy = np.meshgrid(exes, whies)
# X_grid = np.c_[ np.ravel(xx), np.ravel(yy) ]

# zees = zees.reshape(xx.shape)

#plott.contourf(X_Grid,Y_Grid,values,100)

#plott.colorbar()
plott.plot(np.diff(values[0]))
# print(values)
# print(X_Grid)
# print(Y_Grid)
plott.show()
if __name__ == '__main__':
    pass
'''
Created on Dec 22, 2019

@author: andy
'''

import numpy as num
import matplotlib as plotting
from matplotlib import pyplot as plott

low_boundary = 0
high_boundary = 1

number_of_divisions = (high_boundary - low_boundary) * 20 

x_range = num.linspace(low_boundary,high_boundary,number_of_divisions)
y_range = num.linspace(low_boundary,high_boundary,number_of_divisions)

X_Grid,Y_Grid = num.meshgrid(x_range,y_range)
values = num.copy(X_Grid)

f = -1

du = 0

d2udx2 = -f
d2udy2 = -f

dudx = num.copy(X_Grid)
dudy = num.copy(Y_Grid)

values[0][0] = 0 
dudx[0][0] = 0
dudy[0][0] = 0

prev_x_counter = 0
prev_y_counter = 0 

for xcounter in range(number_of_divisions):
    for ycounter in range(number_of_divisions):
        
        if xcounter == 0 and ycounter == 0:
            pass
        
        else:
            dudy[xcounter][ycounter] = dudy[xcounter][ycounter - 1] + d2udy2 * max(num.diff(y_range))
            values[xcounter][ycounter] = values[xcounter][ycounter - 1] + dudy[xcounter][ycounter] * max(num.diff(y_range))
                    
    if xcounter == 0 and ycounter == 0:
        pass
    
    else:
        dudx[xcounter][ycounter] = dudx[xcounter - 1][ycounter] + d2udx2 * max(num.diff(x_range))
        values[xcounter][ycounter] = values[xcounter - 1][ycounter] + dudx[xcounter][ycounter]  * max(num.diff(x_range))


plott.contourf(X_Grid,Y_Grid,values,100)
plott.colorbar()

# print(values)
# print(X_Grid)
# print(Y_Grid)
plott.show()

if __name__ == '__main__':
    pass
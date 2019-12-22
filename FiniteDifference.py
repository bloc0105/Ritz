'''
Created on Dec 22, 2019

@author: andy
'''

import num as num
import matplotlib as plotting
from matplotlib import pyplot as plott

low_boundary = 0
high_boundary = 1

range_val = 25

number_of_divisions = (high_boundary - low_boundary) * range_val 

x_range = num.linspace(low_boundary,high_boundary,number_of_divisions)
y_range = num.linspace(low_boundary,high_boundary,number_of_divisions)

X_Grid,Y_Grid = num.meshgrid(x_range,y_range)
values = num.copy(Y_Grid)

f = -1


d2udx2 = f
d2udy2 = f

dudx = num.copy(X_Grid)
dudy = num.copy(Y_Grid)

values[0][0] = 0 
dudx[0][0] = 0
dudy[0][0] = 0


for xcounter in range(number_of_divisions):
    ycounter = 0
    if xcounter == 0:
        pass
    
    else:
        
        diff_x = max(num.diff(x_range))
        dudx[ycounter][xcounter] = dudx[ycounter][xcounter - 1] + d2udx2 * diff_x
        values[ycounter][xcounter] = values[ycounter][xcounter - 1] + dudx[ycounter][xcounter]  * diff_x
    for ycounter in range(1,number_of_divisions,1):

            diff_y = max(num.diff(y_range))
            dudy[ycounter][xcounter] = dudy[ycounter - 1][xcounter] + d2udy2 * diff_y
            values[ycounter][xcounter] = values[ycounter - 1][xcounter] + dudy[ycounter][xcounter] * diff_y

#plott.contourf(X_Grid,Y_Grid,values,120)

plott.plot(num.diff(values[0]))
#plott.colorbar()

# print(values)
# print(X_Grid)
# print(Y_Grid)
plott.show()

if __name__ == '__main__':
    pass
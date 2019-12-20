import numpy as num
import matplotlib as plotting
from matplotlib import pyplot as plott
import sympy as sym
from array import array

x_range = num.linspace(0,1,20)
y_range = num.linspace(0,1,20)


X_Grid,Y_Grid = num.meshgrid(x_range,y_range)
values = num.copy(X_Grid)

equ_array_length = 4
x,y, f, u, a, phi = sym.symbols('x y f u a phi')
original_function = (1 - x) * (1 - y)
f = 1



array1 = [x**(n) for n in range(equ_array_length)]
array2 = [y**(n) for n in range(equ_array_length)]

trial_functions = [j * a * original_function for a in array1 for j in array2] 

sum_array1 = []
for trial_counter in range(len(trial_functions)):
    phi = 0
    for trial_counter2 in range(len(trial_functions)):
        phi = phi + sym.diff(trial_functions[trial_counter],x) * sym.diff(trial_functions[trial_counter2],x) +  sym.diff(trial_functions[trial_counter],y) * sym.diff(trial_functions[trial_counter2],y)
    sum_array1.append(a * phi + f * trial_functions[trial_counter])
    
solution_array = [sym.solve(2 * sym.integrate(sym.integrate(q,(x,0,1)),(y,0,1)),a) for q in sum_array1]
    
rounded_array  = [sym.N(z[0]) for z in solution_array]

u = 0
for d in range(len(trial_functions)):
    u = u + solution_array[d][0] * trial_functions[d] 
    pass
# q = array1 + array2

# print(solution_array)
# print(rounded_array)
# print(u)

f = sym.lambdify([x,y],u,"numpy")

for counter_x in range(len(x_range)):
    for counter_y in range(len(y_range)):
        values[counter_x][counter_y] = f(X_Grid[counter_x][counter_y],Y_Grid[counter_x][counter_y])

plott.contourf(X_Grid,Y_Grid,values,100)

# print(values)
# print(X_Grid)
# print(Y_Grid)
plott.show()

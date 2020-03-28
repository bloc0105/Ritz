import numpy as num
import matplotlib as plotting
from matplotlib import pyplot as plott
import sympy as sym

low_boundary = 0
high_boundary = 1

number_of_divisions = (high_boundary - low_boundary) * 20 

x_range = num.linspace(low_boundary,high_boundary,number_of_divisions)
y_range = num.linspace(low_boundary,high_boundary,number_of_divisions)


X_Grid,Y_Grid = num.meshgrid(x_range,y_range)
values = num.copy(X_Grid)

equ_array_length = 3
x,y, f, u,a, phi = sym.symbols('x y f u a phi')
original_function = (1 - x) * (1 - y) 
f = -1



array1 = [x**(n) for n in range(equ_array_length)]
array2 = [y**(n) for n in range(equ_array_length)]

trial_coefficients = [sym.symbols('a_' + str(counter)) for counter in range(equ_array_length**2)]


trial_functions = [j * k * original_function for k in array1 for j in array2] 

sum_array1 = []
for functional_counter in range(len(trial_functions)):
    phi = 0
    for summation_counter in range(len(trial_functions)):
        phi = phi + (sym.diff(trial_functions[functional_counter],x) * sym.diff(trial_functions[summation_counter],x) +  sym.diff(trial_functions[functional_counter],y) * sym.diff(trial_functions[summation_counter],y)) * trial_coefficients[summation_counter]
    sum_array1.append(phi + f * trial_functions[functional_counter])
functional_array = []
for counter in range(len(sum_array1)):
    print('Integral - ' + str(counter))
    
    functional_array.append(sym.Eq(2 * sym.N(sym.integrate(sym.integrate(sum_array1[counter],(x,low_boundary,high_boundary)),(y,low_boundary,high_boundary)),7),0))



print(sym.latex(sym.Matrix(functional_array)))

result_set = sym.linsolve(functional_array, trial_coefficients)

u = 0
for d in range(len(result_set.args[0])):
    u = u + result_set.args[0][d] * trial_functions[d] 
    pass

f = sym.lambdify([x,y],u,"numpy")

for counter_x in range(len(x_range)):
    for counter_y in range(len(y_range)):
        values[counter_x][counter_y] = f(X_Grid[counter_x][counter_y],Y_Grid[counter_x][counter_y])
        

plott.contourf(X_Grid,Y_Grid,values, 120)
plott.colorbar()

plott.show()

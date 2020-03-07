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

equ_array_length = 4
x,y, f, u,a, phi = sym.symbols('x y f u a phi')
original_function = (1 - x) * (1 - y)
f = -1



array1 = [x**(n) for n in range(equ_array_length)]
array2 = [y**(n) for n in range(equ_array_length)]

trial_coefficients = [sym.symbols('a_' + str(counter)) for counter in range(equ_array_length**2)]
print(sym.latex(sym.Matrix(trial_coefficients)))

trial_functions = [j * k * original_function for k in array1 for j in array2] 

sum_array1 = []
for functional_counter in range(len(trial_functions)):
    phi = 0
    for summation_counter in range(len(trial_functions)):
        phi = phi + sym.diff(trial_functions[functional_counter],x) * sym.diff(trial_functions[summation_counter],x) +  sym.diff(trial_functions[functional_counter],y) * sym.diff(trial_functions[summation_counter],y)
        sum_array1.append(trial_coefficients[functional_counter] * phi + f * trial_functions[functional_counter])
        trial_function_counter += 1
    
solution_array = [sym.solve(2 * sym.integrate(sym.integrate(q,(x,low_boundary,high_boundary)),(y,low_boundary,high_boundary)),a) for q in sum_array1]
# print(sym.latex(sym.Matrix(solution_array)))
 
rounded_array  = [sym.N(z[0]) for z in solution_array]
print(sym.latex(sym.Matrix(rounded_array)))


u = 0
for d in range(len(trial_functions)):
    u = u + solution_array[d][0] * trial_functions[d] 
    pass

f = sym.lambdify([x,y],u,"numpy")

for counter_x in range(len(x_range)):
    for counter_y in range(len(y_range)):
        values[counter_x][counter_y] = f(X_Grid[counter_x][counter_y],Y_Grid[counter_x][counter_y])
        

# plott.plot(num.diff(values[0]))
plott.contourf(X_Grid,Y_Grid,values, 120)
plott.colorbar()

# print(values)
# print(X_Grid)
# print(Y_Grid)
plott.show()

import numpy as num
import matplotlib as plotting
from matplotlib import pyplot as plott
import sympy as sym

low_boundary = 0
high_boundary = 1

y0,y1,x0,x1 = sym.symbols('y_0 y_1 x_0 x_1')

number_of_divisions = (high_boundary - low_boundary) * 20 

x_range = num.linspace(low_boundary,high_boundary,number_of_divisions)
y_range = num.linspace(low_boundary,high_boundary,number_of_divisions)


X_Grid,Y_Grid = num.meshgrid(x_range,y_range)
values = num.copy(X_Grid)

equ_array_length = 3
x,y, f, u,a, phi = sym.symbols('x y f u a phi')
original_function = (1 - x) * (1 - y) 
phi_solved = sym.Derivative(phi,x,x) + sym.Derivative(phi,y,y)
# print(sym.latex(phi_solved.doit()))
f = -1



array1 = [x**(n) for n in range(equ_array_length)]
array2 = [y**(n) for n in range(equ_array_length)]

trial_coefficients = [sym.symbols('a_' + str(counter)) for counter in range(equ_array_length**2)]
trial_functions = [j * k * original_function for k in array1 for j in array2]
trial_phi = [trial_coefficients[counter] * trial_functions[counter] for counter in range(len(trial_coefficients))]

coefficients = []
# print(sym.latex(sym.Matrix(trial_phi)))
for counter  in range(len(trial_phi)):
    residual = phi_solved.subs(phi,trial_phi[counter]).doit() - f
    current_var = trial_coefficients[counter]
    weighted_average = sym.integrate(sym.integrate(residual * trial_functions[counter],(x,low_boundary,high_boundary)),(y,low_boundary,high_boundary))
#     print(sym.latex(weighted_average))
    ansewr = sym.solve(weighted_average,current_var)
    coefficients.append(ansewr)
    
    pass


# print(sym.latex(sym.Matrix(coefficients)))
# print(sym.latex(sym.Matrix(trial_functions)))

print(len(coefficients))
print(len(trial_functions))

for counter in range(len(coefficients)):
#     print(coefficients[counter][0] * trial_functions[counter])
    print(sym.latex(coefficients[counter][0]))
    print(trial_functions[counter])
f = sym.lambdify([x,y],u,"numpy")
 
for counter_x in range(len(x_range)):
    for counter_y in range(len(y_range)):
        values[counter_x][counter_y] = f(X_Grid[counter_x][counter_y],Y_Grid[counter_x][counter_y])
         
 
plott.contourf(X_Grid,Y_Grid,values, 120)
plott.colorbar()
 
plott.show()

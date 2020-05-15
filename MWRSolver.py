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
# print(sym.latex(phi_solved))
# print(sym.latex(phi_solved.doit()))
f = -1

array1 = [x**(n + 1) for n in range(equ_array_length)]
array2 = [y**(n + 1) for n in range(equ_array_length)]

# array1 = [x**(n) for n in range(equ_array_length)]
# array2 = [y**(n) for n in range(equ_array_length)]

trial_coefficients = [sym.symbols('a_' + str(counter)) for counter in range(equ_array_length**2)]
trial_functions = [j * k * original_function for k in array1 for j in array2]
trial_equations = [trial_coefficients[counter] * trial_functions[counter] for counter in range(len(trial_coefficients))]
trial_phi = sum(trial_equations) 
trial_phi += original_function

# print(sym.latex(sym.Matrix(trial_functions)))
# print(sym.latex(trial_phi))
coefficients = []
eqs = []
residual = phi_solved.subs(phi,trial_phi).doit() - f

# print(sym.latex(residual))
solution_equations = []
for counter  in range(len(trial_equations)):
    print("Computation " + str(counter) + " of " + str(len(trial_equations) - 1))

    innerWight = residual * trial_functions[counter]
#     print(innerWight)
#     print('--------------------------------------------------------------------------------')
    weighted_average = sym.integrate(sym.integrate(innerWight,(x,low_boundary,high_boundary)),(y,low_boundary,high_boundary))
    solution_equations.append(weighted_average)
    
# print(sym.latex(sym.Matrix(solution_equations)))

result_set = sym.linsolve(solution_equations,trial_coefficients)
print(sym.latex(result_set))
# print(sym.latex(sym.Matrix(coefficients)))
# print(sym.latex(sym.Matrix(trial_functions)))

# print(coefficients)
# print(len(trial_functions))

# u = original_function
u = 0
for d in range(len(result_set.args[0])):
    u = u + result_set.args[0][d] * trial_functions[d]  

# print(sym.latex(u))

f = sym.lambdify([x,y],u,"numpy")
 
for counter_x in range(len(x_range)):
    for counter_y in range(len(y_range)):
        values[counter_x][counter_y] = f(X_Grid[counter_x][counter_y],Y_Grid[counter_x][counter_y])
         
 
plott.contourf(X_Grid,Y_Grid,values, 120)
plott.colorbar()
 
plott.show()

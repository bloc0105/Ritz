import numpy as num
import matplotlib as plotting
from matplotlib import pyplot as plott
import sympy as sym

low_boundary = 0 # The lower boundary of the graph
high_boundary = 1 # The upper boundary of the graph
grid_points = 4 # The number of nodes that exist along the axes of the boundary (inclusive)

number_of_divisions = (high_boundary - low_boundary) * grid_points #each quadrant will have the same number of nodes in it 




x_range = num.linspace(low_boundary,high_boundary,number_of_divisions)
y_range = num.linspace(low_boundary,high_boundary,number_of_divisions)


elements_with_coordinates = [ [[[x_range[xcounter],x_range[xcounter + 1],y_range[ycounter],y_range[ycounter + 1]]
                for xcounter in range(number_of_divisions - 1) ] for ycounter in range(number_of_divisions - 1)]]
# print(num.asarray(elements_with_coordinates))
element_phi_symbols = [[sym.symbols('phi_' + str(xcounter) + str(ycounter)) for xcounter in range(number_of_divisions - 1) ] for ycounter in range(number_of_divisions - 1)]


equ_array_length = 5
x,y, f, u, a = sym.symbols('x y f u a')
original_function = (1 - x) * (1 - y)
f = -1



array1 = [x**(n) for n in range(equ_array_length)]
array2 = [y**(n) for n in range(equ_array_length)]

trial_functions = [j * k * original_function for k in array1 for j in array2] 

sum_array1 = []
for xcounter in range(number_of_divisions - 1): 
    for ycounter in range(number_of_divisions - 1): 
        for trial_counter in range(len(trial_functions)):
            phi = 0
            for trial_counter2 in range(len(trial_functions)):
                phi = phi + sym.diff(trial_functions[trial_counter],x) * sym.diff(trial_functions[trial_counter2],x) +  sym.diff(trial_functions[trial_counter],y) * sym.diff(trial_functions[trial_counter2],y)
            sum_array1.append(a * phi + f * trial_functions[trial_counter])
    
solution_array = [sym.solve(2 * sym.integrate(sym.integrate(q,(x,low_boundary,high_boundary)),(y,low_boundary,high_boundary)),a) for q in sum_array1]
 
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

X_Grid,Y_Grid = num.meshgrid(x_range,y_range)
values = num.copy(X_Grid)
for counter_x in range(len(x_range)):
    for counter_y in range(len(y_range)):
        values[counter_x][counter_y] = f(X_Grid[counter_x][counter_y],Y_Grid[counter_x][counter_y])
        

plott.plot(num.diff(values[0]))
#plott.colorbar()

# print(values)
# print(X_Grid)
# print(Y_Grid)
# plott.show()

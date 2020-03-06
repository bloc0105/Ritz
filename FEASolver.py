import numpy as num
import matplotlib as plotting
from matplotlib import pyplot as plott
import sympy as sym

low_boundary = 0 # The lower boundary of the graph
high_boundary = 1 # The upper boundary of the graph
grid_points = 10 # The number of nodes that exist along the axes of the boundary (inclusive)
equ_array_length = 1
sub_grid_size = 5

number_of_divisions = (high_boundary - low_boundary) * grid_points #each quadrant will have the same number of nodes in it 

x_range = num.linspace(low_boundary,high_boundary,number_of_divisions)
y_range = num.linspace(low_boundary,high_boundary,number_of_divisions)


elements_with_coordinates =  [ [[x_range[xcounter],x_range[xcounter + 1],y_range[ycounter],y_range[ycounter + 1]] for xcounter in range(number_of_divisions - 1)] for ycounter in range(number_of_divisions - 1)]
# print(elements_with_coordinates)
solution_matrix = [[0 for xcounter in range(number_of_divisions - 1) ] for ycounter in range(number_of_divisions - 1)]
# print(num.asarray(elements_with_coordinates))
# element_phi_symbols = [[sym.symbols('phi_' + str(xcounter) + str(ycounter)) for xcounter in range(number_of_divisions - 1) ] for ycounter in range(number_of_divisions - 1)]



x,y, f, u, a, g = sym.symbols('x y f u a g')
original_function = (1 - x) * (1 - y)
f = -1



array1 = [x**(n) for n in range(equ_array_length)]
array2 = [y**(n) for n in range(equ_array_length)]

trial_functions = [j * k * original_function for k in array1 for j in array2] 


for xcounter in range(number_of_divisions - 1): 
    for ycounter in range(number_of_divisions - 1): 
        sum_array1 = []
        for trial_counter in range(len(trial_functions)):
            phi = 0
            for trial_counter2 in range(len(trial_functions)):
                phi = phi + sym.diff(trial_functions[trial_counter],x) * sym.diff(trial_functions[trial_counter2],x) +  sym.diff(trial_functions[trial_counter],y) * sym.diff(trial_functions[trial_counter2],y)
            sum_array1.append(a * phi + f * trial_functions[trial_counter])
        print(str(xcounter) + ' - ' + str(ycounter))
        solution_array = [sym.solve(2 * sym.integrate(sym.integrate(q,(x,elements_with_coordinates[xcounter][ycounter][0],elements_with_coordinates[xcounter][ycounter][1])),(y,elements_with_coordinates[xcounter][ycounter][2],elements_with_coordinates[xcounter][ycounter][3])),a) for q in sum_array1]
 
        rounded_array  = [sym.N(z[0]) for z in solution_array]


        u = 0
        for d in range(len(trial_functions)):
            u += solution_array[d][0] * trial_functions[d] 
        solution_matrix[xcounter][ycounter] += u


X_Grid = []
Y_Grid = []
values = []
x_plot = []
vals_plot = []

for xcounter in range(len(elements_with_coordinates)): 
    for ycounter in range(len(elements_with_coordinates[0])):
#         x_lin = num.linspace(elements_with_coordinates[xcounter][ycounter][0], elements_with_coordinates[xcounter][ycounter][1], sub_grid_size, endpoint=(xcounter == len(elements_with_coordinates) - 1))
#         y_lin = num.linspace(elements_with_coordinates[xcounter][ycounter][2], elements_with_coordinates[xcounter][ycounter][3], sub_grid_size, endpoint=(ycounter == len(elements_with_coordinates[0]) - 1))
        x_lin = num.linspace(elements_with_coordinates[xcounter][ycounter][0], elements_with_coordinates[xcounter][ycounter][1], sub_grid_size, endpoint=True)
        y_lin = num.linspace(elements_with_coordinates[xcounter][ycounter][2], elements_with_coordinates[xcounter][ycounter][3], sub_grid_size, endpoint=True)
        g = sym.lambdify([x,y],solution_matrix[xcounter][ycounter],"numpy")
        for list_counter_x in range(len(x_lin)):
            X_list = []
            Y_list = []
            
            
            for list_counter_y in range(len(y_lin)):
                if (y_lin[list_counter_y] == 0):
                    x_plot.append(x_lin[list_counter_x])
                    vals_plot.append(g(x_lin[list_counter_x], 0))
                
                X_list.append(x_lin[list_counter_x])
                Y_list.append(y_lin[list_counter_y])
            X_Grid.append(X_list)
            Y_Grid.append(Y_list)
            print(X_list)
            print(Y_list)
            
            
            vals = [g(X_list[counter], Y_list[counter]) for counter in range(len(X_list))]
            print(vals)
            values.append(vals)    
            print('-----------------------------------------')   
plott.subplot(121)
plott.plot(x_plot,vals_plot, '-o')
plott.subplot(122)
plott.contourf(X_Grid,Y_Grid,values, 120)
plott.colorbar()
plott.show()

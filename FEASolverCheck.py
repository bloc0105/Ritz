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

element_range = number_of_divisions - 1


XGrid,YGrid = num.meshgrid(x_range, y_range)

z_matrix = num.zeros([number_of_divisions,number_of_divisions])


# print(sym.latex(sym.Matrix(XGrid)))
# print(sym.latex(sym.Matrix(YGrid)))
# print(sym.latex(sym.Matrix(z_matrix)))

x,y, f, u = sym.symbols('x y f u')
x0,y0,y1,x1,phi0,phix,phiy = sym.symbols('x_00 y_00 y_y0 x_0x phi_00 phi_0x phi_y0')

m = (phix - phi0)/(x1 - x0)
n = (phiy - phi0)/(y1 - y0)
b = phi0 - m * x0 - n * y0

phi_final = m * x + n * y + b




solution_list = []
var_list= []
for xcounter in range(number_of_divisions): 
    for ycounter in range(number_of_divisions):
        if xcounter < element_range and ycounter < element_range:
            x_max = XGrid[ycounter][xcounter + 1]
            x_min = XGrid[ycounter][xcounter]
        
            y_max = YGrid[ycounter +1][xcounter]
            y_min = YGrid[ycounter][xcounter]
            
            the_plane = phi_final.subs(x0,x_min).subs(x1,x_max).subs(y0,y_min).subs(y1,y_max).subs(phi0,0).subs(phix,1).subs(phiy,2)
            solution_list.append(the_plane)
            
            the_func = sym.lambdify([x,y],the_plane,"numpy")
            z_matrix[ycounter][xcounter] = the_func(XGrid[ycounter][xcounter],YGrid[ycounter][xcounter])
 
print(sym.latex(sym.Matrix(solution_list)))
plott.contourf(XGrid,YGrid,z_matrix,200)
plott.colorbar()
plott.show()

    

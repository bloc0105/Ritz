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

z_matrix = [[sym.symbols('z_x' + str(counterx) + '_y' + str(countery)) for counterx in range(len(XGrid))] for countery in range(len(XGrid[0]))]

# print(sym.latex(sym.Matrix(XGrid)))
# print(sym.latex(sym.Matrix(YGrid)))
# print(sym.latex(sym.Matrix(z_matrix)))

x,y, f, u = sym.symbols('x y f u')
x0,y0,y1,x1,phi0,phix,phiy = sym.symbols('x_0 y_0 y_1 x_1 phi_0 phi_x phi_y')

m = (phix - phi0)/(x1 - x0)
n = (phiy - phi0)/(y1 - y0)
b = phi0 - m * x0 - n * y0

phi_final = m * x + n * y + b

# print(sym.latex(phi_final))


# print(sym.latex(sym.diff(phi_final,x)))
f = -1
phi = (sym.diff(phi_final,x))**2  +  (sym.diff(phi_final,y))**2 + 2 * f * phi_final

print(sym.latex(phi))

phi_integral = sym.integrate(sym.integrate(phi,(x,x0,x1)),(y,y0,y1))

# print(sym.latex(phi_integral))


phi_diff_0 = sym.diff(phi_integral, phi0)
phi_diff_x = sym.diff(phi_integral, phix)
phi_diff_y = sym.diff(phi_integral, phiy)

phi_solves = []

phi_solves.append(sym.simplify(phi_diff_0))
phi_solves.append(sym.simplify(phi_diff_x))
phi_solves.append(sym.simplify(phi_diff_y))

# print(sym.latex(sym.Matrix(phi_solves)))

q = sym.linsolve(phi_solves,[phi0,phix,phiy])

# print(sym.latex(sym.Matrix([phi0,phix,phiy])))

# print(sym.latex(q))
 
# print (sym.latex(phi_final))
# print (sym.latex(phi))
# print (sym.latex(phi_integral))
# print(sym.latex(phi_diff_y))

solution_list = []
var_list= []
for xcounter in range(number_of_divisions): 
    for ycounter in range(number_of_divisions):
        if xcounter < element_range and ycounter < element_range:
            x_max = XGrid[ycounter][xcounter + 1]
            x_min = XGrid[ycounter][xcounter]
        
            y_max = YGrid[ycounter +1][xcounter]
            y_min = YGrid[ycounter][xcounter]
        
            z_y = z_matrix[ycounter + 1][xcounter]
            z_x = z_matrix[ycounter][xcounter + 1]
            z_0 = z_matrix[ycounter][xcounter]
        
            phi_subs_0 = sym.N(phi_diff_0.subs(x0,x_min).subs(x1,x_max).subs(y0,y_min).subs(y1,y_max).subs(phi0,z_0).subs(phix,z_x).subs(phiy,z_y),2)
            phi_subs_x = sym.N(phi_diff_x.subs(x0,x_min).subs(x1,x_max).subs(y0,y_min).subs(y1,y_max).subs(phi0,z_0).subs(phix,z_x).subs(phiy,z_y),2)
            phi_subs_y = sym.N(phi_diff_y.subs(x0,x_min).subs(x1,x_max).subs(y0,y_min).subs(y1,y_max).subs(phi0,z_0).subs(phix,z_x).subs(phiy,z_y),2)
        
            solution_list.append(sym.Eq(phi_subs_0,0))
            solution_list.append(sym.Eq(phi_subs_x,0))
            solution_list.append(sym.Eq(phi_subs_y,0))
            
            var_list.append([x_min,x_max,y_min,y_max, z_0,z_x,z_y])
        else:
            solution_list.append(sym.Eq(z_matrix[ycounter][xcounter],0))
        
#         print(sym.latex(sym.Matrix([x_min,x_max,y_min,y_max,z_0,z_x,z_y])))


z_var_list = [z_matrix[counterx][countery] for counterx in range(len(z_matrix)) for countery in range(len(z_matrix[0]))]       

# print(sym.latex(sym.Matrix(z_var_list))) 
# print(sym.latex(sym.Matrix(var_list)))    
   

# print(len(z_var_list))
# print(sym.latex(sym.Matrix(z_var_list)))
# print(len(solution_list))

print(sym.latex(sym.Matrix(solution_list)))

resultset = sym.linsolve(solution_list,z_var_list)

# print(sym.latex(resultset))


# plott.colorbar()
# plott.show()

    

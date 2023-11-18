import pde
from pde import DiffusionPDE, ScalarField, UnitGrid, MemoryStorage, movie, PDE, CartesianGrid, VectorField, CylindricalSymGrid
import numpy as np

x_length = 100
y_length = 100
z_length = 20
delta_x  = 4
delta_y = 4
delta_z = 2
N_x = int(x_length/delta_x)
N_y = int(y_length/delta_y)
N_z = int(z_length/delta_z)

grid = CartesianGrid([[0, x_length], [0, y_length], [0, z_length]], [N_x, N_y, N_z])  # generate grid
#grid = UnitGrid([16, 16, 16])
Pres = 200 * 10**5
p0 = np.ones((N_x, N_y, N_z))*Pres # начальное пластовое давление случайное

#bc = {'value': Pres}
bc_x_left = {"value": 500*10**5}
bc_x_right = {"value": 200*10**5}
bc_x = [bc_x_left, bc_x_right]
#bc_y = {"derivative": 0}
bc_y = {"value": 200*10**5}
bc_z = {"value": 200*10**5}
bc = [bc_x, bc_y, bc_z]
field = pde.ScalarField(grid, data=p0)
#field.laplace(bc)
#field.insert([100, 100], 5*10**10)
perm = 10**(-14)
k = np.ones((N_x, N_y, N_z))*perm
#k[5:15, 10:12, 7:9] = 100 * perm
k_field = ScalarField(grid, data=k)
fi = 0.2
cf = 10**(-9)
mu = 2*10**(-1)
diff = k/2/fi/mu
#diff = "1+tanh(x)"
diff_number = perm/fi/cf/mu
diff_field = ScalarField(grid, data=diff)
delta_t = 6

courant_number_min = 2*delta_t*np.min(diff)/delta_x**2
courant_number_max = 2*delta_t*np.max(diff)/delta_x**2
print(courant_number_min, courant_number_max)
point = ScalarField(grid, data=0)
point.insert([10, 10, 5], 1) # скважина
#source = f"sin(2*pi*{10}*t)" # давление на скважине как то меняется или...
source = f"30000" # давление на скважине 100 -- или расход? что это??? :)
#eq = PDE({'u': f"diff*laplace(u) + {source}*point"}, bc=bc, consts={'point':point, 'diff': diff_field}) # уравнение диффузии с источником
#eq = PDE({'u': f"diff*laplace(u) "}, bc=bc, bc_ops={'u:laplace': {[100, 100]:10**5}}, consts={'diff': diff_field}) # уравнение диффузии с источником
eq = PDE({'u': f"diff*laplace(u**2)"}, bc=bc, consts={'diff': diff_field}) # уравнение диффузии с источником

gradient_diff = np.gradient(diff)
gradient_diff_field = VectorField(grid, data=gradient_diff)
term_1 = f"diff*laplace(p)"
term_2 = f"dot(gradient_diff_field, gradient(p))"
#eq = DiffusionPDE(0.1)
#eq = PDE({'p': f"{term_1}+{term_2}"}, bc=bc, consts={"diff":diff_field, 'gradient_diff_field':gradient_diff_field}) # уравнение диффузии
#term_2 = f"dot(gradient(diff), gradient(p))"
#eq = PDE({'p': f"{term_1}+{term_2}"}, bc=bc, consts={"diff":diff_field}) # уравнение диффузии
#eq = PDE({'p': f"{term_1}"}, bc=bc, consts={"diff":diff_field}) # уравнение диффузии

#eq = DiffusionPDE(diffusivity= diff, bc =bc)
result = eq.solve(field, t_range=60*60*500, dt = delta_t)
result.plot(cmap='magma')
import numpy as np
import scipy.sparse as sps
import porepy as pp

Nx = Ny = 20
phys_dims = [1,1]
g = pp.CartGrid([Nx, Ny], phys_dims)
g.compute_geometry()

# Permeability
perm = pp.SecondOrderTensor(np.ones(g.num_cells))
print(perm.values)

# Unitary scalar source already integrated in each cell
f = g.cell_volumes*5
print('f=', f)

# Boundary conditions
b_faces = g.tags['domain_boundary_faces'].nonzero()[0]
bc = pp.BoundaryCondition(g, b_faces, ['dir']*b_faces.size)
bc_val = np.ones(g.num_faces)
print(bc_val)

# Collect all parameters in a dictionary
parameters = {"second_order_tensor": perm, "source": f, "bc": bc, "bc_values": bc_val}
print(parameters)

data_key = "flow"
data = pp.initialize_default_data(g, {}, data_key, parameters)

flow_discretization = pp.Tpfa(data_key)
flow_discretization.discretize(g, data)
A, b_flow = flow_discretization.assemble_matrix_rhs(g, data)

rhs_discretization = pp.ScalarSource(data_key)
rhs_discretization.discretize(g, data)
_, b_rhs = rhs_discretization.assemble_matrix_rhs(g, data)

p_tpfa = sps.linalg.spsolve(A, b_flow+b_rhs)
pp.plot_grid(g, p_tpfa, figsize=(15, 12), plot_2d=True)
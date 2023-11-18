import numpy as np
import porepy as pp

# Point coordinates
point_coordinates = np.array([[0, 2, 1, 1],
                              [0, 0, 0, 1]])

# Point connections as a 2 x num_frac array
point_indices = np.array([[0, 2],
                          [1, 3]])

# Define a rectangular domain in terms of range in the two dimensions
domain = {'xmin': -2, 'xmax': 3, 'ymin': -2, 'ymax': 3}
# Define a fracture network in 2d
network_2d = pp.FractureNetwork2d(point_coordinates, point_indices, domain)

# Set preferred mesh size close to the fracture and at the boundary (essentially this is a far-field value)
mesh_args = {'mesh_size_frac': 0.3, 'mesh_size_bound': 1.0}

# Generate a mixed-dimensional grid
mdg = network_2d.mesh(mesh_args)
network_2d.plot()
pp.plot_grid(mdg, figsize=(8,8))
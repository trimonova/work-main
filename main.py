 # This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
import porepy as pp
import inspect

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.
    nx = np.array([3, 2])
    g = pp.CartGrid(nx)
    print(g)
    # Number of cells,
    print(g.num_cells)
    # faces
    print(g.num_faces)
    # and nodes
    print(g.num_nodes)
    # And the grid's dimension
    print(g.dim)
    print(g.nodes)
    g.compute_geometry()
    print(g.cell_centers)
    import matplotlib
    cell_id = np.arange(g.num_cells)
    pp.plot_grid(g, cell_value=cell_id, info='c', alpha=0.5, figsize=(15, 12))


 # Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

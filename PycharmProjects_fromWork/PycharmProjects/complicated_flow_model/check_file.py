from scipy.sparse import coo_matrix, linalg, hstack, vstack
import numpy as np
from scipy.sparse.linalg import spsolve
from shapely.geometry import Polygon
A = np.eye(5)*5
A_coo = coo_matrix(A)
print(hstack([A_coo, A]))

B = np.ones((5,1))
B_coo = coo_matrix(B)

P = np.linalg.solve(A,B)
print(P)

Ainv_coo = linalg.inv(A_coo)
P_coo = Ainv_coo.multiply(B)
print(P_coo)

A_csr = A_coo.tocsr()
B_csr = B_coo.tocsr()
P_csr = spsolve(A_csr, B)
print(P_csr)

all_cols = np.arange(A_csr.shape[1])
print(all_cols)
cols_to_keep = np.where(np.logical_not(np.in1d(all_cols, [1])))[0]
m = A_csr[:, cols_to_keep]
M = m[cols_to_keep, :]
b = np.delete(B, 1)
print(M.toarray())
P = spsolve(M, b)
print(P)

import numpy as np
import matplotlib.pyplot as plt

phi = np.random.randn(20, 20) # initial value for phi
F = 1
dt = 1
it = 100

# for i in range(it):
#     dphi = np.gradient(phi)
#     dphi2 = [i**2 for i in dphi]
#     print(type(dphi))
#     dphi_norm = np.sqrt(np.sum(dphi2, axis=0))
#
#     phi = phi + dt * F * dphi_norm
#
#     # plot the zero level curve of phi
#     plt.contour(phi, 0)
#     plt.show()

dict1 = {(1,4):1, (1,3):2, (2,3):2, (2,4):2}
keys_list = list(dict1.keys())
print(sorted(keys_list))
def func(temp):
    return temp[1], temp[0]
keys_list.sort(key=func)
print(keys_list)
print(np.sin(90))
print(np.sin(np.pi/2))

from scipy.spatial import ConvexHull
x_oil_coords = [0.009652156229659955, 0.008099121039640365, 0.009575555538987226, 0.008034845121081742]
y_oil_coords = [0.008099120732373153, 0.009652156595846758, 0.00803484512108174, 0.009575555538987226]
oil_coords_zip = zip([0.009652156229659955, 0.008099121039640365, 0.009575555538987226, 0.008034845121081742], [0.008099120732373153, 0.009652156595846758, 0.00803484512108174, 0.009575555538987226])
pgon_oil = Polygon(oil_coords_zip)
print(pgon_oil.area)
print(oil_coords_zip)
oil_coords = [(x_oil_coords[i], y_oil_coords[i]) for i in range(len(x_oil_coords))]
hull = ConvexHull(oil_coords).vertices
print(hull)
sort_bound_coords = [oil_coords[i] for i in hull]
print(sort_bound_coords)
print(Polygon(sort_bound_coords).area)
cell_coords = zip([0.009575555538987226, 0.009652159983299123, 0.008099123882050396, 0.008034845121081742, 0.009575555538987226], [0.00803484512108174, 0.008099123882050394, 0.009652159983299123, 0.009575555538987226, 0.00803484512108174])

oil_area_coords = [(0.0125000000000002, 2.7925268031909263), (0.0125, 2.8803352935086544), (0.0125, 2.7925268031909263)]
x_oil_area_coords = [i[0] * np.cos(i[1]) for i in oil_area_coords]
y_oil_area_coords = [i[0] * np.sin(i[1]) for i in oil_area_coords]
xy_oil_area_coords = [(x_oil_area_coords[i], y_oil_area_coords[i]) for i in
                                          range(len(x_oil_area_coords))]
hull = ConvexHull(xy_oil_area_coords).vertices
print(hull)

func_matrix_cell = [11, 12, 10**(-11), 13]
print(func_matrix_cell)
for i in range(len(func_matrix_cell)):
    if abs(func_matrix_cell[i]) < 10 ** (-10):
        func_matrix_cell[i] = 0
print(func_matrix_cell)

a = 2.1729349187329718e-07
print(round(a, 20))

a = [(1,1),(2,1),(3,3),(4,1)]
#a.append(a[0])
print(a)
pgon_oil = Polygon(a)  # Assuming the OP's x,y coordinates
oil_area = pgon_oil.area # площадь многоугольника
print(oil_area)
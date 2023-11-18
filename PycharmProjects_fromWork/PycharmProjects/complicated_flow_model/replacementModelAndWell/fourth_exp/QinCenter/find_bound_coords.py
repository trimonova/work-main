from shapely.geometry import Polygon
import numpy as np

def find_bound_coords(func_matrix, coord_matrix, delta_r_list, delta_fi_list, M_fi_full, N_r_full):
    bound_coords_rad = []
    for m in range(M_fi_full-1):
        for n in range(N_r_full - 1):

            if func_matrix[n][m]>0 and func_matrix[n+1][m]>0 and func_matrix[n][m+1]>0 and func_matrix[n+1][m+1]>0:
                continue
            elif func_matrix[n][m]<0 and func_matrix[n+1][m]<0 and func_matrix[n][m+1]<0 and func_matrix[n+1][m+1]<0:
                continue
            else:

                func_matrix_cell = [func_matrix[n][m], func_matrix[n+1][m], func_matrix[n+1][m+1], func_matrix[n][m+1], func_matrix[n][m]]
                func_matrix_cell_signs = [np.sign(i) for i in func_matrix_cell]
                coord_matrix_cell = [coord_matrix[n][m], coord_matrix[n+1][m], coord_matrix[n+1][m+1], coord_matrix[n][m+1], coord_matrix[n][m]]

                for i in range(4):
                    if func_matrix_cell_signs[i] != func_matrix_cell_signs[i+1] and func_matrix_cell_signs[i] != 0 and func_matrix_cell_signs[i+1] != 0:
                        # значит граница вытеснения пересекает эту сторону. Найдем ее координаты
                        proport_coef = abs(func_matrix_cell[i] / func_matrix_cell[i + 1])
                        if i == 0:
                            bound_coords_rad.append((coord_matrix_cell[i][0]+proport_coef*delta_r_list[n]/(proport_coef+1), coord_matrix_cell[i][1]))
                        if i == 1:
                            bound_coords_rad.append((coord_matrix_cell[i][0], coord_matrix_cell[i][1]+delta_fi_list[m]/(proport_coef+1)))
                        if i == 2:
                            bound_coords_rad.append((coord_matrix_cell[i][0]-proport_coef*delta_r_list[n]/(proport_coef+1), coord_matrix_cell[i][1]))
                        if i == 3:
                            bound_coords_rad.append(
                                (coord_matrix_cell[i][0], coord_matrix_cell[i][1] - delta_fi_list[m] / (proport_coef + 1)))

                for i in range(4):
                    if func_matrix_cell_signs[i] == 0:
                        bound_coords_rad.append((coord_matrix_cell[i]))

    bound_coords_cart = []
    for elem in bound_coords_rad:
        r = elem[0]
        fi = elem[1]
        bound_coords_cart.append((r*np.cos(fi), r*np.sin(fi)))

    return bound_coords_rad, bound_coords_cart





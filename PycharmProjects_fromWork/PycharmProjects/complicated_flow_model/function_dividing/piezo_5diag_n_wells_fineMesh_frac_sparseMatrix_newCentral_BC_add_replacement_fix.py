# Решение уравнения пьезопроводности в цилиндрических координатах, двумерный случай (неявная схема)
# В центре скважина с постоянным давлением P, на границах задается: градиент давления равен нулю.
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve
from start_to_do_replacement import replace_boundary

def PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, P_center, wells_frac_coords,
                         Func_matrix, delta_r_list, delta_fi_list, frac_angle_cell, frac_angle_2_cell):

    B = np.zeros((N_r_full*M_fi_full, 1))
    for m in range(M_fi_full):

        A = np.zeros((N_r_full, N_r_full))

        for n in range(1, N_r_full - 1):
            c1 = 1 / delta_r_list[n] ** 2
            c2 = 1 / 2 / delta_r_list[n]
            A[n][n - 1] = c1 - c2 / (sum(delta_r_list[0:n + 1]))
            if Func_matrix[n][m] > 0:
                A[n][n] = -2 * c1 - c3_oil - 2 / (sum(delta_r_list[0:n + 1])) ** 2 / delta_fi_list[m] ** 2
            else:
                A[n][n] = -2 * c1 - c3_water - 2 / (sum(delta_r_list[0:n + 1])) ** 2 / delta_fi_list[m] ** 2
            A[n][n + 1] = c1 + c2 / (sum(delta_r_list[0:n + 1]))


        c1 = 1 / delta_r_list[0] ** 2
        c2 = 1 / 2 / delta_r_list[0]
        if m == frac_angle_cell or m == frac_angle_2_cell:
            if Func_matrix[n][m] > 0:
                A[0][0] = -2 * c1 - c3_oil - 2 / (delta_r_list[0]) ** 2 / delta_fi_list[m] ** 2
            else:
                A[0][0] = -2 * c1 - c3_water - 2 / (delta_r_list[0]) ** 2 / delta_fi_list[m] ** 2
        else:
            if Func_matrix[n][m] > 0:
                A[0][0] = -2 * c1 - c3_oil - 2/(delta_r_list[0])**2/delta_fi_list[m]**2 + c1 - c2 / (1 * delta_r_list[0])
            else:
                A[0][0] = -2 * c1 - c3_water - 2/(delta_r_list[0])**2/delta_fi_list[m]**2 + c1 - c2 / (1 * delta_r_list[0])
        A[0][1] = c1 + c2 / (1 * delta_r_list[0])

        c1 = 1 / delta_r_list[N_r_full - 1] ** 2
        c2 = 1 / 2 / delta_r_list[N_r_full - 1]
        A[N_r_full - 1][N_r_full - 1] = -2 * c1 - c3_water + c1 + c2 / (sum(delta_r_list[0:N_r_full])) - 2/(sum(delta_r_list[0:N_r_full]))**2/delta_fi_list[m]**2
        A[N_r_full - 1][N_r_full - 2] = c1 - c2 / (sum(delta_r_list[0:N_r_full]))

        c4 = 1 / delta_fi_list[m] ** 2
        A_sym = np.zeros((N_r_full, N_r_full))
        for n in range(0, N_r_full):
            A_sym[n][n] = c4 / (sum(delta_r_list[0:n + 1])) ** 2
        A_sym_coo = coo_matrix(A_sym)

        if m == 0:
            A_line_1 = hstack([A, A_sym_coo, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym_coo])
            A_full = coo_matrix(A_line_1)
        elif m == M_fi_full-1:
            A_line_end = hstack(
                    [A_sym_coo, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym_coo, A])
            A_full = vstack([A_full, A_line_end])
        else:
            A_line = hstack([np.zeros((N_r_full, N_r_full * (m - 1))), A_sym_coo, A, A_sym_coo,
                                 np.zeros((N_r_full, N_r_full * M_fi_full - (3 + (m - 1)) * N_r_full))])
            A_full = vstack([A_full, A_line])

    j = 0
    for m in range(M_fi_full):
        for n in range(N_r_full):
            if Func_matrix[n][m] > 0:
                B[j][0] = -c3_oil * Pres_distrib[n][m]
            else:
                B[j][0] = -c3_water * Pres_distrib[n][m]
            j += 1

    def sort_func(well_coord_couple):
        return (well_coord_couple[1]) * N_r_full + well_coord_couple[0]

    wells_frac_coords.sort(key=sort_func)
    wells_frac_coords_reverse = wells_frac_coords[:: -1]

    A_full = A_full.toarray()

    A_full = coo_matrix(A_full)

    for coord_couple in wells_frac_coords_reverse:
        A_well_column_coo = A_full.getcol((coord_couple[1]-1)*N_r_full + coord_couple[0])
        A_well_column = A_well_column_coo.toarray()
        for cell_number in range(len(A_well_column)):
            if A_well_column[cell_number] != 0:
                B[cell_number] = B[cell_number] - A_well_column[cell_number]*CP_dict[coord_couple]

        A_full = A_full.tocsr()
        all_cols = np.arange(A_full.shape[1])
        cols_to_keep = np.where(np.logical_not(np.in1d(all_cols, [(coord_couple[1]-1)*N_r_full + coord_couple[0]])))[0]
        A_full = A_full[:, cols_to_keep]
        A_full = A_full[cols_to_keep, :]

        B = np.delete(B, (coord_couple[1] - 1) * N_r_full + coord_couple[0], axis=0)

    P_new = spsolve(A_full, B)

    for coord_couple in wells_frac_coords:
        P_new = np.insert(P_new, (coord_couple[1]-1)*N_r_full + coord_couple[0], CP_dict[coord_couple])

    P_new = P_new.reshape(N_r_full*M_fi_full, 1)
    Pres_end = np.zeros((N_r_full, M_fi_full))
    j = 0
    for m in range(M_fi_full):
        for n in range(N_r_full):
            Pres_end[n][m] = P_new[j][0]
            j += 1

    return Pres_end, A_full, B

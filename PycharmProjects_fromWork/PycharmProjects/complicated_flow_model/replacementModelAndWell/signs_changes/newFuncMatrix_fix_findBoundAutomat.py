import numpy as np
from matplotlib import pyplot as plt
import shapely.geometry as geom
from scipy import interpolate

def define_func_matrix(pressure_field, func_matrix, perm, mu_water, mu_oil, delta_r_list, delta_fi_list, t_step, r_well, Pinj, bound_coord, line_old, area_old, Func_coord_dict):
    M_fi_full = len(delta_fi_list)
    N_r_full = len(delta_r_list)

    X = np.zeros((N_r_full, M_fi_full))
    Y = np.zeros((N_r_full, M_fi_full))
    for m in range(M_fi_full):
        for n in range(N_r_full):
            X[n][m] = (r_well + sum(delta_r_list[0:n + 1])) * np.cos(sum(delta_fi_list[0:m]))
            Y[n][m] = (r_well + sum(delta_r_list[0:n + 1])) * np.sin(sum(delta_fi_list[0:m]))

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]

    xi = np.linspace(min(X_list), max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)

    velocity = np.zeros((N_r_full, M_fi_full))
    func_matrix_new = np.zeros((N_r_full, M_fi_full))
    for m in range(0, M_fi_full):
        for n in range(0, N_r_full):
            if n == 0 and m != 0 and m != M_fi_full-1:
                df_dr_right = (func_matrix[n+1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dr_left = df_dr_right
                df_dfi_right = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                df_dfi_left = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n]) + r_well) / \
                               delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n][m] - Pinj) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n][m] - Pinj) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            elif n == N_r_full-1 and m != 0 and m != M_fi_full-1:
                df_dr_left = (func_matrix[n][m] - func_matrix[n-1][m]) / delta_r_list[n]
                df_dr_right = df_dr_left
                df_dfi_right = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                df_dfi_left = (func_matrix[n][m] - func_matrix[n][m - 1]) / (sum(delta_r_list[0:n]) + r_well) / \
                         delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            elif m == 0 and n != 0 and n != N_r_full-1:
                df_dr_right = (func_matrix[n+1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dr_left = (func_matrix[n][m] - func_matrix[n - 1][m]) / delta_r_list[n]
                df_dfi_right = (func_matrix[n][m + 1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / \
                               delta_fi_list[m]
                df_dfi_left = df_dr_right

                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            elif m == M_fi_full-1 and n != 0 and n != N_r_full-1:
                df_dr_right = (func_matrix[n + 1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dr_left = (func_matrix[n][m] - func_matrix[n - 1][m]) / delta_r_list[n]
                df_dfi_left = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well)/delta_fi_list[m]
                df_dfi_right = df_dfi_left
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            elif n == 0 and m == 0:
                df_dr_right = (func_matrix[n + 1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dr_left = df_dr_right
                df_dfi_right = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                df_dfi_left = df_dfi_right
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - Pinj) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - Pinj) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            elif n == 0 and m == M_fi_full-1:
                df_dr_right = (func_matrix[n + 1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dr_left = df_dr_right
                df_dfi_left = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                df_dfi_right = df_dfi_left
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - Pinj) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - Pinj) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            elif n == N_r_full-1 and m == 0:
                df_dr_left = (func_matrix[n][m] - func_matrix[n-1][m]) / delta_r_list[n]
                df_dr_right = df_dr_left
                df_dfi_right = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                df_dfi_left = df_dfi_right
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            elif n == N_r_full-1 and m == M_fi_full-1:
                df_dr_left = (func_matrix[n][m] - func_matrix[n-1][m]) / delta_r_list[n]
                df_dr_right = df_dr_left
                df_dfi_left = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                df_dfi_right = df_dfi_left
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_fi_list[m] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)

            else:
                df_dr_left = (func_matrix[n][m] - func_matrix[n-1][m])/delta_r_list[n]
                df_dr_right = (func_matrix[n+1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dfi_left = (func_matrix[n][m] - func_matrix[n][m-1]) / sum(delta_r_list[0:n])/delta_fi_list[m]
                df_dfi_right = (func_matrix[n][m+1] - func_matrix[n][m]) / sum(delta_r_list[0:n]) / \
                                           delta_fi_list[m]
                df_dr = (func_matrix[n+1][m] - func_matrix[n-1][m])/delta_r_list[n]/2
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well)/delta_fi_list[m]/2

                if func_matrix[n][m] > 0:
                    V_r = -perm/mu_water*(pressure_field[n+1][m]-pressure_field[n][m])/delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m+1] - pressure_field[n][m]) / delta_fi_list[m] / (sum(delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)
            df_dfi = 0
            #velocity[n][m] = V_r*df_dr + V_fi*df_dfi
            velocity[n][m] = max(V_r, 0)*df_dr_left + min(V_r, 0)*df_dr_right + max(V_fi, 0)*df_dfi_left + min(V_fi, 0) * df_dfi_right
            #velocity[n][m] = V_r
            # if V_r < 0:
            #     print(n, m, V_r, pressure_field[n + 1][m], pressure_field[n][m])
            if velocity[n][m] < 0:
                try:
                    print(n, m, V_r, pressure_field[n + 1][m], pressure_field[n][m], df_dr_left, df_dr_right, func_matrix[n][m], func_matrix[n-1][m])
                except IndexError:
                    continue
            #if func_matrix[n][m] > 0:
            #func_matrix_new[n][m] = func_matrix[n][m]
            func_matrix_new[n][m] = func_matrix[n][m] - t_step * velocity[n][m]

            if np.isnan(func_matrix_new[n][m]):
                print(n, m)
                print(func_matrix_new[n][m], func_matrix[n][m])
                print(df_dr, df_dfi)
                print(V_r, V_fi)
    new_bound = []
    new_bound_coord = []
    Func_matrix_list = [l for l in func_matrix_new.flat]
    Func_matrix_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
    surf = plt.contour(xig, yig, Func_matrix_i, 0)
    zero_contour = surf.collections[1].get_paths()[0].vertices
    x_zero = zero_contour[:, 0]
    y_zero = zero_contour[:, 1]
    new_bound_coord = [(x_zero[i], y_zero[i]) for i in range(len(x_zero))]

    # func_matrix_new_2 = func_matrix.copy()
    # Func_matrix_list = [l for l in func_matrix_new_2.flat]
    # Func_matrix_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
    # surf = plt.contour(xig, yig, Func_matrix_i, 0)
    # zero_contour = surf.collections[1].get_paths()[0].vertices
    # x_zero = zero_contour[:, 0]
    # y_zero = zero_contour[:, 1]
    # new_bound_coord_2 = [(x_zero[i], y_zero[i]) for i in range(len(x_zero))]
    # print(new_bound_coord)
    # print(bound_coord)
    # print(new_bound_coord == bound_coord)


    line_1 = geom.LineString(new_bound_coord)
    area_1 = geom.polygon.Polygon(line_1)
    # print(line == line_old)
    # print(area == area_old)

    for i in range(len(delta_r_list)):
        r_const_row = np.zeros((1, len(delta_fi_list)))
        for j in range(len(delta_fi_list)):
            coord_x = Func_coord_dict[(i, j)][0]
            coord_y = Func_coord_dict[(i, j)][1]
            point = geom.Point(coord_x, coord_y)
            min_dist = point.distance(line_1)
            if area_1.contains(point):
                r_const_row[0][j] = -min_dist
            else:
                r_const_row[0][j] = min_dist
        if i == 0:
            Func_matrix_remake = r_const_row
        else:
            Func_matrix_remake = np.vstack((Func_matrix_remake, r_const_row))

    print('proverka')
    print(min(func_matrix.flat), max(func_matrix.flat))
    print(min(func_matrix_new.flat), max(func_matrix_new.flat))
    print(min(Func_matrix_remake.flat), max(Func_matrix_remake.flat))

    Func_matrix_list = [l for l in Func_matrix_remake.flat]
    Func_matrix_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
    surf = plt.contour(xig, yig, Func_matrix_i, 0)
    zero_contour = surf.collections[1].get_paths()[0].vertices
    x_zero = zero_contour[:, 0]
    y_zero = zero_contour[:, 1]
    bound_coord_remake = [(x_zero[i], y_zero[i]) for i in range(len(x_zero))]

    # fig = plt.figure()
    # for coord_pair in new_bound_coord:
    #     plt.scatter(coord_pair[0], coord_pair[1], color='blue')
    # for coord_pair in bound_coord:
    #     plt.scatter(coord_pair[0], coord_pair[1], color='green')
    # # for coord_pair in bound_coord_remake:
    # #     plt.scatter(coord_pair[0], coord_pair[1], color='red')
    # #
    # plt.show()


    return func_matrix_new, velocity, bound_coord_remake, Func_matrix_remake, line_1, area_1

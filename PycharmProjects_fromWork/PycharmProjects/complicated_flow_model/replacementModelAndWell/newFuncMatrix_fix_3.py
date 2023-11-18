import numpy as np
from matplotlib import pyplot as plt

def define_func_matrix(pressure_field, func_matrix, perm, mu_water, mu_oil, delta_r_list, delta_fi_list, t_step, r_well, Pinj, old_bound):
    M_fi_full = len(delta_fi_list)
    N_r_full = len(delta_r_list)
    velocity = np.zeros((N_r_full, M_fi_full))
    func_matrix_new = np.zeros((N_r_full, M_fi_full))
    for m in range(0, M_fi_full):
        for n in range(0, N_r_full):
            if n == 0 and m != 0 and m != M_fi_full-1:
                df_dr_right = (func_matrix[n+1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dr_left = df_dr_right
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m - 1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]/2
                if func_matrix[n][m] < 0:
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
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m - 1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]/2
                if func_matrix[n][m] < 0:
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
                df_dfi = (func_matrix[n][m + 1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / \
                               delta_fi_list[m]
                if func_matrix[n][m] < 0:
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
                df_dfi = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well)/delta_fi_list[m]
                if func_matrix[n][m] < 0:
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
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] < 0:
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
                df_dfi = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] < 0:
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
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] < 0:
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
                df_dfi = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] < 0:
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

                if func_matrix[n][m] < 0:
                    V_r = -perm/mu_water*(pressure_field[n+1][m]-pressure_field[n][m])/delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m+1] - pressure_field[n][m]) / delta_fi_list[m] / (sum(delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_fi_list[
                        m] / (sum(delta_r_list[0:n])+r_well)
            df_dfi = 0
            #velocity[n][m] = V_r*df_dr + V_fi*df_dfi
            velocity[n][m] = max(V_r, 0)*df_dr_left + min(V_r, 0)*df_dr_right

            # if V_r < 0:
            #     print(n, m, V_r, pressure_field[n + 1][m], pressure_field[n][m])
            if velocity[n][m] > 0 and n != N_r_full - 1 and n != 0:
                pass
                #print(n, m, V_r, pressure_field[n + 1][m], pressure_field[n][m],  pressure_field[n - 1][m], df_dr_left, df_dr_right, func_matrix[n][m], func_matrix[n-1][m], func_matrix[n+1][m])
            elif velocity[n][m] > 0 and n == N_r_full - 1:
                pass
                #print(n, m, V_r, pressure_field[n - 1][m], pressure_field[n][m], df_dr_left, df_dr_right, func_matrix[n][m], func_matrix[n-1][m])
            elif velocity[n][m] > 0 and n == 0:
                pass
                #print(n, m, V_r, pressure_field[n + 1][m], pressure_field[n][m], df_dr_left, df_dr_right,
                      #func_matrix[n][m], func_matrix[n + 1][m])

            if df_dr_right > 1:
                print(df_dr_right, n, m)

            #if func_matrix[n][m] > 0:
            func_matrix_new[n][m] = func_matrix[n][m] - t_step * velocity[n][m]


            if np.isnan(func_matrix_new[n][m]):
                print(n, m)
                print(func_matrix_new[n][m], func_matrix[n][m])
                print(df_dr, df_dfi)
                print(V_r, V_fi)
    new_bound = []
    new_bound_coord = []
    for m in range(1, M_fi_full-1):
        for n in range(1, N_r_full-1):
            if func_matrix_new[n][m] > 0:
                if (func_matrix_new[n-1][m] > 0) & (func_matrix_new[n+1][m] > 0) & (func_matrix_new[n][m-1] > 0) & (func_matrix_new[n][m+1] > 0):
                    continue
                else:
                    new_bound.append((n, m))
                    new_bound_coord.append(((r_well+sum(delta_r_list[0:n]))*np.cos(sum(delta_fi_list[0:m])), (r_well+sum(delta_r_list[0:n]))*np.sin(sum(delta_fi_list[0:m]))))
    for n in range(1, N_r_full - 1):
        if (func_matrix_new[n - 1][0] > 0) & (func_matrix_new[n + 1][0] < 0):
            new_bound.append((n, 0))
            new_bound_coord.append(((r_well + sum(delta_r_list[0:n])) * np.cos(sum(delta_fi_list[0:0])),
                                    (r_well + sum(delta_r_list[0:n])) * np.sin(sum(delta_fi_list[0:0]))))
            break
    for n in range(1, N_r_full - 1):
        if (func_matrix_new[n - 1][M_fi_full-1] > 0) & (func_matrix_new[n + 1][M_fi_full-1] < 0):
            new_bound.append((n, M_fi_full-1))
            new_bound_coord.append(((r_well + sum(delta_r_list[0:n])) * np.cos(sum(delta_fi_list[0:M_fi_full-1])),
                                    (r_well + sum(delta_r_list[0:n])) * np.sin(sum(delta_fi_list[0:M_fi_full-1]))))
            break
    print(new_bound_coord)

    # fig = plt.figure()
    # for coord_pair in new_bound_coord:
    #     plt.scatter(coord_pair[0], coord_pair[1], marker='.')
    # plt.show()
    print('hi', set(new_bound))
    print('hi2', set(old_bound))
    if set(new_bound) == set(old_bound):
        func_matrix_remake = func_matrix_new
    else:
        func_matrix_remake = np.zeros((N_r_full, M_fi_full))
        dist_to_bound = []
        for m in range(0, M_fi_full):
            for n in range(0, N_r_full):
                for coord_pair in new_bound_coord:
                    dist_to_bound.append((((r_well+sum(delta_r_list[0:n]))*np.cos(sum(delta_fi_list[0:m])) - coord_pair[0])**2 +
                                          (((r_well+sum(delta_r_list[0:n]))*np.sin(sum(delta_fi_list[0:m])) - coord_pair[1])**2))**0.5)
                min_dist = min(dist_to_bound)
                if func_matrix_new[n][m] > 0:
                    func_matrix_remake[n][m] = min_dist
                else:
                    func_matrix_remake[n][m] = -min_dist
                dist_to_bound = []






    return func_matrix_new, velocity, new_bound, new_bound_coord, func_matrix_remake

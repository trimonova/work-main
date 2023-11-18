import numpy as np
def define_func_matrix(pressure_field, func_matrix, perm, mu_water, mu_oil, delta_r_list, delta_fi_list, t_step, r_well):
    print(np.shape(pressure_field), np.shape(func_matrix))
    M_fi_full = len(delta_fi_list)
    N_r_full = len(delta_r_list)

    func_matrix_new = np.zeros((N_r_full, M_fi_full))
    for m in range(0, M_fi_full):
        for n in range(1, N_r_full):
            if n == 0 and m != 0 and m != M_fi_full-1:
                df_dr = (func_matrix[n+1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m - 1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]/2
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            elif n == N_r_full-1 and m != 0 and m != M_fi_full-1:
                df_dr = (func_matrix[n][m] - func_matrix[n-1][m]) / delta_r_list[n]
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m - 1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]/2
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            elif m == 0 and n != 0 and n != N_r_full-1:
                df_dr = (func_matrix[n+1][m] - func_matrix[n - 1][m]) / delta_r_list[n]/2
                df_dfi = (func_matrix[n][m + 1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / \
                               delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            elif m == M_fi_full-1 and n != 0 and n != N_r_full-1:
                df_dr = (func_matrix[n+1][m] - func_matrix[n - 1][m]) / delta_r_list[n]/2
                df_dfi = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well)/delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            elif n == 0 and m == 0:
                df_dr = (func_matrix[n + 1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            elif n == 0 and m == M_fi_full-1:
                df_dr = (func_matrix[n + 1][m] - func_matrix[n][m]) / delta_r_list[n]
                df_dfi = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            elif n == N_r_full-1 and m == 0:
                df_dr = (func_matrix[n][m] - func_matrix[n-1][m]) / delta_r_list[n]
                df_dfi = (func_matrix[n][m+1] - func_matrix[n][m]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            elif n == N_r_full-1 and m == M_fi_full-1:
                df_dr = (func_matrix[n][m] - func_matrix[n-1][m]) / delta_r_list[n]
                df_dfi = (func_matrix[n][m] - func_matrix[n][m-1]) / (sum(delta_r_list[0:n])+r_well) / delta_fi_list[m]
                if func_matrix[n][m] > 0:
                    V_r = -perm / mu_water * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[n]
                    V_fi = -perm / mu_water * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_r_list[n] / (sum(
                        delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n-1][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m] - pressure_field[n][m-1]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

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
                    V_fi = -perm / mu_water * (pressure_field[n][m+1] - pressure_field[n][m]) / delta_r_list[n] / (sum(delta_r_list[0:n])+r_well)
                else:
                    V_r = -perm / mu_oil * (pressure_field[n + 1][m] - pressure_field[n][m]) / delta_r_list[
                        n]
                    V_fi = -perm / mu_oil * (pressure_field[n][m + 1] - pressure_field[n][m]) / delta_r_list[
                        n] / (sum(delta_r_list[0:n])+r_well)

            #func_matrix_new[n][m] = func_matrix[n][m] - t_step*(V_r*df_dr + V_fi*df_dfi)
            func_matrix_new[n][m] = func_matrix[n][m] + t_step * 0.001

            if np.isnan(func_matrix_new[n][m]):
                print(n, m)
                print(func_matrix_new[n][m], func_matrix[n][m])
                print(df_dr, df_dr, df_dfi, df_dfi)
                print(V_r, V_fi)

    return func_matrix_new

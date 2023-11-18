from piezo_5diag_n_wells_fineMesh_frac_sparseMatrix_newCentral_BC_add_replacement_fix import PorePressure_in_Time
from piezo_5diag_n_wells_fineMesh_sparseMatrix import PorePressure_in_Time as StartPorePressure
import numpy as np
from matplotlib import pyplot as plt
from newFuncMatrix_fix import define_func_matrix
from start_to_do_replacement import replace_boundary
from scipy import interpolate
from plot_graphs import plot_graphs, plot_scatter

perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2*10**(-1)
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_oil*fi*(Cf+Cr)/(perm)
R = 0.215
r_well = 0.0075
Pres = 1*10**5
P_center = 5*10**5

frac_angle = np.pi/6
frac_angle_2 = np.pi/6*7

delta_r = 0.005
delta_r_fine = 0.001
R_for_fine = 0.015
N_r_fine = round(R_for_fine/delta_r_fine)
delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*round((R-r_well-R_for_fine)/delta_r)
N_r_full = len(delta_r_list)
delta_fi = np.pi / 60 # угол-шаг в радианах
delta_fi_fine = np.pi/180
fi_for_fine = np.pi/6
M_fi_fine = round(fi_for_fine / delta_fi_fine)
delta_fi_list_first = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi))
angle_lack = round((2*np.pi - sum(delta_fi_list_first))/delta_fi)
delta_fi_list = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi)+angle_lack)
M_fi_full = len(delta_fi_list)

x_f = 0.01

delta_t = 1

Pres_distrib = np.ones((N_r_full, M_fi_full)) * Pres
c3_oil = k_oil/delta_t
c3_water = k_water/delta_t

T_exp_1 = 10
T_exp = 10
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r**2)/100
wells_coord_real = [(0.17, np.pi/4), (0.17, np.pi/4+np.pi)]
#wells_coord_real = []
wells_angles = []
for well_coord in wells_coord_real:
    for angle_number in range(len(delta_fi_list)):
        if well_coord[1] < sum(delta_fi_list[0:angle_number]):
            wells_angles.append(angle_number-1)
            break
wells_dists = []
for well_coord in wells_coord_real:
    wells_dists.append(N_r_fine + int((well_coord[0]-r_well-R_for_fine)/delta_r))

wells_coord = list(zip(wells_dists, wells_angles))

print(wells_coord)
P_well = [2000000, 100000]

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

CP_dict_wells = CP_dict.copy()

frac_angle_cell = round((frac_angle-fi_for_fine)/delta_fi) + M_fi_fine
frac_angle_2_cell = frac_angle_cell + round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + 2*M_fi_fine
L_N_frac = 20
frac_coord_1 = [i for i in range(L_N_frac)]
frac_coord_2 = [i for i in range(L_N_frac)]
frac_pressure_1 = [5000000]*L_N_frac
frac_pressure_2 = [5000000]*L_N_frac
frac_pressure = frac_pressure_1 + frac_pressure_2
frac_coords = [(i, frac_angle_cell) for i in frac_coord_1] + [(j, frac_angle_2_cell) for j in frac_coord_2]

wells_frac_coords = wells_coord + frac_coords
print(wells_frac_coords)
for i in range(len(frac_coords)):
    CP_dict[frac_coords[i]] = frac_pressure[i]

bound_coord, Func_matrix = replace_boundary(frac_angle, frac_angle_2, r_well, delta_fi_list, delta_r_list, sum(delta_r_list[0:len(frac_coord_1)]))
plot_graphs(delta_r_list, delta_fi_list, r_well, Func_matrix)
plot_scatter(bound_coord)

if __name__ == '__main__':
    pressureInPoint = []
    for t in range(T_exp_1):
        Pres_distrib, A_full, B = StartPorePressure(N_r_full, M_fi_full, Pres_distrib, c3_water/1000, CP_dict_wells,
                                                  wells_coord, delta_r_list, delta_fi_list)
        pressureInPoint.append(Pres_distrib[wells_dists[0]][175])
    plot_graphs(delta_r_list, delta_fi_list, r_well, Pres_distrib)
    for t in range(T_exp):
        print(t)
        Pres_distrib, A_full, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict,
                                                       P_center, wells_frac_coords, Func_matrix, delta_r_list,
                                                       delta_fi_list, frac_angle_cell, frac_angle_2_cell)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))
        Func_matrix = define_func_matrix(Pres_distrib, Func_matrix, perm, mu_water, mu_oil, delta_r_list, delta_fi_list, delta_t, r_well)

        print(Func_matrix[20][10])
        pressureInPoint.append(Pres_distrib[wells_dists[0]][175])


    plt.plot(Pres_distrib[:, 50])
    plt.show()

    plt.plot(pressureInPoint)
    plt.show()

    plot_graphs(delta_r_list, delta_fi_list, r_well, Func_matrix)
    plot_graphs(delta_r_list, delta_fi_list, r_well, Pres_distrib)

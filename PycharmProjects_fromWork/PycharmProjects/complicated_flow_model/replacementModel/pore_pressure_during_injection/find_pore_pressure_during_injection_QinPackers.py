from input_parameters import N_r_full, M_fi_full, T_exp_dir, c3_oil, c3_water, CP_dict, q_coef
from input_parameters import wells_coord, delta_r_list, delta_fi_list, porosity, C_total, perm, \
    delta_t, r_well, frac_angle, frac_angle_2, frac_length_1, frac_length_2, mu_oil, mu_water, \
    coord_matrix_rad, coord_matrix_cart, q, M_1, M_2, N_1, N_2, wells_frac_coords, delta_r, delta_fi, X_matrix, Y_matrix

from QinPackers.newFuncMatrix_fix_5 import define_func_matrix
import copy
from QinPackers.find_viscosity import find_viscosity

from QinPackers.find_bound_coords import find_bound_coords
from QinPackers.find_func_matrix_remake import find_func_matrix_remake
from QinPackers.start_to_do_replacement import replace_boundary
from QinPackers.find_pore_pressure_QinPackers import PorePressure_in_Time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Pres_distrib_in_time = np.load('../pore_pressure_before_injection/Pres_distrib_before_fracturing.npy')
Pres_distrib = Pres_distrib_in_time[-1]

T_exp = 60


Func_matrix_remake = replace_boundary(frac_angle, frac_angle_2, 0.004, frac_length_1, frac_length_2, M_fi_full, N_r_full, coord_matrix_cart)

viscosity_matrix = find_viscosity(mu_oil, mu_water, Func_matrix_remake, coord_matrix_rad, delta_r_list,
                                 delta_fi_list, N_r_full, M_fi_full)
print(max(viscosity_matrix.flat), min(viscosity_matrix.flat))

Pres_distrib_in_Time_2 = []
Func_matrix_in_Time = []
velocity_in_Time = []
bound_coords_rad_in_Time = []
bound_coords_cart_in_Time = []
Func_matrix_remake_in_Time = []
viscosity_in_Time = []

for t in range(T_exp):
    print(t)
    Pres_distrib, A, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, q, wells_frac_coords,
                                              wells_coord, delta_r_list, delta_fi_list, viscosity_matrix, porosity,
                                              C_total, perm, delta_t,
                                              r_well, M_1, M_2, N_1, N_2)
    Func_matrix, velocity = define_func_matrix(Pres_distrib, Func_matrix_remake, perm, delta_r_list, delta_fi_list,
                                               delta_t, r_well, q, viscosity_matrix,  N_r_full, M_fi_full)
    bound_coords_rad_new, bound_coords_cart_new = find_bound_coords(Func_matrix, coord_matrix_rad, delta_r_list,
                                                                    delta_fi_list, M_fi_full, N_r_full)

    Func_matrix_remake = find_func_matrix_remake(coord_matrix_cart, M_fi_full, N_r_full, bound_coords_cart_new)

    viscosity_matrix = find_viscosity(mu_oil, mu_water, Func_matrix_remake, coord_matrix_rad, delta_r_list,
                                      delta_fi_list, N_r_full, M_fi_full)

    Pres_distrib_in_Time_2.append(Pres_distrib)
    Func_matrix_in_Time.append(Func_matrix)
    velocity_in_Time.append(velocity)
    bound_coords_rad_in_Time.append(bound_coords_rad_new)
    bound_coords_cart_in_Time.append(bound_coords_cart_new)
    Func_matrix_remake_in_Time.append(Func_matrix_remake)
    viscosity_in_Time.append(viscosity_matrix)

np.save('result_folder/Pres_distrib_in_Time_2', Pres_distrib_in_Time_2)
np.save('result_folder/viscosity_in_Time', viscosity_in_Time)
np.save('result_folder/Func_matrix_in_Time', Func_matrix_in_Time)
np.save('result_folder/Func_matrix_remake_in_Time', Func_matrix_remake_in_Time)
np.save('result_folder/velocity_in_Time', velocity_in_Time)
np.save('result_folder/bound_coords_cart_in_Time', bound_coords_cart_in_Time, allow_pickle=True)
np.save('result_folder/bound_coords_rad_in_Time', bound_coords_rad_in_Time, allow_pickle=True)




surf = plt.contourf(X_matrix, Y_matrix, Pres_distrib, linewidth=0.2, cmap=plt.get_cmap('jet'))
plt.title('pressure')
plt.colorbar(surf)
plt.show()

surf = plt.contourf(X_matrix, Y_matrix, viscosity_matrix, linewidth=0.2, cmap=plt.get_cmap('jet'))
plt.title('viscosity')
plt.colorbar(surf)
plt.show()

surf = plt.contourf(X_matrix, Y_matrix, Func_matrix_remake, linewidth=0.2, cmap=plt.get_cmap('jet'))
plt.title('Func_matrix')
plt.colorbar(surf)
plt.show()

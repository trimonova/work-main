import numpy as np
from def_find_pore_pressure_QinCenter_equal_0 import PorePressure_in_Time as PorePressure_in_Time_1
from input_parameters import N_r_full, M_fi_full, T_exp_dir, Pres_distrib, c3_oil, c3_water, CP_dict, q_coef
from input_parameters import wells_coord, delta_r_list, delta_fi_list, viscosity_matrix, porosity, C_total, perm, delta_t, r_well

Pres_distrib_in_time = []
for t in range(T_exp_dir):
    print(t)
    Pres_distrib, A, B = PorePressure_in_Time_1(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, q_coef,
                                              wells_coord, delta_r_list, delta_fi_list, viscosity_matrix, porosity,
                                              C_total, perm, delta_t,
                                              r_well)
    Pres_distrib_in_time.append(Pres_distrib)


np.save('Pres_distrib_before_fracturing', Pres_distrib_in_time)
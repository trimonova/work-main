import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib import cm
import copy
from read_mat_file import points_coords_dict_rad as points_r_fi_pair
from read_mat_file import P_filter_0, P_filter_1, P_filter_2, P_filter_3, P_filter_4, P_filter_5
from read_mat_file import P_filter_6, P_filter_7, P_filter_8, P_filter_9, P_filter_10, P_filter_11, P_filter_12, P_filter_13
#print(P_exp)
from QinPackers.find_pore_pressure_QinCenter_equal_0 import PorePressure_in_Time as PorePressure_in_Time_1
from QinPackers.find_pore_pressure_QinPackers import PorePressure_in_Time as PorePressure_in_Time_2
import pandas as pd

perm = 2 * 10 ** (-12)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2 * 10 ** (-3)
porosity = 0.4  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
C_total = (Cf + Cr) * 25
k_water = mu_water * porosity * C_total / perm
k_oil = mu_oil * porosity * C_total / perm
frac_angle = np.pi/4
frac_angle_2 = np.pi/4*5
frac_length_1 = 0.01
frac_length_2 = 0.01
delta_r = 0.0001
delta_r_fine = 0.0001
R_for_fine = 0.02
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine / delta_r_fine)
delta_r_list = [delta_r_fine] * N_r_fine + [delta_r] * round((R - r_well - R_for_fine) / delta_r)
N_r_full = len(delta_r_list)

delta_fi = np.pi / 45  # угол-шаг в радианах
delta_fi_fine = np.pi / 45
fi_for_fine = np.pi / 6
M_fi_fine = round(fi_for_fine / delta_fi_fine)

delta_fi_list_first = [delta_fi] * round((frac_angle - fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine * 2) + [
    delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine * 2) + [
                          delta_fi] * (round((2 * np.pi - frac_angle_2 - fi_for_fine) / delta_fi))
angle_lack = round((2 * np.pi - sum(delta_fi_list_first)) / delta_fi)
# delta_fi_list = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi)+angle_lack)
delta_fi_list = [delta_fi] * round(2 * np.pi / delta_fi)
M_fi_full = len(delta_fi_list)

coord_matrix_rad = []
coord_matrix_cart = []
for n in range(len(delta_r_list)):
    coord_line_rad = []
    coord_line_cart = []
    r = sum(delta_r_list[0:n]) + r_well
    for m in range(len(delta_fi_list)):
        fi = sum(delta_fi_list[0:m])
        coord_line_rad.append((r, fi))
        coord_line_cart.append((r * np.cos(fi), r * np.sin(fi)))
        # coord_matrix_rad[n][m] = (r, fi)
        # coord_matrix_cart[n][m] = (r*np.cos(fi), r*np.sin(fi))
    coord_matrix_rad.append(coord_line_rad)
    coord_matrix_cart.append(coord_line_cart)

delta_t = 0.5
Pres = 1 * 10 ** 5
# P_center = 60*10**5
#Q_center = 0.2 * 10 ** (-8)  # из лаб. данных - 0.2*10**(-6) m3/s
Q_center = 0
s = 4 * 10 ** (-5)  # ширина прорези 2 мм, высота - 1 см
q = Q_center / s  # m/s
q_coef = q * delta_r_fine * mu_oil / perm

Pres_distrib = np.ones((N_r_full, M_fi_full)) * Pres
c3_oil = k_oil / delta_t
c3_water = k_water / delta_t
T_exp_dir = 2
Courant_number_oil = (delta_t / k_oil / delta_fi ** 2 + delta_t / k_oil / delta_r_fine ** 2) / 25
Courant_number_water = (delta_t / k_water / delta_fi ** 2 + delta_t / k_water / delta_r_fine ** 2) / 25
print(Courant_number_water, Courant_number_oil)
# wells_coord_real = [(0.17, np.pi/4), (0.17, np.pi/4 + np.pi)]

wells_coord = [(round(0.1711/delta_r), round(np.pi/4/delta_fi)), (round(0.1711/delta_r), round((np.pi+np.pi/4)/delta_fi))]
P_well = [1000000, 0]
#P_well = []
#wells_coord = []
print(wells_coord)

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]


def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

points_coords = {}
for key in points_r_fi_pair:
    points_coords[key] = ((round(points_r_fi_pair[key][0]/delta_r)), (round(points_r_fi_pair[key][1]/delta_fi)))
#points_coords = [(round(i[0]/delta_r), round(i[1]/delta_fi)) for i in points_r_fi_pair]
print('points_coords', points_coords)

one_point0 = points_coords[0]
one_point1 = points_coords[1]
one_point2 = points_coords[2]
one_point3 = points_coords[3]
one_point4 = points_coords[4]
one_point5 = points_coords[5]
one_point6 = points_coords[6]
one_point7 = points_coords[7]
one_point8 = points_coords[8]
one_point9 = points_coords[9]
one_point10 = points_coords[10]
one_point11 = points_coords[11]
one_point12 = points_coords[12]
one_point13 = points_coords[13]

PinPoint0 = []
PinPoint1 = []
PinPoint2 = []
PinPoint3 = []
PinPoint4 = []
PinPoint5 = []
PinPoint6 = []
PinPoint7 = []
PinPoint8 = []
PinPoint9 = []
PinPoint10 = []
PinPoint11 = []
PinPoint12 = []
PinPoint13 = []
PinPointCenter = []

viscosity_matrix = np.ones((N_r_full, M_fi_full))*mu_water

# for t in range(T_exp_dir):
#     print(t)
#     Pres_distrib, A, B = PorePressure_in_Time_1(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, q_coef,
#                                               wells_coord, delta_r_list, delta_fi_list, viscosity_matrix, porosity,
#                                               C_total, perm, delta_t,
#                                               r_well)
#     print(min(Pres_distrib.flat), max(Pres_distrib.flat))
#     print('center', Pres_distrib[1][11])
#     print(np.shape(Pres_distrib))
#     PinPoint13.append(Pres_distrib[int(one_point13[0])][int(one_point13[1])])
#     print(Pres_distrib[int(one_point13[0])][int(one_point13[1])])
#     print(int(one_point13[0]), int(one_point13[1]))
#     print(Pres_distrib[0][10])
#     print(Pres_distrib[1][10])
#     print(Pres_distrib[0][55])
#     print(Pres_distrib[1][55])
#     PinPoint2.append(Pres_distrib[int(one_point2[0])][int(one_point2[1])])
#     PinPoint3.append(Pres_distrib[int(one_point3[0])][int(one_point3[1])])
#     PinPoint4.append(Pres_distrib[int(one_point4[0])][int(one_point4[1])])
#     PinPoint5.append(Pres_distrib[int(one_point5[0])][int(one_point5[1])])
#     PinPoint6.append(Pres_distrib[int(one_point6[0])][int(one_point6[1])])
#     PinPoint7.append(Pres_distrib[int(one_point7[0])][int(one_point7[1])])
#     PinPoint8.append(Pres_distrib[int(one_point8[0])][int(one_point8[1])])
#     PinPoint9.append(Pres_distrib[int(one_point9[0])][int(one_point9[1])])
#     PinPoint10.append(Pres_distrib[int(one_point10[0])][int(one_point10[1])])
#     PinPoint11.append(Pres_distrib[int(one_point11[0])][int(one_point11[1])])
#     PinPoint12.append(Pres_distrib[int(one_point12[0])][int(one_point12[1])])
#     PinPoint0.append(Pres_distrib[int(one_point0[0])][int(one_point0[1])])
#     PinPoint1.append(Pres_distrib[int(one_point1[0])][int(one_point1[1])])
#     PinPointCenter.append(Pres_distrib[int(one_point0[0])][int(one_point0[1])]) # в центр записываем давление, равное давление на скважине

#np.save('Pres_distrib_before_fracturing', Pres_distrib)
Pres_distrib = np.load('Pres_distrib_before_fracturing.npy')

t_list = [j*delta_t for j in range(T_exp_dir)]
t_list_exp = [i*0.01 for i in range(np.shape(P_filter_1)[0])]

# T_add = 150
# for t in range(int(T_add/delta_t)):
#     PinPoint2.append(PinPoint2[-1])
#     PinPoint3.append(PinPoint3[-1])
#     PinPoint4.append(PinPoint4[-1])
#     PinPoint5.append(PinPoint5[-1])
#     PinPoint6.append(PinPoint6[-1])
#     PinPoint7.append(PinPoint7[-1])
#     PinPoint8.append(PinPoint8[-1])
#     PinPoint9.append(PinPoint9[-1])
#     PinPoint10.append(PinPoint10[-1])
#     PinPoint11.append(PinPoint11[-1])
#     PinPoint12.append(PinPoint12[-1])
#     PinPoint0.append(PinPoint0[-1])
#     PinPoint1.append(PinPoint1[-1])
#     PinPoint13.append(PinPoint13[-1])

t_list = [j*delta_t for j in range(np.shape(PinPoint1)[0])]
t_list_exp = [i*0.01 for i in range(np.shape(P_filter_1)[0])]

fig = plt.figure()
plt.plot(t_list_exp, P_filter_1)
plt.plot(t_list, np.array(PinPoint1)/10**6)
plt.show()

plt.plot(t_list_exp, P_filter_2)
plt.plot(t_list, np.array(PinPoint2)/10**6)
plt.show()

plt.plot(t_list_exp, P_filter_3)
plt.plot(t_list, np.array(PinPoint3)/10**6)
plt.show()

plt.plot(t_list_exp, P_filter_4)
plt.plot(t_list, np.array(PinPoint4)/10**6)
plt.show()

plt.plot(t_list_exp, P_filter_5)
plt.plot(t_list, np.array(PinPoint13)/10**6)
plt.show()

plt.plot(t_list_exp, P_filter_10)
plt.plot(t_list, np.array(PinPoint10)/10**6)
plt.show()

X = np.zeros((N_r_full,M_fi_full))
Y = np.zeros((N_r_full, M_fi_full))
for m in range(M_fi_full):
    for n in range(N_r_full):
        X[n][m] = (r_well+(n+1)*delta_r)*np.cos(delta_fi*m)
        Y[n][m] = (r_well+(n+1)*delta_r)*np.sin(delta_fi*m)

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
P_list = [k for k in Pres_distrib.flat]


CP_list = zip(X_list, Y_list, P_list)

print(min(P_list), max(P_list))

xi = np.linspace(min(X_list),max(X_list), 700)
yi = np.linspace(min(Y_list), max(Y_list), 700)
xig, yig = np.meshgrid(xi, yi)
Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')
levels = list(range(0,2500000,10000))
fig = plt.figure()
surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()




from QinPackers.newFuncMatrix_fix_5 import define_func_matrix
import copy
from QinPackers.find_viscosity import find_viscosity

from QinPackers.find_bound_coords import find_bound_coords
from QinPackers.find_func_matrix_remake import find_func_matrix_remake
from QinPackers.start_to_do_replacement import replace_boundary
#from QinCenter.find_pore_pressure_QinCenter import PorePressure_in_Time
import pandas as pd

angle = 0
M_1 = 0
M_2 = 0
for i in range(len(delta_fi_list) - 1):
    angle = angle + delta_fi_list[i]
    if angle <= frac_angle <= (angle + delta_fi_list[i + 1]):
        print(frac_angle, angle, angle + delta_fi_list[i + 1])
        if abs(frac_angle - angle) < abs(frac_angle - (angle + delta_fi_list[i + 1])):
            M_1 = copy.copy(i)
        else:
            M_1 = copy.copy(i + 1)
    if angle <= frac_angle_2 <= (angle + delta_fi_list[i + 1]):
        print(frac_angle_2, angle, angle + delta_fi_list[i + 1])
        if abs(frac_angle - angle) < abs(frac_angle - (angle + delta_fi_list[i + 1])):
            M_2 = copy.copy(i)
        else:
            M_2 = copy.copy(i + 1)

print(M_1, M_2)
if M_1 == 0 or M_2 == 0:
    print("not found M_1 or M_2")
    raise ValueError

R_i = 0
N_1 = 0
N_2 = 0
for i in range(len(delta_r_list) - 1):
    R_i = R_i + delta_r_list[i]
    if R_i <= frac_length_1 <= (R_i + delta_r_list[i + 1]):
        print(frac_length_1, R_i, R_i + delta_r_list[i + 1])
        if abs(frac_length_1 - R_i) < abs(frac_length_1 - (R_i + delta_r_list[i + 1])):
            N_1 = copy.copy(i)
        else:
            N_1 = copy.copy(i + 1)
    if R_i <= frac_length_2 <= (R_i + delta_r_list[i + 1]):
        print(frac_length_2, R_i, R_i + delta_r_list[i + 1])
        if abs(frac_length_2 - R_i) < abs(frac_length_2 - (R_i + delta_r_list[i + 1])):
            N_2 = copy.copy(i)
        else:
            N_2 = copy.copy(i + 1)

print(N_1, N_2)
if N_1 == 0 or N_2 == 0:
    print("not found N_1 or N_2")
    raise ValueError

frac_coords = []

for i in range(N_1):
    frac_coords.append((copy.copy(i), M_1))
for i in range(N_2):
    frac_coords.append((copy.copy(i), M_2))

for elem in frac_coords:
    #CP_dict[elem] = 10**5
    CP_dict[elem] = 1

wells_frac_coords = wells_coord + frac_coords

mu_oil = 2 * 10 ** (-1)
perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
porosity = 0.4  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
C_total = (Cf + Cr) * 25
k_water = mu_water * porosity * C_total / perm
k_oil = mu_oil * porosity * C_total / perm
c3_oil = k_oil / delta_t
c3_water = k_water / delta_t

# P_center = 60*10**5
Q_center = 0.2 * 10 ** (-8)  # из лаб. данных - 0.2*10**(-6) m3/s
s = 4 * 10 ** (-5)  # ширина прорези 2 мм, высота - 1 см
q = Q_center / s  # m/s
q_coef = q * delta_r_fine * mu_oil / perm

T_exp = 2



Func_matrix_remake = replace_boundary(frac_angle, frac_angle_2, 0.004, frac_length_1, frac_length_2, M_fi_full, N_r_full, coord_matrix_cart)

print(N_r_full, M_fi_full)
viscosity_matrix = find_viscosity(mu_oil, mu_water, Func_matrix_remake, coord_matrix_rad, delta_r_list,
                                 delta_fi_list, N_r_full, M_fi_full)
print(max(viscosity_matrix.flat), min(viscosity_matrix.flat))


velocity_in_point_1 = []
velocity_in_point_2 = []
velocity_in_point_3 = []
pressure_in_point_1 = []
pressure_in_point_2 = []
pressure_in_point_3 = []
dataFrame_to_study = pd.DataFrame()

for t in range(T_exp):
    print(t)
    Pres_distrib, A, B = PorePressure_in_Time_2(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, q, wells_frac_coords,
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

    X = np.zeros((N_r_full, M_fi_full))
    Y = np.zeros((N_r_full, M_fi_full))
    for m in range(M_fi_full):
        for n in range(N_r_full):
            X[n][m] = (r_well + (n + 1) * delta_r) * np.cos(delta_fi * m)
            Y[n][m] = (r_well + (n + 1) * delta_r) * np.sin(delta_fi * m)

    surf = plt.contourf(X, Y, Pres_distrib, linewidth=0.2, cmap=plt.get_cmap('jet'))
    plt.title('pressure')
    plt.colorbar(surf)
    plt.show()

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]
    P_list = [k for k in Pres_distrib.flat]

    CP_list = zip(X_list, Y_list, P_list)

    print(min(P_list), max(P_list))

    xi = np.linspace(min(X_list), max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')
    levels = list(range(0, 2500000, 10000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),
                        linewidth=0.2, levels=levels)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    PinPoint13.append(Pres_distrib[int(one_point13[0])][int(one_point13[1])])
    PinPoint2.append(Pres_distrib[int(one_point2[0])][int(one_point2[1])])
    PinPoint3.append(Pres_distrib[int(one_point3[0])][int(one_point3[1])])
    PinPoint4.append(Pres_distrib[int(one_point4[0])][int(one_point4[1])])
    PinPoint5.append(Pres_distrib[int(one_point5[0])][int(one_point5[1])])
    PinPoint6.append(Pres_distrib[int(one_point6[0])][int(one_point6[1])])
    PinPoint7.append(Pres_distrib[int(one_point7[0])][int(one_point7[1])])
    PinPoint8.append(Pres_distrib[int(one_point8[0])][int(one_point8[1])])
    PinPoint9.append(Pres_distrib[int(one_point9[0])][int(one_point9[1])])
    PinPoint10.append(Pres_distrib[int(one_point10[0])][int(one_point10[1])])
    PinPoint11.append(Pres_distrib[int(one_point11[0])][int(one_point11[1])])
    PinPoint12.append(Pres_distrib[int(one_point12[0])][int(one_point12[1])])
    PinPoint0.append(Pres_distrib[int(one_point0[0])][int(one_point0[1])])
    PinPoint1.append(Pres_distrib[int(one_point1[0])][int(one_point1[1])])
    PinPointCenter.append(Pres_distrib[int(one_point10[0])][int(one_point10[1])] + q_coef)

    dataFrame_to_study['vel ' + str(t)] = velocity[0:N_r_full, int(M_fi_full / 8 * 3)]
    dataFrame_to_study['Pres ' + str(t)] = Pres_distrib[0:N_r_full, int(M_fi_full / 8 * 3)]
    dataFrame_to_study['FMat ' + str(t)] = Func_matrix[0:N_r_full, int(M_fi_full / 8 * 3)]
    dataFrame_to_study['FMatRem ' + str(t)] = Func_matrix_remake[0:N_r_full, int(M_fi_full / 8 * 3)]
    dataFrame_to_study['Visc ' + str(t)] = viscosity_matrix[0:N_r_full, int(M_fi_full / 8 * 3)]

    if t == 9 or t == 10 or t == 19 or t == 20 or t == 49 or t == 50 or t == 89 or t == 90 or t == 91:
        print('print_dataframe')
        dataFrame_vel = pd.DataFrame(velocity)
        dataFrame_visc = pd.DataFrame(viscosity_matrix)
        dataFrame_func_matrix = pd.DataFrame(Func_matrix)
        dataFrame_func_matrix_remake = pd.DataFrame(Func_matrix_remake)
        dataFrame_Pres_distrib = pd.DataFrame(Pres_distrib)

        dataFrame_vel.to_csv(
            'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\\dataframes_QinCenter\\dataframe_vel_t_' + str(
                t) + '.csv')
        dataFrame_visc.to_csv(
            'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\\dataframes_QinCenter\\dataframe_visc_t_' + str(
                t) + '.csv')
        dataFrame_func_matrix.to_csv(
            'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\\dataframes_QinCenter\\dataframe_func_matr_t_' + str(
                t) + '.csv')
        dataFrame_func_matrix_remake.to_csv(
            'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\\dataframes_QinCenter\\dataframe_func_matr_rem_t_' + str(
                t) + '.csv')
        dataFrame_Pres_distrib.to_csv(
            'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\\dataframes_QinCenter\\dataframe_P_res_t_' + str(
                t) + '.csv')

    # old_bound = copy.deepcopy(new_bound)
    velocity_in_point_1.append(velocity[1, int(M_fi_full / 8 * 3)])
    velocity_in_point_2.append(velocity[8, int(M_fi_full / 8 * 3)])
    velocity_in_point_3.append(velocity[9, int(M_fi_full / 8 * 3)])
    pressure_in_point_1.append(Pres_distrib[1, int(M_fi_full / 8 * 3)])
    pressure_in_point_2.append(Pres_distrib[8, int(M_fi_full / 8 * 3)])
    pressure_in_point_3.append(Pres_distrib[9, int(M_fi_full / 8 * 3)])
    print(min(Func_matrix.flat), max(Func_matrix.flat))
    print(min(Func_matrix_remake.flat), max(Func_matrix_remake.flat))
    print(min(velocity.flat), max(velocity.flat))
    print(velocity[:, int(M_fi_full / 8 * 3)])
    pres = Pres_distrib[:, 10]
    vel = velocity[:, 10]
    Func_matr = Func_matrix[:, 10]
    Func_matr_rem = Func_matrix_remake[:, 10]
    visc = viscosity_matrix[:, 10]
    print(min(viscosity_matrix.flat), max(viscosity_matrix.flat))

dataFrame_to_study.to_csv(
    'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\\dataframes_QinCenter\\dataframe.csv')
plt.figure()
plt.plot(velocity[:, int(M_fi_full / 8 * 3)])
plt.show()


# t_list = [9000 + i*delta_t_dir for i in range(T_exp_dir)] + [9000 + T_exp_dir*delta_t_dir + j*delta_t_neu for j in range(T_exp_neu)]
t_list = [delta_t*i for i in range(T_exp_dir)]+ [T_exp_dir*delta_t + j*delta_t for j in range(T_exp)]
t_list_exp = [i*0.01 for i in range(np.shape(P_filter_1)[0])]


fig, axs = plt.subplots(3, 2, figsize=(21, 10))
axs[0, 0].set_title(str((round(points_r_fi_pair[1][0], 2), round(points_r_fi_pair[1][1], 2))), y = 0.75, loc='left')
axs[0, 0].set_xlabel("Time, s")
axs[0, 0].set_ylabel("Pressure, MPa")
axs[0, 0].plot(t_list_exp, P_filter_1)
axs[0, 0].plot(t_list, np.array(PinPoint1)/10**6)

axs[0, 1].set_title(str((round(points_r_fi_pair[2][0], 2), round(points_r_fi_pair[2][1], 2))), y = 0.75, loc='left')
axs[0, 1].set_xlabel("Time, s")
axs[0, 1].set_ylabel("Pressure, MPa")
axs[0, 1].plot(t_list_exp, P_filter_2)
axs[0, 1].plot(t_list, np.array(PinPoint2)/10**6)

axs[2, 1].set_title(str((round(points_r_fi_pair[10][0], 2), round(points_r_fi_pair[10][1], 2))), y = 0.75, loc='left')
axs[2, 1].plot(t_list_exp, P_filter_10)
axs[2, 1].plot(t_list, np.array(PinPoint10)/10**6)
axs[2, 1].plot(t_list, np.array(PinPointCenter)/10**6)

axs[1, 1].set_title(str((round(points_r_fi_pair[6][0], 2), round(points_r_fi_pair[6][1], 2))), y = 0.75, loc='left')
axs[1, 1].set_xlabel("Time, s")
axs[1, 1].set_ylabel("Pressure, MPa")
axs[1, 1].plot(t_list_exp, P_filter_6)
axs[1, 1].plot(t_list, np.array(PinPoint6)/10**6)

axs[1, 0].set_title(str((round(points_r_fi_pair[7][0], 2), round(points_r_fi_pair[7][1], 2))), y = 0.75, loc='left')
axs[1, 0].set_xlabel("Time, s")
axs[1, 0].set_ylabel("Pressure, MPa")
axs[1, 0].plot(t_list_exp, P_filter_7)
axs[1, 0].plot(t_list, np.array(PinPoint7)/10**6)

axs[2, 0].set_title(str((round(points_r_fi_pair[9][0], 2), round(points_r_fi_pair[9][1], 2))), y = 0.75, loc='left')
axs[2, 0].set_xlabel("Time, s")
axs[2, 0].set_ylabel("Pressure, MPa")
axs[2, 0].plot(t_list_exp, P_filter_9)
axs[2, 0].plot(t_list, np.array(PinPoint9)/10**6)

plt.show()

fig, axs = plt.subplots(3, 2, figsize=(21, 10))
axs[0, 0].set_title(str((round(points_r_fi_pair[0][0], 2), round(points_r_fi_pair[0][1], 2))), y = 0.75, loc='left')
axs[0, 0].set_xlabel("Time, s")
axs[0, 0].set_ylabel("Pressure, MPa")
axs[0, 0].plot(t_list_exp, P_filter_0)
axs[0, 0].plot(t_list, np.array(PinPoint0)/10**6)

axs[0, 1].set_title(str((round(points_r_fi_pair[3][0], 2), round(points_r_fi_pair[3][1], 2))), y = 0.75, loc='left')
axs[0, 1].set_xlabel("Time, s")
axs[0, 1].set_ylabel("Pressure, MPa")
axs[0, 1].plot(t_list_exp, P_filter_3)
axs[0, 1].plot(t_list, np.array(PinPoint3)/10**6)

axs[2, 1].set_title(str((round(points_r_fi_pair[4][0], 2), round(points_r_fi_pair[4][1], 2))), y = 0.75, loc='left')
axs[2, 1].plot(t_list_exp, P_filter_4)
axs[2, 1].plot(t_list, np.array(PinPoint4)/10**6)

axs[1, 1].set_title(str((round(points_r_fi_pair[5][0], 2), round(points_r_fi_pair[6][1], 2))), y = 0.75, loc='left')
axs[1, 1].set_xlabel("Time, s")
axs[1, 1].set_ylabel("Pressure, MPa")
axs[1, 1].plot(t_list_exp, P_filter_5)
axs[1, 1].plot(t_list, np.array(PinPoint5)/10**6)

axs[1, 0].set_title(str((round(points_r_fi_pair[8][0], 2), round(points_r_fi_pair[7][1], 2))), y = 0.75, loc='left')
axs[1, 0].set_xlabel("Time, s")
axs[1, 0].set_ylabel("Pressure, MPa")
axs[1, 0].plot(t_list_exp, P_filter_8)
axs[1, 0].plot(t_list, np.array(PinPoint8)/10**6)

axs[2, 0].set_title(str((round(points_r_fi_pair[11][0], 2), round(points_r_fi_pair[11][1], 2))), y = 0.75, loc='left')
axs[2, 0].set_xlabel("Time, s")
axs[2, 0].set_ylabel("Pressure, MPa")
axs[2, 0].plot(t_list_exp, P_filter_11)
axs[2, 0].plot(t_list, np.array(PinPoint11)/10**6)

plt.show()


X = np.zeros((N_r_full,M_fi_full))
Y = np.zeros((N_r_full, M_fi_full))
for m in range(M_fi_full):
    for n in range(N_r_full):
        X[n][m] = (r_well+(n+1)*delta_r)*np.cos(delta_fi*m)
        Y[n][m] = (r_well+(n+1)*delta_r)*np.sin(delta_fi*m)

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
P_list = [k for k in Pres_distrib.flat]


CP_list = zip(X_list, Y_list, P_list)

print(min(P_list), max(P_list))

xi = np.linspace(min(X_list),max(X_list), 700)
yi = np.linspace(min(Y_list), max(Y_list), 700)
xig, yig = np.meshgrid(xi, yi)
Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')
levels = list(range(0,1500000,10000))
fig = plt.figure()
surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

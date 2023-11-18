# Решение уравнения пьезопроводности в цилиндрических координатах, двумерный случай (неявная схема)
# В центре скважина с постоянным давлением P, на границах задается: градиент давления равен нулю.
# Учитывается накопление перемещений в точках, если значение в new_func_matrix  меньше шага по сетке
# ... чтобы можно было шаг по врмение уменьшать
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve
from newFuncMatrix_fix_4 import define_func_matrix
import copy
from find_viscosity import find_viscosity

from find_bound_coords import find_bound_coords
from find_func_matrix_remake import find_func_matrix_remake
from find_pore_pressure_QinCenter import PorePressure_in_Time

perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2*10**(-1)
porosity = 0.4  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
C_total = (Cf + Cr)*25
k_water = mu_water*porosity*C_total/perm
k_oil = mu_oil*porosity*C_total/perm
frac_angle = np.pi/6
frac_angle_2 = np.pi/6*7
delta_r = 0.0001
delta_r_fine = 0.0001
R_for_fine = 0.02
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine/delta_r_fine)
delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*round((R-r_well-R_for_fine)/delta_r)
N_r_full = len(delta_r_list)

delta_fi = np.pi / 18 # угол-шаг в радианах
delta_fi_fine = np.pi/18
fi_for_fine = np.pi/6
M_fi_fine = round(fi_for_fine / delta_fi_fine)

delta_fi_list_first = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi))
angle_lack = round((2*np.pi - sum(delta_fi_list_first))/delta_fi)
#delta_fi_list = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi)+angle_lack)
delta_fi_list = [delta_fi]*round(2*np.pi/delta_fi)
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
        coord_line_cart.append((r*np.cos(fi), r*np.sin(fi)))
        # coord_matrix_rad[n][m] = (r, fi)
        # coord_matrix_cart[n][m] = (r*np.cos(fi), r*np.sin(fi))
    coord_matrix_rad.append(coord_line_rad)
    coord_matrix_cart.append(coord_line_cart)




delta_t = 0.5
Pres = 1*10**5
#P_center = 60*10**5
Q_center = 0.2*10**(-8) # из лаб. данных - 0.2*10**(-6) m3/s
s = 4*10**(-5) # ширина прорези 2 мм, высота - 1 см
q = Q_center/s # m/s
q_coef = q*delta_r_fine*mu_oil/perm

Pres_distrib = np.ones((N_r_full, M_fi_full)) * Pres
c3_oil = k_oil/delta_t
c3_water = k_water/delta_t
T_exp = 100
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r_fine**2)/25
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r_fine**2)/25
print(Courant_number_water, Courant_number_oil)
# wells_coord_real = [(0.17, np.pi/4), (0.17, np.pi/4 + np.pi)]
wells_coord_real = []
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

# P_well = [1000000, 100000]
P_well = []

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]


def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

r_bound_init = 0.005
bound_coords = []
for fi in range(360):
    bound_coords.append(((r_well+r_bound_init)*np.cos(fi*np.pi/180), (r_well + r_bound_init) * np.sin(fi * np.pi / 180)))

Func_matrix = find_func_matrix_remake(coord_matrix_cart, M_fi_full, N_r_full, bound_coords)

for coord_pair in set(bound_coords):
    plt.scatter(coord_pair[0], coord_pair[1])
plt.show()

Func_matrix_remake = Func_matrix

print(min(Func_matrix.flat), max(Func_matrix.flat))
print(N_r_full, M_fi_full)
viscosity_matrix = find_viscosity(mu_oil, mu_water, Func_matrix_remake, coord_matrix_rad, delta_r_list,
                                 delta_fi_list, N_r_full, M_fi_full)
print(viscosity_matrix[49][4])
Func_matr_rem_1 = Func_matrix_remake[:, 4]
Func_matr_rem_2 = Func_matrix_remake[:, 5]
visc_1 = viscosity_matrix[:, 4]
visc_2 = viscosity_matrix[:, 5]

X = np.zeros((N_r_full,M_fi_full))
Y = np.zeros((N_r_full, M_fi_full))
for m in range(M_fi_full):
    for n in range(N_r_full):
        X[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.cos(sum(delta_fi_list[0:m]))
        Y[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.sin(sum(delta_fi_list[0:m]))

X_list = [i for i in X.flat]
Y_list = [i for i in Y.flat]
Func_matrix_list = [l for l in Func_matrix_remake.flat]
viscosity_list = [m for m in viscosity_matrix.flat]

xi = np.linspace(min(X_list),max(X_list), 700)
yi = np.linspace(min(Y_list), max(Y_list), 700)
xig, yig = np.meshgrid(xi, yi)
bound_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
viscosity_i = interpolate.griddata((X_list, Y_list), viscosity_list, (xig, yig), method='cubic')

#fig, axs = plt.subplots(1, 2, figsize=(11, 10))

surf = plt.contourf(xig, yig, viscosity_i, linewidth=0.2, cmap=plt.get_cmap('jet'))
plt.show()
#plt.colorbar(surf, ax=axs[0, 0])

surf_1 = plt.contourf(xig, yig, bound_i, linewidth=0.2, cmap=plt.get_cmap('jet'))
#plt.colorbar(surf_1, ax=axs[0, 1])

plt.show()


print(np.shape(viscosity_matrix))
plt.plot(viscosity_matrix[:, 6])
plt.show()
print(min(viscosity_matrix.flat), max(viscosity_matrix.flat))

plt.plot(Func_matrix_remake[:, 6])
plt.show()

if __name__ == '__main__':
    velocity_in_point_1 = []
    velocity_in_point_2 = []
    velocity_in_point_3 = []
    pressure_in_point_1 = []
    pressure_in_point_2 = []
    pressure_in_point_3 = []
    dataFrame_to_study = pd.DataFrame()

    # dataFrame_delta_r = pd.DataFrame()
    # dataFrame_delta_r['deltaRlist'] = delta_r_list
    # dataFrame_delta_r.to_csv(
    #     'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframe_deltaR.csv')

    for t in range(T_exp):
        print(t)
        Pres_distrib, A, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, q_coef, wells_coord, delta_r_list, delta_fi_list, viscosity_matrix, porosity, C_total, perm, delta_t,
                         r_well)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))
        Func_matrix, velocity = define_func_matrix(Pres_distrib, Func_matrix_remake, perm, delta_r_list, delta_fi_list, delta_t, r_well, q, viscosity_matrix)
        bound_coords_rad_new, bound_coords_cart_new = find_bound_coords(Func_matrix, coord_matrix_rad, delta_r_list, delta_fi_list, M_fi_full, N_r_full)

        # for coord_pair in bound_coords_cart_new:
        #     plt.scatter(coord_pair[0], coord_pair[1])
        # plt.show()

        Func_matrix_remake = find_func_matrix_remake(coord_matrix_cart, M_fi_full, N_r_full, bound_coords_cart_new)
        #bound_coords_rad_new, bound_coords_cart_new = find_bound_coords(Func_matrix_remake, coord_matrix_rad, delta_r_list, delta_fi_list, M_fi_full, N_r_full)

        # for coord_pair in bound_coords_cart_new:
        #     plt.scatter(coord_pair[0], coord_pair[1])
        # plt.show()

        viscosity_matrix = find_viscosity(mu_oil, mu_water, Func_matrix_remake, coord_matrix_rad, delta_r_list,
                                          delta_fi_list, N_r_full, M_fi_full)
        #print(t, Pres_distrib[:, int(M_fi_full / 8 * 3)])
        # plt.plot(Pres_distrib[:, int(M_fi_full / 8 * 3)])
        # plt.show()
        #print(t, velocity[:, int(M_fi_full / 8 * 3)])
        # plt.plot( velocity[:, int(M_fi_full / 8 * 3)])
        # plt.show()
        #print(t, Func_matrix_remake[:, int(M_fi_full / 8 * 3)])
        # plt.plot(Func_matrix_remake[:, int(M_fi_full / 8 * 3)])
        # plt.show()
        #print(t, Func_matrix[:, int(M_fi_full / 8 * 3)])
        #print(t, viscosity_matrix[:, int(M_fi_full / 8 * 3)])
        # plt.plot(viscosity_matrix[:, int(M_fi_full / 8 * 3)])
        # plt.show()
        dataFrame_to_study['vel '+str(t)] = velocity[0:N_r_full, int(M_fi_full / 8 * 3)]
        dataFrame_to_study['Pres ' + str(t)] = Pres_distrib[0:N_r_full, int(M_fi_full / 8 * 3)]
        dataFrame_to_study['FMat ' + str(t)] = Func_matrix[0:N_r_full, int(M_fi_full / 8 * 3)]
        dataFrame_to_study['FMatRem ' + str(t)] = Func_matrix_remake[0:N_r_full, int(M_fi_full / 8 * 3)]
        dataFrame_to_study['Visc ' + str(t)] = viscosity_matrix[0:N_r_full, int(M_fi_full / 8 * 3)]

        if t == 9 or t == 10 or t == 19 or t == 20 or t == 49 or t == 50 or t == 89 or t == 90 or t == 88:
            print('print_dataframe')
            dataFrame_vel = pd.DataFrame(velocity)
            dataFrame_visc = pd.DataFrame(viscosity_matrix)
            dataFrame_func_matrix = pd.DataFrame(Func_matrix)
            dataFrame_func_matrix_remake = pd.DataFrame(Func_matrix_remake)
            dataFrame_Pres_distrib = pd.DataFrame(Pres_distrib)

            dataFrame_vel.to_csv(
                'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframes_QinCenter\\dataframe_vel_t_'+str(t)+'.csv')
            dataFrame_visc.to_csv(
                'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframes_QinCenter\\dataframe_visc_t_'+str(t)+'.csv')
            dataFrame_func_matrix.to_csv(
                'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframes_QinCenter\\dataframe_func_matr_t_'+str(t)+'.csv')
            dataFrame_func_matrix_remake.to_csv(
                'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframes_QinCenter\\dataframe_func_matr_rem_t_'+str(t)+'.csv')
            dataFrame_Pres_distrib.to_csv('C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframes_QinCenter\\dataframe_P_res_t_'+str(t)+'.csv')
            
        #print(dataFrame_to_study)
        # plt.plot(velocity[:, int(M_fi_full / 8 * 3)])
        # plt.show()
        # plt.plot(Func_matrix[:, int(M_fi_full / 8 * 3)])
        # plt.plot(Func_matrix_remake[:, int(M_fi_full / 8 * 3)])
        # plt.show()

        #old_bound = copy.deepcopy(new_bound)
        velocity_in_point_1.append(velocity[1, int(M_fi_full/8*3)])
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

    dataFrame_to_study.to_csv('C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframes_QinCenter\\dataframe.csv')
    plt.figure()
    plt.plot(velocity[:, int(M_fi_full / 8 * 3)])
    plt.show()

    #print(new_bound_G)

    fig1 = plt.figure()
    print(Pres_distrib[:, int(M_fi_full/2)])
    plt.plot(Pres_distrib[:, int(M_fi_full/8*3)])



    plt.xlabel("Номер ячейки")
    plt.ylabel('Давление, Па')
    #plt.ylim(100000, 1200000)
    P_all_center = np.ones((1, M_fi_full)) * P_center
    #Pres_distrib = np.vstack((P_all_center, Pres_distrib))
    plt.show()


    plt.figure()
    plt.plot(velocity_in_point_1)
    plt.show()

    plt.figure()
    plt.plot(velocity_in_point_2)
    plt.show()

    plt.figure()
    plt.plot(velocity_in_point_3)
    plt.show()

    plt.figure()
    plt.plot(pressure_in_point_1)
    plt.show()

    plt.figure()
    plt.plot(pressure_in_point_2)
    plt.show()

    plt.figure()
    plt.plot(pressure_in_point_3)
    plt.show()


    X = np.zeros((N_r_full,M_fi_full))
    Y = np.zeros((N_r_full, M_fi_full))
    for m in range(M_fi_full):
        for n in range(N_r_full):
            X[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.cos(sum(delta_fi_list[0:m]))
            Y[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.sin(sum(delta_fi_list[0:m]))

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]
    P_list = [k for k in Pres_distrib.flat]
    Func_matrix_list = [l for l in Func_matrix.flat]
    velocity_list = [m for m in velocity.flat]
    CP_list = zip(X_list, Y_list, P_list)

    print(min(P_list), max(P_list))
    print(M_fi_full/2)
    xi = np.linspace(min(X_list),max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')
    bound_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
    velocity_i = interpolate.griddata((X_list, Y_list), velocity_list, (xig, yig), method='cubic')

    fig, axs = plt.subplots(2, 2, figsize=(11, 10))
    #axs[0, 0].set_title(str((round(points_r_fi_pair[1][0], 2), round(points_r_fi_pair[1][1], 2))), y=0.75, loc='left')
    #axs[0, 0].set_xlabel("Time, s")
    #axs[0, 0].set_ylabel("Pressure, MPa")
    #axs[0, 0].plot(t_list_exp, P_filter_1)

    #levels = list(range(0,1500000,10000))
    #fig, ax = plt.subplots()
    #surf = plt.contourf(xig, yig, Pi, linewidth=0.2, cmap=plt.get_cmap('jet'), levels=levels)
    surf = axs[0, 0].contourf(xig, yig, Pi, linewidth=0.2, cmap=plt.get_cmap('jet'))
    plt.colorbar(surf, ax=axs[0, 0])
    #axs[0, 0].colorbar(surf)
    #plt.show()

    #fig = plt.figure()
    surf_1 = axs[0, 1].contourf(xig, yig, bound_i, linewidth=0.2, cmap=plt.get_cmap('jet'))
    plt.colorbar(surf_1, ax=axs[0, 1])
    #fig.colorbar(surf)
    #plt.scatter(new_bound_G[0], new_bound_G[1], marker='.')
    #plt.show()

    #levels = [0.00001*i for i in range(-5, 2,1)]
    #fig = plt.figure()
    surf_2 = axs[1, 0].contourf(xig, yig, velocity_i, linewidth=0.2, cmap=plt.get_cmap('jet'))
    plt.colorbar(surf_2, ax=axs[1, 0])
    #fig.colorbar(surf)
    #plt.show()

    # fig = plt.figure()
    # for coord_pair in new_bound_G:
    #     plt.scatter(coord_pair[0]*delta_r_fine*np.cos(coord_pair[1]*delta_fi_fine), coord_pair[0]*delta_r_fine*np.sin(coord_pair[1]*delta_fi_fine), marker='.')
    # plt.show()


    # patches = []
    # for i in range(len(delta_r_list)):
    #     circle = Wedge((0,0), r_well+sum(delta_r_list[0:i]), 0, 360, width=0.001)
    #     patches.append(circle)
    # for i in range(len(delta_fi_list)):
    #     x, y = np.array([[0, (R-r_well)*np.cos(sum(delta_fi_list[0:i]))], [0, (R-r_well)*np.sin(sum(delta_fi_list[0:i]))]])
    #     line = mlines.Line2D(x, y, lw=0.5)
    #     ax.add_line(line)
#
#
#     #    t = np.arange(0, 2 * np.pi, 0.01)
# #    r = 0.215
# #    plt.plot(r * np.sin(t) + Lx/2, r * np.cos(t) + Ly/2)
#     #ax = fig.gca(projection='3d')
#
#     #surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
#
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    # p = PatchCollection(patches)
    # ax.add_collection(p)
    #plt.show()

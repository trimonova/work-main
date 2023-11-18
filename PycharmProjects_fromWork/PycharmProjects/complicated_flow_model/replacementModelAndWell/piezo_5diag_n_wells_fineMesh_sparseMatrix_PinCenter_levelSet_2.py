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
from start_to_do_replacement_for_levelSet import replace_boundary
from newFuncMatrix_fix_3 import define_func_matrix
import copy


perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2*10**(-1)
fi = 0.4  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_oil*fi*(Cf+Cr)/perm
frac_angle = np.pi/6
frac_angle_2 = np.pi/6*7
delta_r = 0.0002
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


delta_t = 0.5
Pres = 1*10**5
P_center = 60*10**5
Pres_distrib = np.ones((N_r_full, M_fi_full)) * Pres
c3_oil = k_oil/delta_t*25
c3_water = k_water/delta_t*25
T_exp = 20
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

bound_coord, Func_matrix, bound_coord_cell, Func_coord_dict = replace_boundary(r_well, delta_fi_list, delta_r_list)
old_bound = bound_coord_cell
print(old_bound)
for coord_pair in set(bound_coord):
    plt.scatter(coord_pair[0], coord_pair[1])
plt.show()


Func_matrix_remake = Func_matrix
print(min(Func_matrix.flat), max(Func_matrix.flat))
print(N_r_full, M_fi_full)
def PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, P_center, wells_frac_coords, delta_r, delta_fi, Func_matrix):

    # пластовое давление во всей области на нулевом временном шаге
    B = np.zeros((N_r_full*M_fi_full, 1))
    for m in range(M_fi_full):

        A = np.zeros((N_r_full, N_r_full))

        for n in range(1, N_r_full - 1):
            c1 = 1 / delta_r_list[n] ** 2
            c2 = 1 / 2 / delta_r_list[n]
            A[n][n - 1] = c1 - c2 / (sum(delta_r_list[0:n + 1]))
            if Func_matrix[n][m] >= 0:
                A[n][n] = -2 * c1 - c3_oil - 2 / (sum(delta_r_list[0:n + 1])) ** 2 / delta_fi_list[m] ** 2
            else:
                A[n][n] = -2 * c1 - c3_water - 2 / (sum(delta_r_list[0:n + 1])) ** 2 / delta_fi_list[m] ** 2
            A[n][n + 1] = c1 + c2 / (sum(delta_r_list[0:n + 1]))


        c1 = 1 / delta_r_list[0] ** 2
        c2 = 1 / 2 / delta_r_list[0]

        A[0][0] = -2 * c1 - c3_oil - 2/(delta_r_list[0])**2/delta_fi_list[m]**2

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
            if n == 0:
                c1 = 1 / delta_r_fine**2
                c2 = 1 / 2 / delta_r_fine
                B[j][0] = -c3_oil * Pres_distrib[n][m] - (c1 - c2 / delta_r_fine)*P_center
            else:
                if Func_matrix[n][m] >= 0:
                    B[j][0] = -c3_oil * Pres_distrib[n][m]
                else:
                    B[j][0] = -c3_water * Pres_distrib[n][m]
                #print(c3_oil, c3_water, Pres_distrib[n][m])
            j += 1

    def sort_func(well_coord_couple):
        return (well_coord_couple[1]) * N_r_full + well_coord_couple[0]

    wells_frac_coords.sort(key=sort_func)
    wells_frac_coords_reverse = wells_frac_coords[:: -1]

    A_full = A_full.toarray()
    # for (n,m) in bound_coord_cell:
    #     print(n,m)
    #     if m != 0 and m != M_fi_full-1:
    #         print(A_full[(m)*N_r_full + n][(m)*N_r_full + n])
    #         print(A_full[(m) * N_r_full + n][(m) * N_r_full + n-1])
    #         print(A_full[(m) * N_r_full + n][(m) * N_r_full + n+1])
    #         print(A_full[(m)*N_r_full + n][(m-1)*N_r_full + n])
    #         print(A_full[(m)*N_r_full + n][(m+1)*N_r_full + n])
    #         print(B[m*N_r_full + n][0])

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
    #print(N_r, M_fi, N_r_oil, np.shape(P_new))

    P_new = P_new.reshape(N_r_full*M_fi_full, 1)
    Pres_end = np.zeros((N_r_full, M_fi_full))
    j = 0
    for m in range(M_fi_full):
        for n in range(N_r_full):
            Pres_end[n][m] = P_new[j][0]
            j += 1

    return Pres_end, A, B

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
        Pres_distrib, A, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, P_center, wells_coord, delta_r, delta_fi, Func_matrix_remake)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))
        #plt.plot(Pres_distrib[:, int(M_fi_full / 8 * 3)])
        #plt.show()
        Func_matrix, velocity, new_bound, new_bound_coord, Func_matrix_remake = define_func_matrix(Pres_distrib, Func_matrix_remake, perm, mu_water, mu_oil, delta_r_list, delta_fi_list, delta_t, r_well, P_center, old_bound)
        print(t, velocity[0:20, int(M_fi_full / 8 * 3)])
        print(t, Func_matrix_remake[0:20, int(M_fi_full / 8 * 3)])
        print(t, Func_matrix[0:20, int(M_fi_full / 8 * 3)])
        dataFrame_to_study['vel '+str(t)] = velocity[0:N_r_full, int(M_fi_full / 8 * 3)]
        dataFrame_to_study['Pres ' + str(t)] = Pres_distrib[0:N_r_full, int(M_fi_full / 8 * 3)]
        dataFrame_to_study['FMat ' + str(t)] = Func_matrix[0:N_r_full, int(M_fi_full / 8 * 3)]
        dataFrame_to_study['FMatRem ' + str(t)] = Func_matrix_remake[0:N_r_full, int(M_fi_full / 8 * 3)]
        print(dataFrame_to_study)
        # plt.plot(velocity[:, int(M_fi_full / 8 * 3)])
        # plt.show()
        # plt.plot(Func_matrix[:, int(M_fi_full / 8 * 3)])
        # plt.show()
        old_bound = copy.deepcopy(new_bound)
        velocity_in_point_1.append(velocity[1, int(M_fi_full/8*3)])
        velocity_in_point_2.append(velocity[8, int(M_fi_full / 8 * 3)])
        velocity_in_point_3.append(velocity[9, int(M_fi_full / 8 * 3)])
        pressure_in_point_1.append(Pres_distrib[1, int(M_fi_full / 8 * 3)])
        pressure_in_point_2.append(Pres_distrib[8, int(M_fi_full / 8 * 3)])
        pressure_in_point_3.append(Pres_distrib[9, int(M_fi_full / 8 * 3)])
        print(min(Func_matrix.flat), max(Func_matrix.flat))
        print(min(Func_matrix_remake.flat), max(Func_matrix_remake.flat))
        print(min(velocity.flat), max(velocity.flat))

    dataFrame_to_study.to_csv('C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\dataframe_2.csv')
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
    plt.plot(velocity[:, int(M_fi_full / 8 * 3)])
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

    for i in range(len(bound_coord)):
        plt.scatter(bound_coord[i][0], bound_coord[i][1], marker='.')
        plt.scatter(new_bound_G[i][0], new_bound_G[i][1], marker='*')
    plt.show()

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

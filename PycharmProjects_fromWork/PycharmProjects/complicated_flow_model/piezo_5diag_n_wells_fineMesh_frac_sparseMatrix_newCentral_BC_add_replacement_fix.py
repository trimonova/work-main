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


perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2*10**(-1)
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_oil*fi*(Cf+Cr)/perm
frac_angle = np.pi/6
frac_angle_2 = np.pi/6*7
delta_r = 0.005
delta_r_fine = 0.005
R_for_fine = 0.015
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine/delta_r_fine)
delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*round((R-r_well-R_for_fine)/delta_r)
N_r_full = len(delta_r_list)

delta_fi = np.pi / 60 # угол-шаг в радианах
delta_fi_fine = np.pi/180
fi_for_fine = np.pi/6
M_fi_fine = round(fi_for_fine / delta_fi_fine)
print(int((frac_angle-fi_for_fine)/delta_fi))
print(int(M_fi_fine*2))
print(int((np.pi - 2 * fi_for_fine) / delta_fi))
print(int((np.pi - frac_angle - fi_for_fine)/delta_fi))
delta_fi_list_first = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi))
angle_lack = round((2*np.pi - sum(delta_fi_list_first))/delta_fi)
delta_fi_list = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi)+angle_lack)

M_fi_full = len(delta_fi_list)
print(sum(delta_fi_list))
print((2*np.pi - sum(delta_fi_list))/delta_fi)


x_f = 0.01


delta_t = 1
Pres = 1*10**5
P_center = 5*10**5
Pres_distrib = np.ones((N_r_full, M_fi_full)) * Pres
c3_oil = k_oil/delta_t
c3_water = k_water/delta_t
T_exp = 20
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r**2)/100
wells_coord_real = [(0.1, np.pi/6+np.pi), (0.1, np.pi/2), (0.07, np.pi/3*5)]
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
P_well = [1500000, 700000, 20000]

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]


frac_angle_cell = round((frac_angle-fi_for_fine)/delta_fi) + M_fi_fine
frac_angle_2_cell = frac_angle_cell + round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + 2*M_fi_fine
L_N_frac = 20
frac_coord_1 = [i for i in range(L_N_frac)]
frac_coord_2 = [i for i in range(L_N_frac)]
frac_pressure_1 = [1500000]*L_N_frac
frac_pressure_2 = [1500000]*L_N_frac
frac_pressure = frac_pressure_1 + frac_pressure_2
frac_coords = [(i, frac_angle_cell) for i in frac_coord_1] + [(j, frac_angle_2_cell) for j in frac_coord_2]

wells_frac_coords = wells_coord + frac_coords
print(wells_frac_coords)
for i in range(len(frac_coords)):
    CP_dict[frac_coords[i]] = frac_pressure[i]

def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

bound_coord, Func_matrix, bound_coord_cell = replace_boundary(frac_angle, frac_angle_2, r_well, delta_fi_list, delta_r_list, sum(delta_r_list[0:len(frac_coord_1)]))

def PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, P_center, wells_frac_coords):

    # пластовое давление во всей области на нулевом временном шаге
    print(N_r_full, M_fi_full)
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

        A_sym_right = A_sym.copy()
        A_sym_left = A_sym.copy()

        # sign = 0
        # for coord_cell in bound_coord_cell:
        #     if coord_cell[1] == m:
        #         n = coord_cell[0]
        #         if n == 0:
        #             raise TypeError
        #         coef_1 = 1/delta_r_list[n]
        #         coef_2 = perm/mu_oil
        #         coef_3 = perm/mu_water
        #         coef_4 = 1/delta_fi_list[m]/sum(delta_r_list[0:n+1])
        #         A[n][n] = -(coef_1+coef_4)*(coef_2+coef_3)
        #         A[n][n-1] = coef_2*coef_1
        #         A[n][n+1] = coef_3*coef_1
        #
        #         A_sym_right[n][n] = coef_3 * coef_4
        #         A_sym_left[n][n] = coef_2 * coef_4
        #
        #         sign = 1

        A_sym_coo = coo_matrix(A_sym)
        # A_sym_right_coo = coo_matrix(A_sym_right)
        # A_sym_left_coo = coo_matrix(A_sym_left)

        if m == 0:
            # if sign == 1:
            #     A_line_1 = hstack([A, A_sym_right_coo, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym_left_coo])
            #     A_full = coo_matrix(A_line_1)
            # else:
            A_line_1 = hstack([A, A_sym_coo, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym_coo])
            A_full = coo_matrix(A_line_1)
        elif m == M_fi_full-1:
            # if sign == 1:
            #     A_line_end = hstack([A_sym_right_coo, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym_left_coo, A])
            #     A_full = vstack([A_full, A_line_end])
            # else:
            A_line_end = hstack(
                    [A_sym_coo, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym_coo, A])
            A_full = vstack([A_full, A_line_end])
        else:
            # if sign == 1:
            #     A_line = hstack([np.zeros((N_r_full,N_r_full*(m-1))), A_sym_left_coo, A, A_sym_right_coo, np.zeros((N_r_full, N_r_full * M_fi_full - (3+(m-1)) * N_r_full))])
            #     A_full = vstack([A_full, A_line])
            # else:
            A_line = hstack([np.zeros((N_r_full, N_r_full * (m - 1))), A_sym_coo, A, A_sym_coo,
                                 np.zeros((N_r_full, N_r_full * M_fi_full - (3 + (m - 1)) * N_r_full))])
            A_full = vstack([A_full, A_line])

    # A_full = A_full.toarray()
    # for (n,m) in bound_coord_cell:
    #     print(n,m)
    #     if m != 0 and m != M_fi_full-1:
    #         coef_4 = 1 / delta_fi_list[m] / sum(delta_r_list[0:n])
    #         print(A_full[(m)*N_r_full + n][(m)*N_r_full + n])
    #         print(A_full[(m)*N_r_full + n][(m)*N_r_full + n])
    #         A_full[(m)*N_r_full + n][(m-1)*N_r_full + n] = coef_2*coef_4
    #         A_full[(m)*N_r_full + n][(m+1)*N_r_full + n] = coef_3* coef_4
    #         print(A_full[(m) * N_r_full + n][(m - 1) * N_r_full + n])
    #         print(A_full[(m) * N_r_full + n][(m+1) * N_r_full + n])
    #     elif m == 0:
    #         coef_4 = 1 / delta_fi_list[m] / sum(delta_r_list[0:n])
    #         A_full[(m) * N_r_full + n][(M_fi_full - 1) * N_r_full + n] = coef_2 * coef_4
    #         A_full[(m) * N_r_full + n][(m + 1) * N_r_full + n] = coef_3 * coef_4
    #     elif m == M_fi_full-1:
    #         coef_4 = 1 / delta_fi_list[m] / sum(delta_r_list[0:n])
    #         A_full[(m) * N_r_full + n][(m - 1) * N_r_full + n] = coef_2 * coef_4
    #         A_full[(m) * N_r_full + n][n] = coef_3 * coef_4
    #
    # A_full = coo_matrix(A_full)
    j = 0
    for m in range(M_fi_full):
        for n in range(N_r_full):
            # if n == 0 and (m == frac_angle_cell or m == frac_angle_2_cell) :
            c1 = 1 / delta_r_list[n] ** 2
            c2 = 1 / 2 / delta_r_list[n]
            if Func_matrix[n][m] > 0:
                B[j][0] = -c3_oil * Pres_distrib[n][m]
            else:
                B[j][0] = -c3_water * Pres_distrib[n][m]
            j += 1

    def sort_func(well_coord_couple):
        return (well_coord_couple[1]) * N_r_full + well_coord_couple[0]

    wells_frac_coords.sort(key=sort_func)
    print(wells_frac_coords)
    wells_frac_coords_reverse = wells_frac_coords[:: -1]
    print(wells_frac_coords_reverse)

    # for coord_cell in bound_coord_cell:
    #     #if coord_cell[1] != 0 and coord_cell[1] != M_fi_full-1:
    #     B[coord_cell[1] * N_r_full + coord_cell[0]][0] = 0

    A_full = A_full.toarray()
    for (n,m) in bound_coord_cell:
        print(n,m)
        if m != 0 and m != M_fi_full-1:
            print(A_full[(m)*N_r_full + n][(m)*N_r_full + n])
            print(A_full[(m) * N_r_full + n][(m) * N_r_full + n-1])
            print(A_full[(m) * N_r_full + n][(m) * N_r_full + n+1])
            print(A_full[(m)*N_r_full + n][(m-1)*N_r_full + n])
            print(A_full[(m)*N_r_full + n][(m+1)*N_r_full + n])
            print(B[m*N_r_full + n][0])

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

        # A_full = np.delete(A_full, (coord_couple[1]-1)*N_r + coord_couple[0], axis=0)
        # A_full = np.delete(A_full,
        #                    (coord_couple[1] - 1) * N_r + coord_couple[0],axis=1)
        B = np.delete(B, (coord_couple[1] - 1) * N_r_full + coord_couple[0], axis=0)
    P_new = spsolve(A_full, B)
    #print(np.shape(A_full), np.shape(B))
    #P_new = np.linalg.solve(A_full, B)
    #P_new = P_new_csr.toarray()
    for coord_couple in wells_frac_coords:
        P_new = np.insert(P_new, (coord_couple[1]-1)*N_r_full + coord_couple[0], CP_dict[coord_couple])
    #print(N_r, M_fi, N_r_oil, np.shape(P_new))
    print(N_r_full, M_fi_full)
    print(np.shape(P_new))
    P_new = P_new.reshape(N_r_full*M_fi_full, 1)
    Pres_end = np.zeros((N_r_full, M_fi_full))
    j = 0
    for m in range(M_fi_full):
        for n in range(N_r_full):
            print(P_new[j][0])
            Pres_end[n][m] = P_new[j][0]
            j += 1

    return Pres_end, A_full, B

if __name__ == '__main__':
    for t in range(T_exp):
        print(t)
        Pres_distrib, A_full, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, P_center, wells_frac_coords)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))

    fig1 = plt.figure()
    print(np.shape(Pres_distrib))
    print(Pres_distrib[:, int(M_fi_full/2)])
    plt.plot(Pres_distrib[:, int(M_fi_full/2)])
    P_all_center = np.ones((1, M_fi_full)) * P_center
    Pres_distrib = np.vstack((P_all_center, Pres_distrib))
    #plt.plot(Pres_distrib[:,30])
    plt.show()
    X = np.zeros((N_r_full+1,M_fi_full))
    Y = np.zeros((N_r_full+1, M_fi_full))
    for m in range(M_fi_full):
        for n in range(N_r_full+1):
            X[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.cos(sum(delta_fi_list[0:m]))
            Y[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.sin(sum(delta_fi_list[0:m]))

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
    fig, ax = plt.subplots()
    surf = plt.contourf(xig, yig, Pi, linewidth=0.2, cmap=plt.get_cmap('jet'), levels=levels)

    # for coord_pair in bound_coord:
    #     plt.scatter(coord_pair[0], coord_pair[1], marker='.')
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
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # p = PatchCollection(patches)
    # ax.add_collection(p)
    plt.show()

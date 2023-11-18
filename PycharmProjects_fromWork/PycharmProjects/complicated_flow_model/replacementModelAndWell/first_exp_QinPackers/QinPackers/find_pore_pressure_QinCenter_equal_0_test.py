import numpy as np
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve

# Решение уравнения пьезопроводности в цилиндрических координатах, двумерный случай (неявная схема)
# В центре скважина с постоянным давлением P, на границах задается: градиент давления равен нулю.
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
import copy

if __name__ == '__main__':

    perm = 2 * 10 ** (-15)  # м2 проницаемость
    mu_water = 1 * 10 ** (-3)  # Па*с вязкость
    mu_oil = 1 * 10 ** (-3)
    porosity = 0.19  # пористость
    Cf = 10 ** (-9)  # сжимаемость флюида
    Cr = 5 * 10 ** (-10)  # сжимаемость скелета
    C_total = (Cf + Cr) * 25
    k_water = mu_water * porosity * C_total / perm
    k_oil = mu_oil * porosity * C_total / perm
    frac_angle = np.pi / 4
    frac_angle_2 = np.pi / 4 * 5
    frac_length_1 = 0.01
    frac_length_2 = 0.01
    delta_r = 0.005
    delta_r_fine = 0.005
    R_for_fine = 0.02
    R = 0.215
    r_well = 0.0075
    N_r_fine = round(R_for_fine / delta_r_fine)
    delta_r_list = [delta_r_fine] * N_r_fine + [delta_r] * round((R - r_well - R_for_fine) / delta_r)
    N_r_full = len(delta_r_list)

    delta_fi = np.pi / 30  # угол-шаг в радианах
    delta_fi_fine = np.pi / 30
    fi_for_fine = np.pi / 6
    M_fi_fine = round(fi_for_fine / delta_fi_fine)

    delta_fi_list_first = [delta_fi] * round((frac_angle - fi_for_fine) / delta_fi) + [delta_fi_fine] * (
                M_fi_fine * 2) + [
                              delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [
                              delta_fi_fine] * (M_fi_fine * 2) + [
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

    delta_t = 10
    Pres = 1 * 10 ** 5
    # P_center = 60*10**5
    # Q_center = 0.2 * 10 ** (-8)  # из лаб. данных - 0.2*10**(-6) m3/s
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

    wells_coord = [(round(0.1711 / delta_r), round(np.pi / 4 / delta_fi)),
                   (round(0.1711 / delta_r), round((np.pi + np.pi / 4) / delta_fi))]
    P_well = [450000, 140000]
    print(wells_coord)

    CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

    for i in range(len(wells_coord)):
        CP_dict[wells_coord[i]] = P_well[i]


    def sortByRad(inputSet):
        return inputSet[0]


    def sortByAngle(inputSet):
        return inputSet[1]




    viscosity_matrix = np.ones((N_r_full, M_fi_full)) * mu_water


def PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, q_coef, wells_frac_coords,
                         delta_r_list, delta_fi_list, viscosity_matrix, fi, C_total, perm, delta_t,
                         r_well):

    # пластовое давление во всей области на нулевом временном шаге
    B = np.zeros((N_r_full*M_fi_full, 1))
    for m in range(M_fi_full):

        A = np.zeros((N_r_full, N_r_full))

        for n in range(1, N_r_full - 1):
            k_fluid = viscosity_matrix[n][m] * fi * C_total / perm
            c3_fluid = k_fluid/delta_t
            c1 = 1 / delta_r_list[n] ** 2
            c2 = 1 / 2 / delta_r_list[n]
            A[n][n - 1] = c1 - c2 / (sum(delta_r_list[0:n+1])+r_well)

            A[n][n] = -2 * c1 - c3_fluid - 2 / (sum(delta_r_list[0:n+1])+r_well) ** 2 / delta_fi_list[m] ** 2

            A[n][n + 1] = c1 + c2 / (sum(delta_r_list[0:n+1])+r_well)


        c1 = 1 / delta_r_list[0] ** 2
        c2 = 1 / 2 / delta_r_list[0]

        A[0][0] = -c1 - c3_oil - 2/(delta_r_list[0]+r_well)**2/delta_fi_list[m]**2 - c2/(delta_r_list[0]+r_well)

        A[0][1] = c1 + c2 / (delta_r_list[0]+r_well)

        c1 = 1 / delta_r_list[N_r_full - 1] ** 2
        c2 = 1 / 2 / delta_r_list[N_r_full - 1]
        A[N_r_full - 1][N_r_full - 1] = -2 * c1 - c3_water + c1 + c2 / (sum(delta_r_list[0:N_r_full])+r_well) - 2/(sum(delta_r_list[0:N_r_full])+r_well)**2/delta_fi_list[m]**2
        A[N_r_full - 1][N_r_full - 2] = c1 - c2 / (sum(delta_r_list[0:N_r_full])+r_well)

        c4 = 1 / delta_fi_list[m] ** 2
        A_sym = np.zeros((N_r_full, N_r_full))
        for n in range(0, N_r_full):
            A_sym[n][n] = c4 / (sum(delta_r_list[0:n+1])+r_well) ** 2


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
                c1 = 1 / delta_r_list[n]**2
                c2 = 1 / 2 / delta_r_list[n]
                B[j][0] = -c3_oil * Pres_distrib[n][m] - (c1 - c2 / (delta_r_list[n] + r_well))*q_coef
            else:
                k_fluid = viscosity_matrix[n][m] * fi * C_total / perm
                c3_fluid = k_fluid / delta_t
                #print(c3_fluid, Pres_distrib[n][m], viscosity_matrix[n][m], fi, C_total, perm)
                B[j][0] = -c3_fluid * Pres_distrib[n][m]

            j += 1

    def sort_func(well_coord_couple):
        return (well_coord_couple[1]) * N_r_full + well_coord_couple[0]

    wells_frac_coords.sort(key=sort_func)
    wells_frac_coords_reverse = wells_frac_coords[:: -1]

    #A_full = A_full.toarray()
    # for (n,m) in bound_coord_cell:
    #     print(n,m)
    #     if m != 0 and m != M_fi_full-1:
    #         print(A_full[(m)*N_r_full + n][(m)*N_r_full + n])
    #         print(A_full[(m) * N_r_full + n][(m) * N_r_full + n-1])
    #         print(A_full[(m) * N_r_full + n][(m) * N_r_full + n+1])
    #         print(A_full[(m)*N_r_full + n][(m-1)*N_r_full + n])
    #         print(A_full[(m)*N_r_full + n][(m+1)*N_r_full + n])
    #         print(B[m*N_r_full + n][0])

    #A_full = coo_matrix(A_full)
    #A_array = A_full.toarray()
    # for coord_couple in wells_frac_coords_reverse:
    #     print(A_full[(coord_couple[1]-1)*N_r_full + coord_couple[0]])

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
    print(np.min(A_full), np.max(A_full))
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
    PinPoint = []
    for t in range(T_exp_dir):
        print(t)
        Pres_distrib, A_full, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_oil, c3_water, CP_dict, q_coef, wells_coord,
                         delta_r_list, delta_fi_list, viscosity_matrix, porosity, C_total, perm, delta_t, r_well)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))
        print('center', Pres_distrib[1][11])
        PinPoint.append(Pres_distrib[2][2])

    fig1 = plt.figure()
    print(np.shape(Pres_distrib))
    print(Pres_distrib[:, int(M_fi_full/2)])
    plt.plot(Pres_distrib[:,wells_coord[0][1]])
    # P_all_center = np.ones((1, M_fi)) * P_center
    # Pres_distrib = np.vstack((P_all_center, Pres_distrib))
    plt.show()
    print(PinPoint)
    fig2 = plt.figure()
    plt.plot(PinPoint)
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
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)

#    t = np.arange(0, 2 * np.pi, 0.01)
#    r = 0.215
#    plt.plot(r * np.sin(t) + Lx/2, r * np.cos(t) + Ly/2)
    #ax = fig.gca(projection='3d')

    #surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

#  Можно задавать много точек с источниками и давлениями. Если задавать только давления (граничные условия) то задача устойчива при любых шагах времени и координаты,
# если задавать еще источники, то задача устойчива при каком-то соотношении t_step и hx, Q задается вроде бы в м3/с
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve


if __name__ == '__main__':
    perm = 2 * 10 ** (-15)  # м2 проницаемость
    mu = 2 * 10 ** (-3)  # Па*с вязкость
    fi = 0.2  # пористость
    Cf = 10 ** (-9)  # сжимаемость флюида
    Cr = 5 * 10 ** (-10)  # сжимаемость скелета
    k = mu*fi*(Cf+Cr)/perm

    hx = 0.05
    hy = 0.05

    t_step = 0.1
    T_exp = 40
    Lx = 2
    Ly = 2

    Courant_number = t_step/k/hx**2 + t_step/k/hy**2
    print(Courant_number)

    N = int(Lx/hx) # количество ячеек вдоль оси х
    M = int(Ly/hy)
    #wells_with_Q = {(int(0.430/hx),int(0.430/hy)): -0.000003}
    wells_with_Q = {}
    #wells_with_Q = {(int(0.309/hx),int(0.309/hy)): -0.00001, (int(0.309/hx), int(0.551/hy)): -0.00001, (int(0.551/hx),int(0.309/hy)): -0.00001}
    wells_coord = [(int(round(0.5/hx)+15), int(round(1.1489/hy+5)))]
    P_well = [500000]

    # wells_coord = []
    # P_well = []
    print(wells_coord)

    CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

    for i in range(len(wells_coord)):
        CP_dict[wells_coord[i]] = P_well[i]

    frac_pressure = [1000000, 1100000, 1200000, 1100000, 1000000]
    frac_coords = [(18, 20), (19, 20), (20, 20), (21, 20), (22, 20)]
    # frac_pressure = []
    # frac_coords = []
    wells_frac_coords = wells_coord + frac_coords
    print(wells_frac_coords)
    for i in range(len(frac_coords)):
        CP_dict[frac_coords[i]] = frac_pressure[i]

    Pres = 1*10**5 # давление в пласте
    Pbound = 1*10**5 #  давление на границе\

    Pres_distrib = np.ones((N, M)) * Pres
    P_add_hor = np.ones((N, 1)) * Pres
    P_add_vert = np.ones((1, M + 2)) * Pres



def PorePressure_in_Time(Pres_distrib, N, M, hx, hy, k, t_step, Pres, wells_frac_coords, CP_dict):

    A = np.zeros((N, N))
    B = np.zeros((N * M, 1))


    for n in range(1, N-1):
        A[n][n-1] = 1/hx**2
        A[n][n] = (-2/hx**2 - k/t_step - 2/hy**2)
        A[n][n+1] = 1/hx**2

    A[0][0] = (-2/hx**2 - k/t_step - 2/hy**2)
    A[0][1] = 1/hx**2
    A[N-1][N-1] = A[0][0]
    A[N-1][N-2] = A[0][1]

    A_sym = np.zeros((N, N))
    for n in range(0, N):
        A_sym[n][n] = 1/hy**2
    A_sym_coo = coo_matrix(A_sym)
    print(np.shape(A), np.shape(A_sym_coo), np.shape(np.zeros((N, N * M - 2 * N))))
    A_line_1 = hstack([A, A_sym_coo, np.zeros((N, N * M - 2 * N))])
    A_full = coo_matrix(A_line_1)

    for m in range(1, M - 1):
        A_line = hstack(
            (np.zeros((N, N * (m - 1))), A_sym_coo, A, A_sym_coo, np.zeros((N, N * M - (3 + (m - 1)) * N))))
        A_full = vstack((A_full, A_line))

    A_line_end = hstack((np.zeros((N, N * M - 2 * N)), A_sym_coo, A))
    A_full = vstack((A_full, A_line_end))

    j = 0
    for m in range(M):
        for n in range(N):
            if m == 0 or m == M-1:
                B[j][0] = -k/ t_step * Pres_distrib[n][m] - 1/hy**2 * Pres
                if n == 0 or n == N-1:
                    B[j][0] = B[j][0] - 1/hx**2 * Pres
            elif n == 0 or n == N-1:
                B[j][0] = -k / t_step * Pres_distrib[n][m] - 1 / hx ** 2 * Pres
            else:
                B[j][0] = -k/t_step * Pres_distrib[n][m]
            j += 1

    wells_frac_coords.sort()
    wells_coord_reverse = wells_frac_coords[:: -1]
    for well_coord_couple in wells_coord_reverse:
        A_well_column_coo = A_full.getcol((well_coord_couple[1] - 1) * N + well_coord_couple[0])
        A_well_column = A_well_column_coo.toarray()
        for cell_number in range(len(A_well_column)):
            if A_well_column[cell_number] != 0:
                B[cell_number] = B[cell_number] - A_well_column[cell_number] * CP_dict[well_coord_couple]

        A_full = A_full.tocsr()
        all_cols = np.arange(A_full.shape[1])
        cols_to_keep = np.where(np.logical_not(np.in1d(all_cols, [(well_coord_couple[1] - 1) * N + well_coord_couple[0]])))[0]
        A_full = A_full[:, cols_to_keep]
        A_full = A_full[cols_to_keep, :]
        # A_full = np.delete(A_full, (well_coord_couple[1] - 1) * N + well_coord_couple[0], axis=0)
        # A_full = np.delete(A_full,
        #                    (well_coord_couple[1] - 1) * N + well_coord_couple[0], axis=1)
        B = np.delete(B, (well_coord_couple[1] - 1) * N + well_coord_couple[0], axis=0)

    # print(np.shape(A_full), np.shape(B))
    # P_new = np.linalg.solve(A_full, B)
    P_new = spsolve(A_full, B)
    print(min(P_new), max(P_new))
    for well_coord_couple in wells_frac_coords:
        print((well_coord_couple[1] - 1) * N + well_coord_couple[0], CP_dict[well_coord_couple])
        print(well_coord_couple[1], N, well_coord_couple[0])
        P_new = np.insert(P_new, (well_coord_couple[1] - 1) * N + well_coord_couple[0], CP_dict[well_coord_couple])
    # print(N_r, M_fi, N_r_oil, np.shape(P_new))
    P_new = P_new.reshape(N * M, 1)
    Pres_end = np.zeros((N, M))
    j = 0
    for m in range(M):
        for n in range(N):
            Pres_end[n][m] = P_new[j][0]
            j += 1

                    # for coord_key in wells_with_Q:
        #     if (n,m) == coord_key:
        #         B[n-1][0] = -V*beta/t_step*Pres_distrib[n][m]- alpha*coeff_1*(Pres_distrib[n][m-1] - 2*Pres_distrib[n][m] + Pres_distrib[n][m+1]) + wells_with_Q[coord_key]


    #P_total = np.hstack((P_total,P_new))


    # P_total = np.hstack((P_add_hor, Pres_end, P_add_hor))
    # P_total = np.vstack((P_add_vert, P_total, P_add_vert))
    #print(np.shape(P_total))

    # Pres_end = np.array(P_total.copy())

    return Pres_end


#----------------------------------------------------------------------------
if __name__ == '__main__':
    for t in range(T_exp):
        print(t)
        Pres_distrib = PorePressure_in_Time(Pres_distrib,  N, M, hx, hy, k, t_step, Pres, wells_frac_coords)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))

    X = np.zeros((N,M))
    Y = np.zeros((N, M))
    for m in range(M):
        for n in range(N):
            X[n][m] = n*hx
            Y[n][m] = m*hy

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]
    P_list = [k for k in Pres_distrib.flat]


    CP_list = zip(X_list, Y_list, P_list)

    print(min(P_list), max(P_list))

    xi = np.linspace(min(X_list),max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')

    levels = list(range(0,1500000,50000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)
    #ax = fig.gca(projection='3d')

    #surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

import numpy as np
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy import interpolate
from matplotlib import cm
import matplotlib.pyplot as plt
perm = 2 * 10 ** (-12)  # м2 проницаемость
mu_water = 2 * 10 ** (-5)  # Па*с вязкость

porosity = 0.4  # пористость
k_water = perm/2/mu_water/porosity
L_x = 0.2
L_y = 0.4
delta_x = 0.01
delta_y = 0.01

N_r_full = round(L_x/delta_x)

M_fi_full = round(L_y/delta_y)

delta_t = 0.5
Pres = 1*10**5
Pres_distrib = np.ones((N_r_full, M_fi_full)) * Pres
P_bottom = 20*10**5
P_top = 1*10**5

T = 100

viscosity_matrix = np.ones((N_r_full, M_fi_full))*mu_water

wells_coord = []
P_well = []
CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений
for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

def PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, P_bottom, P_top, CP_dict, wells_frac_coords,
                         viscosity_matrix, fi, perm, delta_t):

    # пластовое давление во всей области на нулевом временном шаге
    B = np.zeros((N_r_full*M_fi_full, 1))
    for m in range(M_fi_full):

        A = np.zeros((N_r_full, N_r_full))

        for n in range(1, N_r_full - 1):
            alpha = perm/2/fi/viscosity_matrix[n][m]

            A[n][n - 1] = Pres_distrib[n][m]/delta_x**2
            A[n][n] = -2 * Pres_distrib[n][m]*(1/delta_x**2 + 1/delta_y**2) - 1/2/delta_t/alpha
            A[n][n + 1] = A[n][n - 1]

        alpha = perm/2/fi/viscosity_matrix[0][m]
        A[0][0] = -2 * Pres_distrib[0][m]*(1/delta_x**2 + 1/delta_y**2) - 1/2/delta_t/alpha + Pres_distrib[0][m]/delta_x**2
        A[0][1] = Pres_distrib[0][m]/delta_x**2

        alpha = perm/2/fi/viscosity_matrix[N_r_full-1][m]
        A[N_r_full - 1][N_r_full - 1] = -2 * Pres_distrib[N_r_full-1][m]*(1/delta_x**2 + 1/delta_y**2) - 1/2/delta_t/alpha + Pres_distrib[N_r_full-1][m]/delta_x**2
        A[N_r_full - 1][N_r_full - 2] = Pres_distrib[N_r_full-1][m]/delta_x**2

        A_sym = np.zeros((N_r_full, N_r_full))
        for n in range(0, N_r_full):
            A_sym[n][n] = Pres_distrib[n][m] / delta_y ** 2


        A_sym_coo = coo_matrix(A_sym)

        if m == 0:
            A_line_1 = hstack([A, A_sym_coo, np.zeros((N_r_full, N_r_full * M_fi_full - 2 * N_r_full))])
            A_full = coo_matrix(A_line_1)
        elif m == M_fi_full-1:
            A_line_end = hstack(
                    [np.zeros((N_r_full, N_r_full * M_fi_full - 2 * N_r_full)), A_sym_coo, A])
            A_full = vstack([A_full, A_line_end])
        else:
            A_line = hstack([np.zeros((N_r_full, N_r_full * (m - 1))), A_sym_coo, A, A_sym_coo,
                                 np.zeros((N_r_full, N_r_full * M_fi_full - (3 + (m - 1)) * N_r_full))])
            A_full = vstack([A_full, A_line])

    j = 0
    for m in range(0, M_fi_full):
        for n in range(0, N_r_full):
            alpha = perm / 2 / fi / viscosity_matrix[n][m]
            if m == 0 and n == 0:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n + 1][m] - Pres_distrib[n][m]) / delta_x) ** 2 - (
                                      (Pres_distrib[n][m + 1] - Pres_distrib[n][m]) / delta_y) ** 2 - Pres_distrib[n][m]/delta_y**2*P_bottom
            elif m == 0 and n == N_r_full-1:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n][m] - Pres_distrib[n - 1][m]) / delta_x) ** 2 - (
                                      (Pres_distrib[n][m + 1] - Pres_distrib[n][m]) / delta_y) ** 2 - Pres_distrib[n][m]/delta_y**2*P_bottom
            elif m == M_fi_full - 1 and n == 0:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n + 1][m] - Pres_distrib[n][m]) / delta_x) ** 2 - (
                                      (Pres_distrib[n][m] - Pres_distrib[n][m - 1]) / delta_y) ** 2 - Pres_distrib[n][m]/delta_y**2*P_top
            elif m == M_fi_full -1 and n == N_r_full-1:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n][m] - Pres_distrib[n - 1][m]) / delta_x) ** 2 - (
                                      (Pres_distrib[n][m] - Pres_distrib[n][m - 1]) / delta_y) ** 2 - Pres_distrib[n][m]/delta_y**2*P_top
            elif m == 0:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n + 1][m] - Pres_distrib[n - 1][m]) / 2 / delta_x) ** 2 - (
                                      (Pres_distrib[n][m + 1] - Pres_distrib[n][m]) / delta_y) ** 2 - Pres_distrib[n][m]/delta_y**2*P_bottom
            elif m == M_fi_full-1:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n + 1][m] - Pres_distrib[n - 1][m]) / 2 / delta_x) ** 2 - (
                                      (Pres_distrib[n][m] - Pres_distrib[n][m - 1]) / delta_y) ** 2 - Pres_distrib[n][m]/delta_y**2*P_top
            elif n == 0:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n + 1][m] - Pres_distrib[n][m]) / delta_x) ** 2 - (
                                      (Pres_distrib[n][m + 1] - Pres_distrib[n][m - 1]) / 2 / delta_y) ** 2
            elif n == N_r_full-1:
                B[j][0] = -1 / 2 / alpha / delta_t * Pres_distrib[n][m] - (
                            (Pres_distrib[n][m] - Pres_distrib[n - 1][m]) / delta_x) ** 2 - (
                                      (Pres_distrib[n][m + 1] - Pres_distrib[n][m - 1]) / 2 / delta_y) ** 2
            else:
                B[j][0] = -1/2/alpha/delta_t*Pres_distrib[n][m] - ((Pres_distrib[n+1][m]-Pres_distrib[n-1][m])/2/delta_x)**2 - ((Pres_distrib[n][m+1]-Pres_distrib[n][m-1])/2/delta_y)**2
            j += 1


    def sort_func(well_coord_couple):
        return (well_coord_couple[1]) * N_r_full + well_coord_couple[0]

    wells_frac_coords.sort(key=sort_func)
    wells_frac_coords_reverse = wells_frac_coords[:: -1]


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

for t in range(T):
    print(t)
    Pres_distrib, A, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, P_bottom, P_top, CP_dict, wells_coord,
                         viscosity_matrix, porosity, perm, delta_t)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))

X = np.zeros((N_r_full,M_fi_full))
Y = np.zeros((N_r_full, M_fi_full))
for m in range(M_fi_full):
    for n in range(N_r_full):
        X[n][m] = delta_x
        Y[n][m] = delta_y

X_list = [delta_x*i for i in range(N_r_full)]
Y_list = [delta_y*i for i in range(M_fi_full)]
P_list = [k for k in Pres_distrib.flat]

print(np.shape(X_list), np.shape(Y_list), np.shape(P_list))
CP_list = zip(X_list, Y_list, P_list)

print(min(P_list), max(P_list))


plt.plot(Pres_distrib[int(N_r_full/2), :])
plt.show()

xi = np.linspace(min(X_list),max(X_list), 1000)
yi = np.linspace(min(Y_list), max(Y_list), 1000)
xig, yig = np.meshgrid(xi, yi)
#Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')
#levels = list(range(0,max(P_list),10000))
fig = plt.figure()
#surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)
print(len(X_list), len(Y_list), np.shape(Pres_distrib))
surf = plt.contourf(Y_list, X_list, Pres_distrib, cmap=cm.jet, antialiased=True, vmin=np.nanmin(P_list), vmax=np.nanmax(P_list),linewidth=0.2)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
# Решение уравнения пьезопроводности в цилиндрических координатах, двумерный случай (неявная схема)
# В центре скважина с постоянным давлением P, на границах задается: градиент давления равен нулю.
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve

perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2*10**(-2)
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_oil*fi*(Cf+Cr)/perm
delta_r = 0.002
delta_fi = np.pi / 36 # угол-шаг в радианах
R = 0.215
r_well = 0.0075
x_f = 0.05 
N_r_oil = int(x_f/delta_r)
N_r = int((R-r_well)/delta_r)
M_fi = int(2*np.pi/delta_fi)
delta_t = 0.1
Pres = 1*10**5
P_center = 20*10**5
Q = 2*10**(-5)
s = 0.003*0.01
q = Q/s/10000/2.5
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3_oil = k_oil/delta_t
c3_water = k_water/delta_t
c4 = 1/delta_fi**2
T_exp = 400
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r**2)/100
print(Courant_number_oil, Courant_number_water)
wells_coord = [(int(0.17/delta_r), int(np.pi/4/delta_fi)), (int(0.17/delta_r), int(5*np.pi/4/delta_fi))]
P_well = [2000000, 100000]

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

#for i in range(len(wells_coord)):
#    Pres_distrib[wells_coord[i][0]][wells_coord[i][1]] = P_well[i]

def PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3_oil, c3_water, c4, wells_coord, CP_dict, delta_r, P_center, N_r_oil):

    # пластовое давление во всей области на нулевом временном шаге
    P_total = np.ones((N_r, 1)) * Pres
    print(N_r, M_fi)
    A = np.zeros((N_r, N_r))
    B = np.zeros((N_r*M_fi, 1))
    for n in range(1, N_r_oil - 1):
        A[n][n - 1] = c1 - c2 / ((n + 1) * delta_r)
        A[n][n] = -2 * c1 - c3_oil - 2/((n+1)*delta_r)**2/delta_fi**2
        A[n][n + 1] = c1 + c2 / ((n + 1) * delta_r)

    for n in range(N_r_oil, N_r - 1):
        A[n][n - 1] = c1 - c2 / ((n + 1) * delta_r)
        A[n][n] = -2 * c1 - c3_water - 2/((n+1)*delta_r)**2/delta_fi**2
        A[n][n + 1] = c1 + c2 / ((n + 1) * delta_r)

    A[N_r_oil - 1][N_r_oil - 1] = 1 / delta_r / mu_oil + 1 / delta_r / mu_water
    A[N_r_oil - 1][N_r_oil] = -1 / delta_r / mu_water
    A[N_r_oil - 1][N_r_oil - 2] = -1 / delta_r / mu_oil
    A[0][0] = -2 * c1 - c3_oil - 2/(delta_r)**2/delta_fi**2 + c1 - c2 / delta_r
    A[0][1] = c1 + c2 / (1 * delta_r)
    A[N_r - 1][N_r - 1] = -2 * c1 - c3_water + c1 + c2 / ((N_r) * delta_r) - 2/((N_r)*delta_r)**2/delta_fi**2
    A[N_r - 1][N_r - 2] = c1 - c2 / ((N_r) * delta_r)

    A_sym = np.zeros((N_r, N_r))
    for n in range(0,N_r):
        A_sym[n][n] = c4/((n+1)*delta_r)**2
    A_sym[N_r_oil-1][N_r_oil-1] = 0

    A_sym_coo = coo_matrix(A_sym)
    A_line_1 = hstack((A, A_sym_coo, np.zeros((N_r, N_r*M_fi-3*N_r)), A_sym_coo))
    A_full = coo_matrix(A_line_1)

    for m in range(1, M_fi-1):
        A_line = hstack((np.zeros((N_r,N_r*(m-1))), A_sym_coo, A, A_sym_coo, np.zeros((N_r, N_r * M_fi - (3+(m-1)) * N_r))))
        A_full = vstack((A_full, A_line))

    A_line_end = hstack((A_sym_coo, np.zeros((N_r, N_r*M_fi-3*N_r)), A_sym_coo, A))
    A_full = vstack((A_full, A_line_end))

    j = 0
    for m in range(M_fi):
        for n in range(N_r):
            if n == 0:
                #print(j, n, m)
                B[j][0] = -c3_oil * Pres_distrib[n][m] - (c1 - c2/(delta_r))*q*mu_oil*delta_r/perm
            elif n < N_r_oil-1:
                B[j][0] = -c3_oil * Pres_distrib[n][m]
            elif n > N_r_oil-1:
                B[j][0] = -c3_water * Pres_distrib[n][m]
            elif n == N_r_oil-1:
                B[j][0] = 0
            j += 1

    def sort_func(well_coord_couple):
        return (well_coord_couple[1]) * N_r + well_coord_couple[0]

    wells_coord.sort(key=sort_func)
    wells_coord_reverse = wells_coord[:: -1]
    for well_coord_couple in wells_coord_reverse:
        A_well_column_coo = A_full.getcol((well_coord_couple[1] - 1) * N_r + well_coord_couple[0])
        A_well_column = A_well_column_coo.toarray()
        for cell_number in range(len(A_well_column)):
            if A_well_column[cell_number] != 0:
                B[cell_number] = B[cell_number] - A_well_column[cell_number] * CP_dict[well_coord_couple]

        A_full = A_full.tocsr()
        all_cols = np.arange(A_full.shape[1])
        cols_to_keep = \
        np.where(np.logical_not(np.in1d(all_cols, [(well_coord_couple[1] - 1) * N_r + well_coord_couple[0]])))[0]
        A_full = A_full[:, cols_to_keep]
        A_full = A_full[cols_to_keep, :]
        B = np.delete(B, (well_coord_couple[1] - 1) * N_r + well_coord_couple[0], axis=0)
    P_new = spsolve(A_full, B)

    #print(np.shape(A_full), np.shape(B))
    for well_coord_couple in wells_coord:
        P_new = np.insert(P_new, (well_coord_couple[1]-1)*N_r + well_coord_couple[0], CP_dict[well_coord_couple])
    #print(N_r, M_fi, N_r_oil, np.shape(P_new))
    P_new = P_new.reshape(N_r*M_fi, 1)
    Pres_end = np.zeros((N_r, M_fi))
    j = 0
    for m in range(M_fi):
        for n in range(N_r):
            Pres_end[n][m] = P_new[j][0]
            j += 1

    return Pres_end, A_full, B

if __name__ == '__main__':
    PinPoint = []
    for t in range(T_exp):
        print(t)
        Pres_distrib, A_full, B = PorePressure_in_Time(N_r, M_fi, Pres_distrib, c1, c2, c3_oil, c3_water, c4, wells_coord, CP_dict, delta_r, P_center, N_r_oil)
        PinPoint.append(Pres_distrib[0, int(M_fi/4)]+q*mu_oil*delta_r/perm)
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))

    fig1 = plt.figure()
    plt.plot(Pres_distrib[:,int(M_fi/2)])
    P_all_center = np.ones((1, M_fi)) * P_center
    Pres_distrib = np.vstack((P_all_center, Pres_distrib))
    plt.plot(Pres_distrib[:,30])
    plt.show()
    plt.plot(PinPoint)
    plt.show()
    X = np.zeros((N_r+1,M_fi))
    Y = np.zeros((N_r+1, M_fi))
    for m in range(M_fi):
        for n in range(N_r+1):
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

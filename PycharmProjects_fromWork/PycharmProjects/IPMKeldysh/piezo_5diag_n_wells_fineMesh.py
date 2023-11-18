# Решение уравнения пьезопроводности в цилиндрических координатах, двумерный случай (неявная схема)
# В центре условие непротикания, на границах тоже: градиент давления равен нулю.
# модель однофазного течения, когда во вспомогательную скважину (0.057, 0.127) в декартовых, в цилиндрич: (0.1392, 1.1489)
#  подается постоянное давление до установления
# стац. режима. Потом поток перекрывается. Смотрятся кривые падения давления.
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection



perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2*10**(-3)
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_oil*fi*(Cf+Cr)/perm
frac_angle = 1.1489
delta_r = 0.005
delta_r_fine = 0.005
R_for_fine = 0.04
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine/delta_r_fine)
N_r = round((R-r_well-R_for_fine)/delta_r)
delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*N_r
N_r_full = len(delta_r_list)

delta_fi = np.pi / 30 # угол-шаг в радианах
delta_fi_fine = np.pi/30
fi_for_fine = np.pi/18
M_fi_fine = round(fi_for_fine / delta_fi_fine)
print(int((frac_angle-fi_for_fine)/delta_fi))
print(int(M_fi_fine*2))
print(int((np.pi - 2 * fi_for_fine) / delta_fi))
print(int((np.pi - frac_angle - fi_for_fine)/delta_fi))
delta_fi_list_first = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((np.pi - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((np.pi - frac_angle - fi_for_fine)/delta_fi))
angle_lack = round((2*np.pi - sum(delta_fi_list_first))/delta_fi)
delta_fi_list = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((np.pi - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((np.pi - frac_angle - fi_for_fine)/delta_fi)+angle_lack)
M_fi_full = len(delta_fi_list)
print(sum(delta_fi_list))
print((2*np.pi - sum(delta_fi_list))/delta_fi)


# x_f = 0.02
# N_r_oil = round(x_f/delta_r_fine)


delta_t = 1
Pres = 1*10**5
# P_center = 15*10**5
Pres_distrib = np.ones((N_r_full, M_fi_full)) * Pres
# c1 = 1/delta_r**2
# c2 = 1/2/delta_r
c3_oil = k_oil/delta_t
c3_water = k_water/delta_t
#c4 = 1/delta_fi**2
T_exp = 10
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r**2)/100
wells_coord_real = [(0.1392, 1.1489), (0.1392, np.pi+1.1489)]
#wells_coord = [(int(0.15/delta_r), int(np.pi/4/delta_fi)+5), (int(0.15/delta_r), int(5*np.pi/4/delta_fi)+5), (int(0.17/delta_r), int(3*np.pi/4/delta_fi)+5)]
wells_angles = []
for well_coord in wells_coord_real:
    for angle_number in range(len(delta_fi_list)):
        if well_coord[1] < sum(delta_fi_list[0:angle_number]):
            wells_angles.append(angle_number)
            break
wells_dists = []
for well_coord in wells_coord_real:
    wells_dists.append(N_r_fine + int((well_coord[0]-r_well-R_for_fine)/delta_r)+2)

wells_coord = list(zip( wells_dists, wells_angles))
print(well_coord)

print(wells_coord)
P_well = [1450000, 100000]

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

#for i in range(len(wells_coord)):
#    Pres_distrib[wells_coord[i][0]][wells_coord[i][1]] = P_well[i]

def PorePressure_in_Time(N_r, M_fi, Pres_distrib, c3_water, wells_coord, CP_dict, delta_r):

    # пластовое давление во всей области на нулевом временном шаге
    print(N_r_full, M_fi_full)
    B = np.zeros((N_r_full*M_fi_full, 1))
    for m in range(M_fi_full):

        A = np.zeros((N_r_full, N_r_full))

        for n in range(1, N_r - 1):
            c1 = 1 / delta_r_list[n] ** 2
            c2 = 1 / 2 / delta_r_list[n]
            A[n][n - 1] = c1 - c2 / (sum(delta_r_list[0:n+1]))
            A[n][n] = -2 * c1 - c3_water - 2/(sum(delta_r_list[0:n+1]))**2/delta_fi_list[m]**2
            A[n][n + 1] = c1 + c2 / (sum(delta_r_list[0:n+1]))

        # for n in range(N_r_oil, N_r - 1):
        #     c1 = 1 / delta_r_list[n] ** 2
        #     c2 = 1 / 2 / delta_r_list[n]
        #     A[n][n - 1] = c1 - c2 / (sum(delta_r_list[0:n + 1]))
        #     A[n][n] = -2 * c1 - c3_water - 2 / (sum(delta_r_list[0:n + 1])) ** 2 / delta_fi_list[m] ** 2
        #     A[n][n + 1] = c1 + c2 / (sum(delta_r_list[0:n + 1]))

        # A[N_r_oil - 1][N_r_oil - 1] = 1 / delta_r_list[N_r_oil-1] / mu_oil + 1 / delta_r_list[N_r_oil] / mu_water
        # A[N_r_oil - 1][N_r_oil] = -1 / delta_r_list[N_r_oil] / mu_water
        # A[N_r_oil - 1][N_r_oil - 2] = -1 / delta_r_list[N_r_oil-1] / mu_oil

        c1 = 1 / delta_r_list[0] ** 2
        c2 = 1 / 2 / delta_r_list[0]
        A[0][0] = -2 * c1 - c3_water - 2/(delta_r_list[0])**2/delta_fi_list[m]**2 + (c1 - c2/(delta_r_list[0]))
        A[0][1] = c1 + c2 / (1 * delta_r_list[0])

        c1 = 1 / delta_r_list[N_r_full - 1] ** 2
        c2 = 1 / 2 / delta_r_list[N_r_full - 1]
        A[N_r_full - 1][N_r_full - 1] = -2 * c1 - c3_water + c1 + c2 / (sum(delta_r_list[0:N_r_full])) - 2/(sum(delta_r_list[0:N_r_full]))**2/delta_fi_list[m]**2
        A[N_r_full - 1][N_r_full - 2] = c1 - c2 / (sum(delta_r_list[0:N_r_full]))

        c4 = 1 / delta_fi_list[m] ** 2
        A_sym = np.zeros((N_r_full, N_r_full))
        for n in range(0,N_r_full):
            A_sym[n][n] = c4/(sum(delta_r_list[0:n + 1]))**2

        if m == 0:
            A_line_1 = np.hstack((A, A_sym, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym))
            A_full = A_line_1.copy()
        elif m == M_fi_full-1:
            A_line_end = np.hstack((A_sym, np.zeros((N_r_full, N_r_full * M_fi_full - 3 * N_r_full)), A_sym, A))
            A_full = np.vstack((A_full, A_line_end))
        else:
            A_line = np.hstack((np.zeros((N_r_full,N_r_full*(m-1))), A_sym, A, A_sym, np.zeros((N_r_full, N_r_full * M_fi_full - (3+(m-1)) * N_r_full))))
            A_full = np.vstack((A_full, A_line))



    j = 0
    for m in range(M_fi):
        for n in range(N_r):
            if n == 0:
                #print(j, n, m)
                c1 = 1 / delta_r_list[n] ** 2
                c2 = 1 / 2 / delta_r_list[n]
                B[j][0] = -c3_water * Pres_distrib[n][m]
            else:
                B[j][0] = -c3_water * Pres_distrib[n][m]
            j += 1

    wells_coord.sort(key=sortByAngle)
    wells_coord_reverse = wells_coord[:: -1]
    for well_coord_couple in wells_coord_reverse:
        A_well_column = A_full[:][(well_coord_couple[1]-1)*N_r + well_coord_couple[0]]
        for cell_number in range(len(A_well_column)):
            if A_well_column[cell_number] != 0:
                B[cell_number] = B[cell_number] - A_well_column[cell_number]*CP_dict[well_coord_couple]
        A_full = np.delete(A_full, (well_coord_couple[1]-1)*N_r + well_coord_couple[0], axis=0)
        A_full = np.delete(A_full,
                           (well_coord_couple[1] - 1) * N_r + well_coord_couple[0],axis=1)
        B = np.delete(B, (well_coord_couple[1] - 1) * N_r + well_coord_couple[0], axis=0)

    #print(np.shape(A_full), np.shape(B))
    P_new = np.linalg.solve(A_full, B)
    for well_coord_couple in wells_coord:
        P_new = np.insert(P_new, (well_coord_couple[1]-1)*N_r + well_coord_couple[0], CP_dict[well_coord_couple])
    #print(N_r, M_fi, N_r_oil, np.shape(P_new))
    P_new = P_new.reshape(N_r*M_fi, 1)
    Pres_end = np.zeros((N_r, M_fi))
    j = 0
    for m in range(M_fi):
        for n in range(N_r):
            print(P_new[j][0])
            Pres_end[n][m] = P_new[j][0]
            j += 1

    return Pres_end, A_full, B

if __name__ == '__main__':
    PinPoint = []
    for t in range(T_exp):
        print(t)
        Pres_distrib, A_full, B = PorePressure_in_Time(N_r_full, M_fi_full, Pres_distrib, c3_water, wells_coord, CP_dict, delta_r)
        PinPoint.append(Pres_distrib[2][2])
        print(min(Pres_distrib.flat), max(Pres_distrib.flat))

    fig1 = plt.figure()
    plt.plot(Pres_distrib[:, wells_coord[0][1]])
    plt.plot(Pres_distrib[:, wells_coord[1][1]])
    # P_all_center = np.ones((1, M_fi_full)) * P_center
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
            X[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.cos(sum(delta_fi_list[0:m]))
            Y[n][m] = (r_well+sum(delta_r_list[0:n+1]))*np.sin(sum(delta_fi_list[0:m]))

    print(np.shape(X), np.shape(Y), np.shape(Pres_distrib))

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
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)

    # patches = []
    # for i in range(len(delta_r_list)):
    #     circle = Wedge((0,0), r_well+sum(delta_r_list[0:i]), 0, 360, width=0.001)
    #     patches.append(circle)
    for i in range(len(delta_fi_list)):
        x, y = np.array([[0, (R-r_well)*np.cos(sum(delta_fi_list[0:i]))], [0, (R-r_well)*np.sin(sum(delta_fi_list[0:i]))]])
        line = mlines.Line2D(x, y, lw=0.5)
        ax.add_line(line)
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

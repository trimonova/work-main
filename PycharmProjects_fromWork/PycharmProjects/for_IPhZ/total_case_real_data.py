from flow_in_frac_simple_case import Pressure_in_frac
from piezo_for_IPhZ_sparseMatrix import PorePressure_in_Time
from calculationK1c import calcK1c
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm


perm = 0.5 * 10 ** (-15)  # м2 проницаемость
mu = 2.45 * 10 ** (-3)  # Па*с вязкость
fi = 0.1  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k = mu*fi*(Cf+Cr)/perm
hx = 4
hy = 4
t_step = 60
T_exp = 10
Tswitch = 5
Lx = 1000
Ly = 1000
K1c = 112380

Courant_number = t_step/k/hx**2 + t_step/k/hy**2
print(Courant_number)

N = int(Lx/hx) # количество ячеек вдоль оси х
M = int(Ly/hy)
wells_with_Q = {}

wells_coord = []
P_well = []
print(N, M)

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

Pres = 272*10**5 # давление в пласте
Pbound = 272*10**5 #  давление на границе\
Pres_distrib = np.ones((N, M)) * Pres
# P_add_hor = np.ones((N, 1)) * Pres
# P_add_vert = np.ones((1, M + 2)) * Pres


l_fr0 = 16 # половина длины трещины
N_fr = int(l_fr0/hx)
nu = 0.22
H = 56.7
E = 27*10**9
G = E/2/(1+nu)
k = 4*(1-nu)*H/3.14/G

alpha = 1/12/mu/k
#Pinj = 50*10**5
Sh = 4.95*10**7
#w0 = k*(Pinj - Sh) # 2*10**(-4)
#coef = -perm/mu*(Pinj/2-Pres)/hy/2
#q = np.ones((N_fr-1, 1))*coef
q = np.zeros((N_fr-1, 1))
Qinj = 0.0000007 # м3/с скорость закачки жидкости в скважину с ГРП
w_perf = 0.03 # ширина перфорации
S = H*w_perf
qinj = Qinj/S
w = np.ones((N_fr - 1, 1))*k*(Sh - Sh)
bc_left = k*qinj*mu*hx/perm
w0 = bc_left

# frac_pressure = [Sh for i in range(N_fr)]

l_fr_list = [l_fr0]

X = np.zeros((N,M))
Y = np.zeros((N, M))
for m in range(M):
    for n in range(N):
        X[n][m] = n*hx
        Y[n][m] = m*hy

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]

for t in range(T_exp):
    if t < Tswitch:
        for iter in range(30):
            P_new, w_new = Pressure_in_frac(N_fr, t_step, alpha, w0, q, w, k, Sh, hx, bc_left) # находим давление в трещине
            w = w_new
    else:
        for iter in range(1):
            P_new, w_new = Pressure_in_frac(N_fr, t_step, alpha, w0, q, w, k, Sh, hx, bc_left=0) # находим давление в трещине
            w = w_new
    print(P_new)
    P_new_list = list(P_new)
    P_new_list_reverse = P_new_list[::-1]
    frac_pressure = P_new_list_reverse + P_new_list

    #if t == 50:
    fig = plt.figure()
    surf = plt.plot(frac_pressure)
    plt.show()

    frac_coords = [(int(round(N / 2) - len(frac_pressure)/2 + i), int(round(M / 2))) for i in range(len(frac_pressure))]

    wells_frac_coords = wells_coord + frac_coords
    for i in range(len(frac_coords)):
        CP_dict[frac_coords[i]] = frac_pressure[i]

    Pres_distrib = PorePressure_in_Time(Pres_distrib, N, M, hx, hy, k, t_step, Pres, wells_frac_coords, CP_dict) # распределение давление везде

    P_list = [k for k in Pres_distrib.flat]

    CP_list = zip(X_list, Y_list, P_list)

    print(min(P_list), max(P_list))
    #if t == 50:
    xi = np.linspace(min(X_list), max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list, Y_list), P_list, (xig, yig), method='cubic')

    levels = list(range(Pres, 70000000, 1000000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),
                        linewidth=0.2, levels=levels)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

    print(Pres_distrib[:][round(M/2)+1])
    q = [100*perm/mu/hy * (frac_pressure[i]-Pres_distrib[round(N/2)-len(frac_pressure)/2+i][round(M/2)+1]) for i in range(len(frac_pressure))] # считаем утечки

    K1 = calcK1c(l_fr0, hx, N_fr-1, P_new, Sh)
    print(K1)
    if K1 > K1c:
        l_fr0 += hx
        N_fr += 1
        w = np.vstack((w, np.array(0)))
        q = np.vstack((q, np.array(0)))

    l_fr_list.append(l_fr0)


fig = plt.figure()
surf = plt.plot(l_fr_list)
plt.show()









from piezo_5diag_n_wells_DirichletBC import PorePressure_in_Time as PorePressure_in_Time_Dirichlet
from piezo_5diag_n_wells_NeumanBC import PorePressure_in_Time as PorePressure_in_Time_Neuman
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib import cm

perm = 2 * 10 ** (-15)  # м2 проницаемость
mu_water = 1 * 10 ** (-3)  # Па*с вязкость
mu_oil = 1*10**(-3)
fi = 0.19  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k_water = mu_water*fi*(Cf+Cr)/perm
k_oil = mu_water*fi*(Cf+Cr)/perm
delta_r = 0.005
delta_fi = np.pi / 30 # угол-шаг в радианах
R = 0.215
r_well = 0.0075
# x_f = 0.05
# N_r_oil = int(x_f/delta_r)
N_r = int((R-r_well)/delta_r)
M_fi = int(2*np.pi/delta_fi)
delta_t = 1
Pres = 1*10**5
# P_center = 5*10**5
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
c3_oil = k_oil/delta_t
c3_water = k_water/delta_t
c4 = 1/delta_fi**2
T_exp_dir = 10
T_exp_neu = 10
Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t/k_water/delta_fi**2 + delta_t/k_water/delta_r**2)/100

wells_coord = [(round(0.1392/delta_r), round(1.1489/delta_fi)), (round(0.1392/delta_r), round((np.pi+1.1489)/delta_fi))]
P_well = [1450000, 100000]
print(wells_coord)

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

PinPoint = []
for t in range(T_exp_dir):
    print(t)
    Pres_distrib, A_full, B = PorePressure_in_Time_Dirichlet(N_r, M_fi, Pres_distrib, c1, c2, c3_water, c4, wells_coord, CP_dict, delta_r, delta_fi)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))
    PinPoint.append(Pres_distrib[2][2])

P_well = ['neuman', 100000]
BC_type = ['neuman', 'dirichlet']

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

for t in range(T_exp_neu):
    print(t)
    Pres_distrib, A_full, B = PorePressure_in_Time_Neuman(N_r, M_fi, Pres_distrib, c1, c2, c3_water, c4, wells_coord, CP_dict, delta_r, delta_fi)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))
    PinPoint.append(Pres_distrib[2][2])

fig2 = plt.figure()
plt.plot(PinPoint)
plt.show()
X = np.zeros((N_r,M_fi))
Y = np.zeros((N_r, M_fi))
for m in range(M_fi):
    for n in range(N_r):
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
levels = list(range(0,int(max(P_list)),10000))
fig = plt.figure()
surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

from piezo_5diag_n_wells_DirichletBC import PorePressure_in_Time as PorePressure_in_Time_Dirichlet
from piezo_5diag_n_wells_NeumanBC import PorePressure_in_Time as PorePressure_in_Time_Neuman
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib import cm
from read_mat_file import P as P_exp, r_fi_pair3 as points_r_fi_pair
from read_mat_file import P_filter_0, P_filter_1, P_filter_2, P_filter_3, P_filter_4, P_filter_5
from read_mat_file import P_filter_6, P_filter_7, P_filter_8, P_filter_9, P_filter_10, P_filter_11, P_filter_12
#print(P_exp)
print(points_r_fi_pair)
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
delta_t_dir = 1
delta_t_neu = 2
Pres = 1*10**5*0
# P_center = 5*10**5
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
#c3_oil = k_oil/delta_t
c3_water_dir = k_water/delta_t_dir
c3_water_neu = k_water/delta_t_neu*75
c4 = 1/delta_fi**2
T_exp_dir = 6
T_exp_neu = 2000
#Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t_dir/k_water/delta_fi**2 + delta_t_dir/k_water/delta_r**2)/100
Courant_number_water_2 = (delta_t_neu/k_water/delta_fi**2 + delta_t_neu/k_water/delta_r**2)/100
wells_coord = [(round(0.1392/delta_r), round(1.1489/delta_fi)), (round(0.1392/delta_r), round((np.pi+1.1489)/delta_fi))]
P_well = [1450000, 100000]
print(wells_coord)

points_coords = [(round(i[0]/1000/delta_r), round(i[1]/delta_fi)) for i in points_r_fi_pair]
print('points_coords', points_coords)

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

one_point = points_coords[1]
one_point2 = points_coords[2]
one_point3 = points_coords[3]
one_point4 = points_coords[4]
one_point5 = points_coords[5]
one_point6 = points_coords[6]
one_point7 = points_coords[7]
one_point8 = points_coords[8]
one_point9 = points_coords[9]
one_point10 = points_coords[10]
one_point11 = points_coords[11]
one_point12 = points_coords[12]
one_point0 = points_coords[0]

PinPoint = []
PinPoint2 = []
PinPoint3 = []
PinPoint4 = []
PinPoint5 = []
PinPoint6 = []
PinPoint7 = []
PinPoint8 = []
PinPoint9 = []
PinPoint10 = []
PinPoint11 = []
PinPoint12 = []
PinPoint0 = []

for t in range(T_exp_dir):
    print(t)
    Pres_distrib, A_full, B = PorePressure_in_Time_Dirichlet(N_r, M_fi, Pres_distrib, c1, c2, c3_water_dir, c4, wells_coord, CP_dict, delta_r, delta_fi)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))
    print(np.shape(Pres_distrib))
    # PinPoint.append(Pres_distrib[int(one_point[0])][int(one_point[1])])
    # PinPoint2.append(Pres_distrib[int(one_point2[0])][int(one_point2[1])])
    # PinPoint3.append(Pres_distrib[int(one_point3[0])][int(one_point3[1])])
    # PinPoint4.append(Pres_distrib[int(one_point4[0])][int(one_point4[1])])
    # PinPoint5.append(Pres_distrib[int(one_point5[0])][int(one_point5[1])])
    # PinPoint6.append(Pres_distrib[int(one_point6[0])][int(one_point6[1])])
    # PinPoint7.append(Pres_distrib[int(one_point7[0])][int(one_point7[1])])
    # PinPoint8.append(Pres_distrib[int(one_point8[0])][int(one_point8[1])])
    # PinPoint9.append(Pres_distrib[int(one_point9[0])][int(one_point9[1])])
    # PinPoint10.append(Pres_distrib[int(one_point10[0])][int(one_point10[1])])
    # PinPoint11.append(Pres_distrib[int(one_point11[0])][int(one_point11[1])])
    # PinPoint12.append(Pres_distrib[int(one_point12[0])][int(one_point12[1])])
    # PinPoint0.append(Pres_distrib[int(one_point0[0])][int(one_point0[1])])

P_well = ['neuman', 130000]
BC_type = ['neuman', 'dirichlet']

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]


for t in range(T_exp_neu):
    print(t)
    Pres_distrib, A_full, B = PorePressure_in_Time_Neuman(N_r, M_fi, Pres_distrib, c1, c2, c3_water_neu, c4, wells_coord, CP_dict, delta_r, delta_fi)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))
    PinPoint.append(Pres_distrib[int(one_point[0])][int(one_point[1])])
    PinPoint2.append(Pres_distrib[int(one_point2[0])][int(one_point2[1])])
    PinPoint3.append(Pres_distrib[int(one_point3[0])][int(one_point3[1])])
    PinPoint4.append(Pres_distrib[int(one_point4[0])][int(one_point4[1])])
    PinPoint5.append(Pres_distrib[int(one_point5[0])][int(one_point5[1])])
    PinPoint6.append(Pres_distrib[int(one_point6[0])][int(one_point6[1])])
    PinPoint7.append(Pres_distrib[int(one_point7[0])][int(one_point7[1])])
    PinPoint8.append(Pres_distrib[int(one_point8[0])][int(one_point8[1])])
    PinPoint9.append(Pres_distrib[int(one_point9[0])][int(one_point9[1])])
    PinPoint10.append(Pres_distrib[int(one_point10[0])][int(one_point10[1])])
    PinPoint11.append(Pres_distrib[int(one_point11[0])][int(one_point11[1])])
    PinPoint12.append(Pres_distrib[int(one_point12[0])][int(one_point12[1])])
    PinPoint0.append(Pres_distrib[int(one_point0[0])][int(one_point0[1])])

# t_list = [9000 + i*delta_t_dir for i in range(T_exp_dir)] + [9000 + T_exp_dir*delta_t_dir + j*delta_t_neu for j in range(T_exp_neu)]
t_list = [9000 + T_exp_dir*delta_t_dir + j*delta_t_neu for j in range(T_exp_neu)]
t_list_exp = [i*0.01 for i in range(np.shape(P_filter_1)[0])]


fig, axs = plt.subplots(3, 2, figsize=(21, 10))
axs[0, 0].set_title(str((round(points_r_fi_pair[1][0], 2), round(points_r_fi_pair[1][1], 2))), y = 0.75, loc='left')
axs[0, 0].set_xlabel("Time, s")
axs[0, 0].set_ylabel("Pressure, MPa")
axs[0, 0].plot(t_list_exp, P_filter_1)
axs[0, 0].plot(t_list, np.array(PinPoint)/10**6)
axs[0, 1].set_title(str((round(points_r_fi_pair[2][0], 2), round(points_r_fi_pair[2][1], 2))), y = 0.75, loc='left')
axs[0, 1].set_xlabel("Time, s")
axs[0, 1].set_ylabel("Pressure, MPa")
axs[0, 1].plot(t_list_exp, P_filter_2)
axs[0, 1].plot(t_list, np.array(PinPoint2)/10**6)
# axs[1, 0].set_title(str((round(points_r_fi_pair[3][0], 2), round(points_r_fi_pair[3][1], 2))), y = 0.75, loc='left')
# axs[1, 0].plot(t_list_exp, P_filter_3*10**6)
# axs[1, 0].plot(t_list, PinPoint3)
# axs[1, 1].set_title(str((round(points_r_fi_pair[4][0], 2), round(points_r_fi_pair[4][1], 2))), y = 0.75, loc='left')
# axs[1, 1].plot(t_list_exp, P_filter_4*10**6)
# axs[1, 1].plot(t_list, PinPoint4)
# axs[2, 0].set_title(str((round(points_r_fi_pair[5][0], 2), round(points_r_fi_pair[5][1], 2))), y = 0.75, loc='left')
# axs[2, 0].plot(t_list_exp, P_filter_5*10**6)
# axs[2, 0].plot(t_list, PinPoint5)
axs[1, 1].set_title(str((round(points_r_fi_pair[6][0], 2), round(points_r_fi_pair[6][1], 2))), y = 0.75, loc='left')
axs[1, 1].set_xlabel("Time, s")
axs[1, 1].set_ylabel("Pressure, MPa")
axs[1, 1].plot(t_list_exp, P_filter_6)
axs[1, 1].plot(t_list, np.array(PinPoint6)/10**6)
axs[1, 0].set_title(str((round(points_r_fi_pair[7][0], 2), round(points_r_fi_pair[7][1], 2))), y = 0.75, loc='left')
axs[1, 0].set_xlabel("Time, s")
axs[1, 0].set_ylabel("Pressure, MPa")
axs[1, 0].plot(t_list_exp, P_filter_7)
axs[1, 0].plot(t_list, np.array(PinPoint7)/10**6)
# axs[1, 2].set_title(str((round(points_r_fi_pair[8][0], 2), round(points_r_fi_pair[8][1], 2))), y = 0.75, loc='left')
# axs[1, 2].plot(t_list_exp, P_filter_8*10**6)
# axs[1, 2].plot(t_list, PinPoint8)
axs[2, 0].set_title(str((round(points_r_fi_pair[9][0], 2), round(points_r_fi_pair[9][1], 2))), y = 0.75, loc='left')
axs[2, 0].set_xlabel("Time, s")
axs[2, 0].set_ylabel("Pressure, MPa")
axs[2, 0].plot(t_list_exp, P_filter_9)
axs[2, 0].plot(t_list, np.array(PinPoint9)/10**6)
# axs[0, 3].set_title(str((round(points_r_fi_pair[10][0], 2), round(points_r_fi_pair[10][1], 2))), y = 0.75, loc='left')
# axs[0, 3].plot(t_list_exp, P_filter_10*10**6)
# axs[0, 3].plot(t_list, PinPoint10)
# axs[1, 3].set_title(str((round(points_r_fi_pair[11][0], 2), round(points_r_fi_pair[11][1], 2))), y = 0.75, loc='left')
# axs[1, 3].plot(t_list_exp, P_filter_11*10**6)
# axs[1, 3].plot(t_list, PinPoint11)
# axs[2, 3].set_title(str((round(points_r_fi_pair[12][0], 2), round(points_r_fi_pair[12][1], 2))), y = 0.75, loc='left')
# axs[2, 3].plot(t_list_exp, P_filter_12*10**6)
# axs[2, 3].plot(t_list, PinPoint12)

# fig1 = plt.figure()
# #plt.plot(PinPoint)
# plt.plot(t_list_exp, P_filter_1*10**6)
#
# plt.plot(t_list, PinPoint)
# #plt.plot(P_filter_0)
# plt.show()
#
# fig2 = plt.figure()
# plt.plot(t_list_exp, P_filter_2*10**6)
#
# plt.plot(t_list, PinPoint2)
# plt.show()
#
# fig3 = plt.figure()
# plt.plot(t_list_exp, P_filter_3*10**6)
#
# plt.plot(t_list, PinPoint3)
# plt.show()

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

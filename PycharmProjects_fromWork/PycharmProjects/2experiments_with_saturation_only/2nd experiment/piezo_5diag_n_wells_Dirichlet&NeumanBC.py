from piezo_5diag_n_wells_DirichletBC import PorePressure_in_Time as PorePressure_in_Time_Dirichlet
from piezo_5diag_n_wells_NeumanBC import PorePressure_in_Time as PorePressure_in_Time_Neuman
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from matplotlib import cm
from read_mat_file import pressure as pressure_exp
from read_mat_file import points_coords_dict_rad

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
N_r = round((R-r_well)/delta_r)
M_fi = int(2*np.pi/delta_fi)
delta_t_dir = 10
delta_t_neu = 10
Pres = 1*10**5
# P_center = 5*10**5
Pres_distrib = np.ones((N_r, M_fi)) * Pres
c1 = 1/delta_r**2
c2 = 1/2/delta_r
#c3_oil = k_oil/delta_t
c3_water_dir = k_water/delta_t_dir*25
c3_water_neu = k_water/delta_t_neu*25
c4 = 1/delta_fi**2
T_exp_dir = 480
T_exp_neu = 200
#Courant_number_oil = (delta_t/k_oil/delta_fi**2 + delta_t/k_oil/delta_r**2)/100
Courant_number_water = (delta_t_dir/k_water/delta_fi**2 + delta_t_dir/k_water/delta_r**2)/100
Courant_number_water_2 = (delta_t_neu/k_water/delta_fi**2 + delta_t_neu/k_water/delta_r**2)/100
wells_coord = [(round(0.171/delta_r), round((np.pi+np.pi/4)/delta_fi)), (round(0.171/delta_r), round((np.pi/4)/delta_fi))]
P_well = [450000, 140000]
print(wells_coord)

points_coord = {} # координаты датчиков в шагах сетки
for elem in points_coords_dict_rad:
    points_coord[elem] = [round(points_coords_dict_rad[elem][0]/delta_r), round(points_coords_dict_rad[elem][1]/delta_fi)]

print('points_coords', points_coord)

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

def sortByRad(inputSet):
    return inputSet[0]

def sortByAngle(inputSet):
    return inputSet[1]

column_in_pressure_6 = 6
column_in_pressure_3 = 3
column_in_pressure_4 = 4
column_in_pressure_8 = 8
column_in_pressure_11 = 11
column_in_pressure_13 = 13
column_in_pressure_15 = 15
column_in_pressure_17 = 17
column_in_pressure_19 = 19
column_in_pressure_21 = 21
column_in_pressure_23 = 23
one_point_6 = points_coord[column_in_pressure_6] # координаты датчкика
one_point_3 = points_coord[column_in_pressure_3]
one_point_4 = points_coord[column_in_pressure_4]
one_point_8 = points_coord[column_in_pressure_8]
one_point_11 = points_coord[column_in_pressure_11]
one_point_13 = points_coord[column_in_pressure_13]
one_point_15 = points_coord[column_in_pressure_15]
one_point_17 = points_coord[column_in_pressure_17]
one_point_19 = points_coord[column_in_pressure_19]
one_point_21 = points_coord[column_in_pressure_21]
one_point_23 = points_coord[column_in_pressure_23]

PinPoint_6 = []
PinPoint_3 = []
PinPoint_4 = []
PinPoint_8 = []
PinPoint_11 = []
PinPoint_13 = []
PinPoint_15 = []
PinPoint_17 = []
PinPoint_19 = []
PinPoint_21 = []
PinPoint_23 = []

for t in range(T_exp_dir):
    print(t)
    Pres_distrib, A_full, B = PorePressure_in_Time_Dirichlet(N_r, M_fi, Pres_distrib, c1, c2, c3_water_dir, c4, wells_coord, CP_dict, delta_r, delta_fi, r_well)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))
    print(np.shape(Pres_distrib))
    PinPoint_6.append(Pres_distrib[int(one_point_6[0])][int(one_point_6[1])])
    PinPoint_3.append(Pres_distrib[int(one_point_3[0])][int(one_point_3[1])])
    PinPoint_4.append(Pres_distrib[int(one_point_4[0])][int(one_point_4[1])])
    PinPoint_8.append(Pres_distrib[int(one_point_8[0])][int(one_point_8[1])])
    PinPoint_11.append(Pres_distrib[int(one_point_11[0])][int(one_point_11[1])])
    PinPoint_13.append(Pres_distrib[int(one_point_13[0])][int(one_point_13[1])])
    PinPoint_15.append(Pres_distrib[int(one_point_15[0])][int(one_point_15[1])])
    PinPoint_17.append(Pres_distrib[int(one_point_17[0])][int(one_point_17[1])])
    PinPoint_19.append(Pres_distrib[int(one_point_19[0])][int(one_point_19[1])])
    PinPoint_21.append(Pres_distrib[int(one_point_21[0])][int(one_point_21[1])])
    PinPoint_23.append(Pres_distrib[int(one_point_23[0])][int(one_point_23[1])])

P_well = ['neuman', 130000]
BC_type = ['neuman', 'dirichlet']

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]


for t in range(T_exp_neu):
    print(t)
    Pres_distrib, A_full, B = PorePressure_in_Time_Neuman(N_r, M_fi, Pres_distrib, c1, c2, c3_water_neu, c4, wells_coord, CP_dict, delta_r, delta_fi, r_well)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))
    PinPoint_6.append(Pres_distrib[int(one_point_6[0])][int(one_point_6[1])])
    PinPoint_3.append(Pres_distrib[int(one_point_3[0])][int(one_point_3[1])])
    PinPoint_4.append(Pres_distrib[int(one_point_4[0])][int(one_point_4[1])])
    PinPoint_8.append(Pres_distrib[int(one_point_8[0])][int(one_point_8[1])])
    PinPoint_11.append(Pres_distrib[int(one_point_11[0])][int(one_point_11[1])])
    PinPoint_13.append(Pres_distrib[int(one_point_13[0])][int(one_point_13[1])])
    PinPoint_15.append(Pres_distrib[int(one_point_15[0])][int(one_point_15[1])])
    PinPoint_17.append(Pres_distrib[int(one_point_17[0])][int(one_point_17[1])])
    PinPoint_19.append(Pres_distrib[int(one_point_19[0])][int(one_point_19[1])])
    PinPoint_21.append(Pres_distrib[int(one_point_21[0])][int(one_point_21[1])])
    PinPoint_23.append(Pres_distrib[int(one_point_23[0])][int(one_point_23[1])])

t_list = [i*delta_t_dir for i in range(T_exp_dir)] + [T_exp_dir*delta_t_dir + j*delta_t_neu for j in range(T_exp_neu)]
t_list_exp = [i*0.01 for i in range(np.shape(pressure_exp)[0])]

# fig2 = plt.figure()
# plt.plot(t_list_exp, pressure_exp[:, column_in_pressure_6])
# plt.plot(t_list, PinPoint_6 )
# plt.show()
fig, axs = plt.subplots(4, 2, figsize=(21, 10))
axs[0, 0].set_title(str((round(one_point_6[0]*delta_r*1000, 2), round(one_point_6[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[0, 0].set(xlim=(0, 7000), ylim=(0, 1))
axs[0, 0].set_xlabel("Time, s")
axs[0, 0].set_ylabel("Pressure, MPa")
axs[0, 0].plot(t_list_exp, pressure_exp[:, column_in_pressure_6]/10**6)
axs[0, 0].plot(t_list, np.array(PinPoint_6)/10**6)
# axs[0, 1].set_title(str((round(one_point_3[0], 2), round(one_point_3[1], 2))), y = 0.75, loc='left')
# axs[0, 1].plot(t_list_exp, pressure_exp[:, column_in_pressure_3])
# axs[0, 1].plot(t_list, PinPoint_3)
# axs[0, 2].set_title(str((round(one_point_4[0], 2), round(one_point_4[1], 2))), y = 0.75, loc='left')
# axs[0, 2].plot(t_list_exp, pressure_exp[:, column_in_pressure_4])
# axs[0, 2].plot(t_list, PinPoint_4)
axs[0, 1].set_title(str((round(one_point_8[0]*delta_r*1000, 2), round(one_point_8[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[0, 1].set(xlim=(0, 7000), ylim=(0, 1))
axs[0, 1].set_xlabel("Time, s")
axs[0, 1].set_ylabel("Pressure, MPa")
axs[0, 1].plot(t_list_exp, pressure_exp[:, column_in_pressure_8]/10**6)
axs[0, 1].plot(t_list, np.array(PinPoint_8)/10**6)
axs[1, 0].set_title(str((round(one_point_11[0]*delta_r*1000, 2), round(one_point_11[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[1, 0].set(xlim=(0, 7000), ylim=(0, 1))
axs[1, 0].set_xlabel("Time, s")
axs[1, 0].set_ylabel("Pressure, MPa")
axs[1, 0].plot(t_list_exp, pressure_exp[:, column_in_pressure_11]/10**6)
axs[1, 0].plot(t_list, np.array(PinPoint_11)/10**6)
axs[1, 1].set_title(str((round(one_point_13[0]*delta_r*1000, 2), round(one_point_13[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[1, 1].set(xlim=(0, 7000), ylim=(0, 1))
axs[1, 1].set_xlabel("Time, s")
axs[1, 1].set_ylabel("Pressure, MPa")
axs[1, 1].plot(t_list_exp, pressure_exp[:, column_in_pressure_13]/10**6)
axs[1, 1].plot(t_list, np.array(PinPoint_13)/10**6)
axs[2, 0].set_title(str((round(one_point_15[0]*delta_r, 2), round(one_point_15[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[2, 0].set(xlim=(0, 7000), ylim=(0, 1))
axs[2, 0].set_xlabel("Time, s")
axs[2, 0].set_ylabel("Pressure, MPa")
axs[2, 0].plot(t_list_exp, pressure_exp[:, column_in_pressure_15]/10**6)
axs[2, 0].plot(t_list, np.array(PinPoint_15)/10**6)
axs[2, 1].set_title(str((round(one_point_17[0]*delta_r, 2), round(one_point_17[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[2, 1].set(xlim=(0, 7000), ylim=(0, 1))
axs[2, 1].set_xlabel("Time, s")
axs[2, 1].set_ylabel("Pressure, MPa")
axs[2, 1].plot(t_list_exp, pressure_exp[:, column_in_pressure_17]/10**6)
axs[2, 1].plot(t_list, np.array(PinPoint_17)/10**6)
# axs[2, 0].set_title(str((round(one_point_19[0], 2), round(one_point_19[1], 2))), y = 0.75, loc='left')
# axs[2, 0].plot(t_list_exp, pressure_exp[:, column_in_pressure_19])
# axs[2, 0].plot(t_list, PinPoint_19)
axs[3, 0].set_title(str((round(one_point_21[0]*delta_r, 2), round(one_point_21[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[3, 0].set(xlim=(0, 7000), ylim=(0, 1))
axs[3, 0].set_xlabel("Time, s")
axs[3, 0].set_ylabel("Pressure, MPa")
axs[3, 0].plot(t_list_exp, pressure_exp[:, column_in_pressure_21]/10**6)
axs[3, 0].plot(t_list, np.array(PinPoint_21)/10**6)
axs[3, 1].set_title(str((round(one_point_23[0]*delta_r, 2), round(one_point_23[1]*delta_fi, 2))), y = 0.75, loc='left')
axs[3, 1].set(xlim=(0, 7000), ylim=(0, 1))
axs[3, 1].set_xlabel("Time, s")
axs[3, 1].set_ylabel("Pressure, MPa")
axs[3, 1].plot(t_list_exp, pressure_exp[:, column_in_pressure_23]/10**6)
axs[3, 1].plot(t_list, np.array(PinPoint_23)/10**6)

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

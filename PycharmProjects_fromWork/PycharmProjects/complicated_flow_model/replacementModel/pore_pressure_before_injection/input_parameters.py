import numpy as np

perm = 2 * 10 ** (-12)  # м2 проницаемость
mu_water = 2 * 10 ** (-3)  # Па*с вязкость
mu_oil = 2 * 10 ** (-3)
porosity = 0.4  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
C_total = (Cf + Cr) * 25
k_water = mu_water * porosity * C_total / perm
k_oil = mu_oil * porosity * C_total / perm
frac_angle = np.pi/4
frac_angle_2 = np.pi/4*5
frac_length_1 = 0.01
frac_length_2 = 0.01
delta_r = 0.0001
delta_r_fine = 0.0001
R_for_fine = 0.02
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine / delta_r_fine)
delta_r_list = [delta_r_fine] * N_r_fine + [delta_r] * round((R - r_well - R_for_fine) / delta_r)
N_r_full = len(delta_r_list)

delta_fi = np.pi / 45  # угол-шаг в радианах
delta_fi_fine = np.pi / 45
fi_for_fine = np.pi / 6
M_fi_fine = round(fi_for_fine / delta_fi_fine)

delta_fi_list_first = [delta_fi] * round((frac_angle - fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine * 2) + [
    delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine * 2) + [
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

delta_t = 0.5
Pres = 0 * 10 ** 5
#Q_center = 0.2 * 10 ** (-8)  # из лаб. данных - 0.2*10**(-6) m3/s
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

wells_coord = [(round(0.1711/delta_r), round(np.pi/4/delta_fi)), (round(0.1711/delta_r), round((np.pi+np.pi/4)/delta_fi))]
P_well = [1000000, 0]
#P_well = []
#wells_coord = []
print(wells_coord)

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

viscosity_matrix = np.ones((N_r_full, M_fi_full))*mu_water

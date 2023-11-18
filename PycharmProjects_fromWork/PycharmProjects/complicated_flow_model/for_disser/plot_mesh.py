import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection

frac_angle = np.pi/6
frac_angle_2 = np.pi/6*7

delta_r = 0.01
delta_r_fine = 0.005
R_for_fine = 0.05
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine/delta_r_fine)
delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*round((R-r_well-R_for_fine)/delta_r)

delta_fi = np.pi / 18 # угол-шаг в радианах
delta_fi_list = [delta_fi]*round(2*np.pi/delta_fi)
# delta_fi_fine = np.pi/18
# fi_for_fine = np.pi/6
# M_fi_fine = round(fi_for_fine / delta_fi_fine)
#
# delta_fi_list_first = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi))
# angle_lack = round((2*np.pi - sum(delta_fi_list_first))/delta_fi)
# delta_fi_list = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi)+angle_lack)


fig, ax = plt.subplots()
patches = []
for i in range(len(delta_r_list)):
    circle = Wedge((0,0), r_well+sum(delta_r_list[0:i]), 0, 360, width=0.001)
    patches.append(circle)
for i in range(len(delta_fi_list)):
    x, y = np.array([[0, (R-r_well)*np.cos(sum(delta_fi_list[0:i]))], [0, (R-r_well)*np.sin(sum(delta_fi_list[0:i]))]])
    line = mlines.Line2D(x, y, lw=0.5)
    ax.add_line(line)

plt.xlim(-0.25, 0.25)
plt.ylim(-0.25, 0.25)
plt.xlabel("X, м")
plt.ylabel("Y, м")

# t = np.arange(0, 2 * np.pi, 0.01)
# r = 0.215
# plt.plot(r * np.sin(t), r * np.cos(t))
# ax = fig.gca(projection='3d')

p = PatchCollection(patches)
ax.add_collection(p)
plt.show()
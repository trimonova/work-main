from scipy.signal import savgol_filter
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

f = h5py.File('Давление/data2.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

x = pressure_dict['xp']/1000
y = pressure_dict['yp']/1000
pt = pressure_dict['p']*10**6

t = pressure_dict['t']
fig = plt.figure()
#plt.axis([15000,38000, 0, 15000000])
#time = [i/100 for i in range(8000,20000)]
#plt.plot(pt[10])
print(np.shape(pt[10]))
p_drop = pt[10][30820:-1]
p_drop_filter = savgol_filter(p_drop, 201, 3)
print(np.shape(p_drop_filter))
plt.plot(p_drop_filter)
plt.show()

# print(np.shape(y))
# print(np.shape(pt))
#
# x_list = []
# y_list = []
#
# for elem in x[0]:
#     x_list.append(elem)
#
# for elem in y[0]:
#     y_list.append(elem)
#
# set_coord = list(zip(x_list, y_list))
# print(set_coord)
#
# dist_from_0 = []
# for couple in set_coord:
#     dist_from_0.append((couple[0]**2 + couple[1]**2)**0.5)
#
# x = []
# y = []
# i = 0
#
# fig = plt.figure()
# plt.axis([15000,38000, 0, 15000000])
# time = [i/100 for i in range(8000,20000)]
# plt.plot(pt[10])
#
#
# plt.xlabel('Time, *0.01 s')
# plt.ylabel('Pressure, Pa')
# plt.show()


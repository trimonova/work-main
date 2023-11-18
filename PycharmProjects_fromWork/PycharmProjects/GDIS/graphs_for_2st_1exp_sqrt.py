from scipy.signal import savgol_filter
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import diff

f = h5py.File('Давление/data1.mat', 'r')

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
#plt.plot(pt[9])
print(np.shape(pt[9]))
p_drop = pt[9][8768:-1]


tp = 0
p_drop_filter = savgol_filter(p_drop, 3001, 3)

dP = abs(np.diff(p_drop_filter))

time = np.array([0.01+0.01*i for i in range(np.shape(p_drop_filter)[0])])
sqrt_time = np.array([np.sqrt(i) for i in time])

d_sqrt_time = np.diff(sqrt_time)
dp_dsqrt_t = dP/d_sqrt_time
dP_func = sqrt_time[0:-1]*dp_dsqrt_t

plt.plot(sqrt_time, p_drop_filter)
plt.plot( sqrt_time[0:-1], dp_dsqrt_t)
plt.plot( sqrt_time[0:-1], dP_func, 'r')
plt.show()

#plt.plot(time, log_time)


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


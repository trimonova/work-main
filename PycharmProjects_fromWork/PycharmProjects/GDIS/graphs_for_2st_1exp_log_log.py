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
plt.plot(pt[9])
plt.show()
print(np.shape(pt[9]))
p_drop = pt[9][8768:-1]
p_drop_filter = savgol_filter(p_drop, 3001, 3)
p_drop_filter_diff = np.diff(p_drop_filter)
print(type(p_drop_filter_diff))
print(np.shape(p_drop))
print(np.shape(p_drop_filter))
print(np.shape(p_drop_filter_diff))
print(np.shape(p_drop_filter_diff)[0])
#plt.plot(p_drop_filter)
plt.plot(p_drop_filter)
plt.show()

tp = 0
time = np.array([0.01+0.01*i for i in range(np.shape(p_drop_filter)[0])])
time_diff = np.diff(time)
p_drop_true_diff = p_drop_filter_diff/time_diff
log_time = np.array([np.log10(i) for i in time])
under_diff = time[0:-1]*p_drop_true_diff
print(time[0], p_drop_filter_diff[0], under_diff[0])
log_y = np.array([np.log10(abs(i)) for i in under_diff])
#plt.plot(time, log_time)

log_y2 = np.array([np.log10(i) for i in p_drop_filter])

print(np.shape(log_y), np.shape(log_time))
plt.plot(log_time[0:-1], log_y)
plt.plot(log_time, log_y2)
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


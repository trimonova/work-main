import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import scipy.io
f = scipy.io.loadmat('flowPlusflow_smooth.mat')
#f = h5py.File('flowPlusflow_smooth.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value
    #print(elem)
pressure = pressure_dict['flow_smooth']*10**6
plt.plot(pressure[:, 6])
plt.show()

points_coords_dict = {3:[0, 0], 4:[-121, -121], 6:[-127,0], 8: [127, 0], 11: [57,	-127], 13: [-57, 127], 15: [0, 70],
17: [0, -185], 19: [65,	65], 21: [-121, 121], 23: [57, 127]}
points_coords_dict_rad = {}
for elem in points_coords_dict:
    x = points_coords_dict[elem][0]/1000
    y = points_coords_dict[elem][1]/1000
    if np.arctan2(y, x) >= 0:
        points_coords_dict_rad[elem] = [(x**2+y**2)**0.5, np.arctan2(y,x)]
    else:
        points_coords_dict_rad[elem] = [(x ** 2 + y ** 2) ** 0.5, 2*np.pi + np.arctan2(y, x)]
# pressure = pressure_dict['flow_smooth'].transpose()*10**6
# P_filter_6 = savgol_filter(pressure[:,6], polyorder=3, window_length=151)
# plt.plot(pressure[:, 6])
# plt.plot(P_filter_6)
# plt.axis([0, 700000, 0, 4*10**5])
# plt.show()
#
# P_filter_13 = savgol_filter(pressure[:,13], polyorder=3, window_length=5051)
# plt.plot(pressure[:, 13])
# plt.plot(P_filter_13)
# plt.axis([0, 700000, 0, 4*10**5])
# plt.show()
#
# P_filter_11 = savgol_filter(pressure[:,11], polyorder=3, window_length=151)
# plt.plot(pressure[:, 11])
# #plt.plot(P_filter_11)
# plt.axis([0, 700000, 0, 4*10**5])
# plt.show()
#
# P_filter_21 = savgol_filter(pressure[:,21], polyorder=3, window_length=151)
# plt.plot(pressure[:, 21])
# #plt.plot(P_filter_21)
# plt.axis([0, 700000, 0, 4*10**5])
# plt.show()
#
# #print(pressure_dict['flow'])
# #print(points_coords_dict_rad)
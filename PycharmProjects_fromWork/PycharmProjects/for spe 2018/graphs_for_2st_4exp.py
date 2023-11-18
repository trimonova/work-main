from scipy import interpolate
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd

f = h5py.File('Давление\data4.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

x = pressure_dict['xp']/1000
y = pressure_dict['yp']/1000
pt = pressure_dict['p']
xhole = pressure_dict['xhole']/1000
yhole = pressure_dict['yhole']/1000
Pini = pressure_dict['Pini'].transpose()*10**6
t = pressure_dict['t']
print(np.shape(y))
print(np.shape(pt))

x_list = []
y_list = []

for elem in x[0]:
    x_list.append(elem)

for elem in y[0]:
    y_list.append(elem)

set_coord = list(zip(x_list, y_list))
print(set_coord)

dist_from_0 = []
for couple in set_coord:
    dist_from_0.append((couple[0]**2 + couple[1]**2)**0.5)

x = []
y = []
i = 0

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
#plt.axis([4500, 20000, 0, 15000000])

for dist in dist_from_0:
    if i != 11:
        ax2.plot(t, pt[i])
        i += 1
    else:
        ax1.plot(t, pt[i])
        i += 1
ax1.set_xlabel('Время, секунды')
ax1.set_ylabel('Давление в центральной скважине, МПа')
ax2.set_ylabel('Давление в образце, МПа')
ax2.axis([40, 160, 0, 3])


plt.show()
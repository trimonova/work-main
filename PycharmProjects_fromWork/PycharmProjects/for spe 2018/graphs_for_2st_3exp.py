from scipy import interpolate
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

f = h5py.File('Давление\data3.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value

x = pressure_dict['xp']/1000
y = pressure_dict['yp']/1000
pt = pressure_dict['p']

t = pressure_dict['t']*0.01

x_list = []
y_list = []

for elem in x[0]:
    x_list.append(elem)

for elem in y[0]:
    y_list.append(elem)

set_coord = list(zip(x_list, y_list))


dist_from_0 = []
for couple in set_coord:
    dist_from_0.append((couple[0]**2 + couple[1]**2)**0.5)

x = []
y = []
i = 0
if __name__ == '__main__':
    fig = plt.figure()
    #plt.axis([8000, 25000, 0, 15000000])
    time = [i/100 for i in range(8000,20000)]
    plt.plot(t, pt.transpose())

    print(np.shape(y))
    print(np.shape(pt))
    print(set_coord)
    plt.xlabel('Time, *0.01 s')
    plt.ylabel('Pressure, Pa')
    plt.show()


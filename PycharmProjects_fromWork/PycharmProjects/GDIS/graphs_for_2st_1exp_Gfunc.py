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



p_drop_filter = savgol_filter(p_drop, 3001, 3)

dP = abs(np.diff(p_drop_filter))
tp = 20
time = np.array([20 + 0.01+0.01*i for i in range(np.shape(p_drop_filter)[0])])
delta_time = (time - tp)/tp
g_time = 4/3*((1+delta_time)**1.5 - delta_time**1.5)
G_time = 4/np.pi*(g_time - 4/3)

G_time_diff = np.diff(G_time)

dP_dG = dP/G_time_diff
G_dP_dG = G_time[0:-1]*dP_dG

print(G_time)
#plt.plot(G_time, p_drop_filter)
plt.plot(G_time[0:-1], dP_dG)
plt.plot(G_time[0:-1], G_dP_dG)
plt.show()

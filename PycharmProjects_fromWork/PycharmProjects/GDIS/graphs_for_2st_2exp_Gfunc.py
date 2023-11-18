from scipy.signal import savgol_filter
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import diff

f = h5py.File('pressure/data2.mat', 'r')

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
plt.plot(pt[10])

print(np.shape(pt[10]))
p_drop = pt[10][30820:-1]

p_drop_filter = savgol_filter(p_drop, 301, 3)
plt.plot(p_drop_filter)
plt.show()
dP = abs(np.diff(p_drop_filter))
tp = 308.2
time = np.array([0.01+0.01*i for i in range(np.shape(p_drop_filter)[0])])
delta_time = time/tp
g_time = 4/3*((1+delta_time)**1.5 - delta_time**1.5)
G_time = 4/np.pi*(g_time - 4/3)

G_time_diff = np.diff(G_time)

dP_dG = dP/G_time_diff
G_dP_dG = G_time[0:-1]*dP_dG

dP_dG_filter = savgol_filter(dP_dG, 301, 3)
G_dP_dG_filter = savgol_filter(G_dP_dG, 301, 3)

print(G_time)
#plt.plot(G_time, p_drop_filter)
fig, ax1 = plt.subplots()
ax1.set_ylabel('pressure', color='red')
ax1.set_xlabel('G-func')
ax1.plot(G_time[0:-1], p_drop_filter[0:-1], color='red')
ax2 = ax1.twinx()
ax2.set_ylabel('derivatives', color='blue')
ax2.plot(G_time[0:-1], dP_dG_filter, color='green')
ax2.plot(G_time[0:-1], G_dP_dG_filter, color='blue')
plt.show()

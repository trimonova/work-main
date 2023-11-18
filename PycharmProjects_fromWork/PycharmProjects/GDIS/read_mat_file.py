import h5py
import numpy as np
import matplotlib.pyplot as plt
f = h5py.File('pressure for RN/Pressure Exp8-1.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())
print(keys_list)
#print(f['Pini'])
for elem in keys_list:
    if elem != '#subsystem#':
        k_value = np.array(f.get(elem))
        pressure_dict[elem] = k_value
        print(elem)



#x = pressure_dict['xp']/1000
#y = pressure_dict['yp']/1000
P = pressure_dict['P'].transpose()
plt.plot(P[:, 0])
plt.show()
print(np.shape(P))
print(pressure_dict['xp'])
print(pressure_dict['yp'])

for j in range(np.shape(P)[1]):
    for i in range(np.shape(P)[0]-1):
        if P[i][j] > 0.8:
            k = 1
            while P[i+k][j] > 0.8:
                k = k+1
            P[i][j] = P[i+k][j]
        if P[i][j] < 0:
            P[i][j] = 0
        if abs(P[i+1][j]-P[i][j]) > 0.15:
            P[i + 1][j] = P[i][j]

P = np.delete(P, [0], 1)
print(np.shape(P))
print(np.shape(P[:, 1]))
# fig = plt.figure()
# plt.plot(P[:,2])
# plt.show()

from scipy.signal import savgol_filter

P_filter_0 = savgol_filter(P[:,0], polyorder=3, window_length=151)
P_filter_1 = savgol_filter(P[:,1], polyorder=3, window_length=151) # window size 51, polynomial order 3
P_filter_2 = savgol_filter(P[:,2], polyorder=3, window_length=151)
P_filter_3 = savgol_filter(P[:,3], polyorder=3, window_length=151)
P_filter_4 = savgol_filter(P[:,4], polyorder=3, window_length=151)
P_filter_5 = savgol_filter(P[:,5], polyorder=3, window_length=151) # window size 51, polynomial order 3
P_filter_6 = savgol_filter(P[:,6], polyorder=3, window_length=151)
P_filter_7 = savgol_filter(P[:,7], polyorder=3, window_length=151)
P_filter_8 = savgol_filter(P[:,8], polyorder=3, window_length=151)
P_filter_9 = savgol_filter(P[:,9], polyorder=3, window_length=151) # window size 51, polynomial order 3
P_filter_10 = savgol_filter(P[:,10], polyorder=3, window_length=151)
P_filter_11 = savgol_filter(P[:,11], polyorder=3, window_length=151)
P_filter_12 = savgol_filter(P[:,12], polyorder=3, window_length=151)

print(np.shape(P_filter_0))

# fig = plt.figure()
# plt.plot(P[:, 2])
# plt.plot(P_filter_2)
# plt.show()


xp = pressure_dict['xp'][0]
yp = pressure_dict['yp'][0]
# r_fi_pair = [((xp[i]**2+yp[i]**2)**0.5, np.arctan(yp[i]/xp[i])) for i in range(len(xp))]
# r_fi_pair2 = [((xp[i]**2+yp[i]**2)**0.5, np.arctan2(yp[i],xp[i])) for i in range(len(xp))]

r_fi_pair3 = []
for i in range(len(xp)):
    if np.arctan2(yp[i],xp[i]) >= 0:
        r_fi_pair3.append(((xp[i]**2+yp[i]**2)**0.5, np.arctan2(yp[i],xp[i])))
    else:
        r_fi_pair3.append(((xp[i] ** 2 + yp[i] ** 2) ** 0.5, 2*np.pi + np.arctan2(yp[i], xp[i])))

# print(r_fi_pair)
# print(r_fi_pair2)
# print(r_fi_pair3)
from scipy.signal import savgol_filter
import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('flow.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())
pressure = np.array(f.get('flow'))
print(np.shape(pressure))

p_one = pressure[0, :]
fig = plt.figure()
plt.plot(p_one)

p_filter = savgol_filter(p_one, 3001, 5)
print(np.shape(p_filter))
plt.plot(p_filter, 'r')
plt.show()

# for i in range(25):
#     plt.plot(pressure[i,:])
#     #plt.axis([0, 620000, 0, 0.1])
#     plt.show()

#
# fig = plt.figure()
# plt.plot(pt)

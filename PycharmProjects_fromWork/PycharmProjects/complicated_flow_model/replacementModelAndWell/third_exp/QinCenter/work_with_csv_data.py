import openpyxl
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import pyplot as plt
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter


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
plt.plot(pt[5])
plt.show()
print(np.shape(pt[10]))
p_drop = pt[10][30820:-1]
p_drop_filter = savgol_filter(pt[5], 3001, 3)
# p_drop_filter_diff = np.diff(p_drop_filter)
# print(type(p_drop_filter_diff))
# print(np.shape(p_drop))
# print(np.shape(p_drop_filter))
# print(np.shape(p_drop_filter_diff))
# print(np.shape(p_drop_filter_diff)[0])
#plt.plot(p_drop_filter)
plt.plot(p_drop_filter)
plt.show()

df = pd.read_csv('C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects/complicated_flow_model/replacementModelAndWell/dataframes_QinCenter/Q0_2_10_-8/dataframe_visc_t_10.csv')
print(df['0'].to_list())
print(df.columns)
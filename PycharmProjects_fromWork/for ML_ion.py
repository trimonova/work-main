import numpy as np
import h5py
#f = h5py.File('C:\disk_D\ионосферные штуки\Data_for_Masha_2017_09_06.mat', 'r')
#from mat4py import loadmat

#data = loadmat('C:\disk_D\ионосферные штуки\Data_for_Masha_2017_09_06.mat')
#data = f.get()
#data = np.array(data)
import numpy as np
import scipy.io
import pandas as pd
f = scipy.io.loadmat('C:\disk_D\ионосферные штуки\Data_for_Masha_2017_09_06.mat')
wl = f['wl']
M = f['M']
el = f['el']
I = f['I']
print(np.shape(wl), np.shape(M), np.shape(el), np.shape(I))

for station_number in range(np.shape(el)[1]):
    I_temp = I * el[0][station_number]
    M_temp = M[:, station_number].reshape(8,1)
    #print(np.shape(M_temp))
    if station_number == 0:
        I_total = I_temp
        M_total = M_temp
    else:
        I_total = np.vstack((I_total, I_temp))
        M_total = np.vstack((M_total, M_temp))
print(np.shape(I_total), np.shape(M_total))
total_data = np.hstack((M_total, I_total))
print(np.shape(total_data))
frame = pd.DataFrame(total_data)
frame.to_csv('C:\disk_D\ионосферные штуки\Data_for_Masha_2017_09_06.csv',index=False)


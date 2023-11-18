import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib as mpl
import matplotlib.pyplot as plt

f = h5py.File('pressure/p.mat', 'r')

pressure_dict = {}
keys_list = list(f.keys())

for elem in keys_list:
    k_value = np.array(f.get(elem))
    pressure_dict[elem] = k_value
    print(elem)



#x = pressure_dict['xp']/1000
#y = pressure_dict['yp']/1000
pt = pressure_dict['p'].transpose()*10**6
fig = plt.figure()
plt.plot(pt)
plt.show()
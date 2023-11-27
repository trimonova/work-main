import numpy as np
import matplotlib.pyplot as plt
from input_parameters import N_r_full, M_fi_full, r_well, delta_r, delta_fi

X = np.zeros((N_r_full, M_fi_full))
Y = np.zeros((N_r_full, M_fi_full))
for m in range(M_fi_full):
    for n in range(N_r_full):
        X[n][m] = (r_well + (n + 1) * delta_r) * np.cos(delta_fi * m)
        Y[n][m] = (r_well + (n + 1) * delta_r) * np.sin(delta_fi * m)

Pres_distrib_in_time = np.load('Pres_distrib_before_fracturing.npy')
surf = plt.contourf(X, Y, Pres_distrib_in_time[-1], linewidth=0.2, cmap=plt.get_cmap('jet'))
plt.title('pressure')
plt.colorbar(surf)
plt.show()

for time_step in range(len(Pres_distrib_in_time)):
    plt.scatter(time_step, Pres_distrib_in_time[time_step][int(N_r_full/3), int(M_fi_full/3)])
plt.show()
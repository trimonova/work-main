import os
import numpy as np
from matplotlib import pyplot as plt

directory = '/Users/trimonovds/PycharmProjects/pythonProject3.9/PycharmProjects_fromWork/PycharmProjects/gazogidraty/different 2D cases for article'

t_list = [100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000, 500000, 600000, 650000, 700000, 750000, 800000, 850000, 900000, 950000, 1000000]
t_list_in_days = [elem/60/24 for elem in t_list]

min_value_list = []
max_value_list = []
for entry in os.scandir(directory):
    if entry.is_file() and entry.name.endswith('.npy'):
        #print(entry.path)
        #print(entry.name)
        min_value = np.min(np.load(entry.name))
        min_value_list.append(min_value/10**6)

        max_value = np.max(np.load(entry.name))
        max_value_list.append(max_value/10**6)
        #print(min_value)

min_value_list.sort()
max_value_list.sort()

for element in min_value_list:
    print(element)

L_x = 100
L_y = 100
delta_x = 5
delta_y = 5

N_r_full = round(L_x/delta_x)

M_fi_full = round(L_y/delta_y)
Pres = 1*10**5
ro = 1700 # 1750
g = 9.8
P_lit_distrib = np.zeros((N_r_full, M_fi_full))
for n in range(N_r_full):
    for m in range(M_fi_full):
        P_lit_distrib[n][-1-m] = (m+1)*delta_y*g*ro + Pres

P_lit_min = np.min(P_lit_distrib)
P_lit_max = np.max(P_lit_distrib)

P_lit_max_list = [P_lit_max/10**6]*len(t_list)
P_lit_min_list = [P_lit_min/10**6]*len(t_list)

#plt.plot(t_list, min_value_list)
plt.plot(t_list_in_days, max_value_list)
#plt.plot(t_list_in_days, P_lit_max_list)
plt.xlabel('Time, days')
plt.ylabel('Max pressure, MPa')
#plt.show()
plt.savefig('figures/MaxPressure_vs_time_fromArticle.png', dpi=1000, bbox_inches='tight', pad_inches=0.1)

plt.figure()
plt.plot(t_list_in_days, min_value_list)
#plt.plot(t_list_in_days, P_lit_min_list)
plt.xlabel('Time, days')
plt.ylabel('Min pressure, MPa')
#plt.show()
plt.savefig('figures/MinPressure_vs_time_fromArticle.png', dpi=1000, bbox_inches='tight', pad_inches=0.1)






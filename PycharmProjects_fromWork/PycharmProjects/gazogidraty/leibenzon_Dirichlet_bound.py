import numpy as np
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import csv
import pandas as pd


def PorePressure_in_Time(N_r_full, Pres_distrib, delta_r, mu, fi, perm, delta_t, P_bound_start, P_bound_end):

    # пластовое давление во всей области на нулевом временном шаге
    alpha = perm/2/fi/mu
    B = np.zeros((N_r_full, 1))
    A = np.zeros((N_r_full, N_r_full))

    for n in range(1, N_r_full - 1):

        c1 = 2*alpha*Pres_distrib[n]/delta_r**2
        c2 = 1/delta_t + 4*alpha*Pres_distrib[n]/delta_r**2
        A[n][n - 1] = c1

        A[n][n] = -c2

        A[n][n + 1] = c1

        B[n] = -Pres_distrib[n]/delta_t - 2*alpha*((Pres_distrib[n+1]-Pres_distrib[n-1])/2/delta_r)**2

    c1 = 2 * alpha * Pres_distrib[0] / delta_r ** 2
    c2 = 1 / delta_t + 4 * alpha * Pres_distrib[0] / delta_r ** 2
    A[0][0] = -c2
    A[0][1] = c1
    B[0] = -c1*P_bound_start - Pres_distrib[0]/delta_t - 2*alpha*((Pres_distrib[1]-P_bound_start)/2/delta_r)**2

    c1 = 2 * alpha * Pres_distrib[N_r_full-1] / delta_r ** 2
    c2 = 1 / delta_t + 4 * alpha * Pres_distrib[N_r_full-1] / delta_r ** 2
    A[N_r_full - 1][N_r_full - 1] = -c2
    A[N_r_full - 1][N_r_full - 2] = c1
    B[N_r_full-1] = -c1*P_bound_end - Pres_distrib[N_r_full-1]/delta_t - 2*alpha*((P_bound_end-Pres_distrib[N_r_full-2])/2/delta_r)**2

    P_new = spsolve(A, B)

    return P_new, A, B

mu_0 = 1028*10**(-8)
C = 198
T = 280
#mu = mu_0*(273+C)/(T+C)*(T/273)**1.5
mu = 1.8*10**(-5)
l = 0.2
delta_r = 0.005
N_r_full = round(l/delta_r)
P_res = 1.013*10**(5)
Pres_distrib = B = np.ones((N_r_full, 1))*P_res
fi = 0.4
perm = 3*10**(-10)
delta_t = 1
P_bound_end = 1.013*10**(5)
P_bound_start = 1.0476*10**(5)
P_bound_start_1 = 1.013*10**(5)
P_bound_start_2 = 1.0476*10**(5)
P_bound_start_2_exp = 1.013*10**(5) + 0.02335*10**5

P_bound_distrib = []
Q = 10**(-5)
S = 1
T_exp = 150
T_1 = 64

for i in range(T_exp):
    P_bound_distrib.append((P_bound_start_2_exp - P_bound_start_1)/T_1*i + P_bound_start_1)

P_bound_distrib_MPa = [i/1000000 for i in P_bound_distrib]
my_df = pd.DataFrame(P_bound_distrib_MPa)
my_df.to_csv('bound_pressure.csv', index=False, header=False)

plt.plot(P_bound_distrib)
plt.title("Зависимость давления на нижней границе от времени")
plt.show()

P_bound_distrib.append(1.0476*10**(5))
P_bound_distrib.append(1.0476*10**(5))


alpha = perm/2/fi/mu
delta_t_max = perm*delta_r**2/2/alpha/P_res
print(delta_t_max)
ro = 1750
g = 9.8
P_lit = [i*delta_r*g*ro + P_res for i in range(N_r_full+2)]
P_lit_reverse = P_lit[: : -1]
P_lit_der = [g*ro for i in range(N_r_full+2)]
plt.plot(P_lit)
plt.show()

coords_list = [delta_r*i for i in range(N_r_full+2)]
print(coords_list)
pressure_time_depend = [coords_list]
pressure_time_depend.append(P_lit)
P_in_point = []
for t in range(105):
    P_new, A, B = PorePressure_in_Time(N_r_full, Pres_distrib, delta_r, mu, fi, perm, delta_t, P_bound_distrib[t], P_bound_end)
    Pres_distrib = P_new
    P_total = [P_bound_distrib[t]] + list(P_new) + [P_bound_end]
    P_in_point.append(P_new[0])
    der_P = [-(P_total[i+1]-P_total[i])/delta_r for i in range(N_r_full+1)]
    for i in range(N_r_full+2):
        if P_total[i] > P_lit_reverse[i]:
            print('critical P', i, t)
            # plt.plot(P_total)
            # plt.plot(P_lit_reverse)
            # plt.show()
            break
        # if der_P[i] > P_lit_der[i]:
        #     print('critical derP', i, t)
        #     plt.plot(der_P)
        #     plt.plot(P_lit_der)
        #     plt.show()
        #     break
    if t == 10:
        pressure_time_depend.append(P_total[: : -1])
    if t == 20:
        pressure_time_depend.append(P_total[: : -1])
    if t == 30:
        pressure_time_depend.append(P_total[: : -1])
    if t == 40:
        pressure_time_depend.append(P_total[: : -1])
    if t == 50:
        pressure_time_depend.append(P_total[: : -1])
    if t == 64:
        pressure_time_depend.append(P_total[: : -1])
    if t == 70:
        pressure_time_depend.append(P_total[: : -1])
    if t == 80:
        pressure_time_depend.append(P_total[: : -1])
    if t == 90:
        pressure_time_depend.append(P_total[: : -1])
    if t == 95:
        pressure_time_depend.append(P_total[: : -1])
    if t == 100:
        pressure_time_depend.append(P_total[: : -1])

    plt.plot(coords_list, P_total[: : -1])


plt.plot(coords_list, P_lit, color='r')
my_df = pd.DataFrame(pressure_time_depend)
my_df.to_csv('my_csv_3.csv', index=False, header=False)

    #plt.show()
# plt.plot(P_in_point)
plt.show()
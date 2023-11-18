from calculationK1c import calcK1c
from graphs_for_2st_3exp import pt
from flow_in_frac_simple_case import Pressure_in_frac
import numpy as np
import matplotlib.pyplot as plt
l_fr0 = 0.05
lag = 0.01
hx = 0.005
t_step = 0.5
N_fr = int((l_fr0-lag)/hx)
nu = 0.2
H = 0.07
E = 3*10**9
G = E/2/(1+nu)
k = 4*(1-nu)*H/3.14/G
mu = 0.1
perm = 2*10**(-15)
alpha = 1/12/mu/k
Pinj = 50*10**5
Sh = 5*10**5
Pres = 1*10**5
#coef = perm/mu*(Pinj/2-Pres)/hx/2 #-5*10**(-6)
#coef = 5*10**(-5)
coef = 0
q = np.ones((N_fr-1, 1))*coef
#q = np.zeros((N_fr-1, 1))
K1c = 21019374
T_exp = 20
w = np.zeros((N_fr - 1, 1))
l_fr_list = [l_fr0]
for t in range(T_exp):
    if t == 0:
        for iter in range(9):
            w0 = k * (pt[11][11000+50*t] - Sh)
            P_new, w_new = Pressure_in_frac(N_fr, t_step, alpha, w0, q, w, k, Sh, hx)
            w = w_new
    else:
        w0 = k * (pt[11][11000 + 50 * t] - Sh)
        P_new, w_new = Pressure_in_frac(N_fr, t_step, alpha, w0, q, w, k, Sh, hx)
        w = w_new
    K1 = calcK1c(l_fr0, hx, N_fr-1, P_new)
    # print(P_new)
    # print(w_new)
    fig = plt.figure()
    surf = plt.plot(P_new)
    plt.show()

    if K1 > K1c:
        l_fr0 += hx
        N_fr += 1
        w = np.vstack((w, np.array(0)))
        q = np.vstack((q, np.array(0)))

    l_fr_list.append(l_fr0)

fig = plt.figure()
surf = plt.plot(l_fr_list)
plt.show()

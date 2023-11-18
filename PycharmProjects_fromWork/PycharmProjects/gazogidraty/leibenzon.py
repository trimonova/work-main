import numpy as np
from scipy.sparse import coo_matrix, linalg, hstack, vstack, csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def PorePressure_in_Time(N_r_full, Pres_distrib, delta_r, mu, fi, perm, delta_t, P_bound, Q, S):

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
    Q_coef = Q*mu*delta_r/S/perm
    A[0][0] = -c2 + c1
    A[0][1] = c1
    B[0] = -Pres_distrib[0]/delta_t - 2*alpha*((Pres_distrib[1]-Pres_distrib[0]-Q_coef)/2/delta_r)**2 - c1*Q_coef

    c1 = 2 * alpha * Pres_distrib[N_r_full-1] / delta_r ** 2
    c2 = 1 / delta_t + 4 * alpha * Pres_distrib[N_r_full-1] / delta_r ** 2
    A[N_r_full - 1][N_r_full - 1] = -c2
    A[N_r_full - 1][N_r_full - 2] = c1
    B[N_r_full-1] = -c1*P_bound - Pres_distrib[N_r_full-1]/delta_t - 2*alpha*((P_bound-Pres_distrib[N_r_full-2])/2/delta_r)**2

    P_new = spsolve(A, B)

    return P_new, A, B

mu_0 = 1028*10**(-8)
C = 198
T = 280
mu = mu_0*(273+C)/(T+C)*(T/273)**1.5
l = 5
delta_r = 0.05
N_r_full = round(l/delta_r)
P_res = 10**(5)
Pres_distrib = B = np.ones((N_r_full, 1))*P_res
fi = 0.2
perm = 7.5*10**(-15)
delta_t = 20
P_bound = 10**(5)
Q = 10**(-5)
S = 1
T_exp = 2500
alpha = perm/2/fi/mu
delta_t_max = perm*delta_r**2/2/alpha/P_res
print(delta_t_max)
ro = 1800
g = 9.8
P_lit = [i*delta_r*g*ro + P_res for i in range(N_r_full+2)]
P_lit_reverse = P_lit[: : -1]
P_lit_der = [g*ro for i in range(N_r_full+1)]
plt.plot(P_lit)
plt.show()

P_in_point = []
for t in range(T_exp):
    P_new, A, B = PorePressure_in_Time(N_r_full, Pres_distrib, delta_r, mu, fi, perm, delta_t, P_bound, Q, S)
    Pres_distrib = P_new
    P_total = [Q*mu*delta_r/S/perm + P_new[0]] + list(P_new) + [P_bound]
    P_in_point.append(P_new[0])
    der_P = [-(P_total[i+1]-P_total[i])/delta_r for i in range(N_r_full+1)]
    plt.plot(P_total)
    for i in range(N_r_full+2):
        if P_total[i] > P_lit_reverse[i]:
            print('critical P', i, t)
            plt.plot(P_total)
            plt.plot(P_lit_reverse)
            plt.title('Pressure')
            plt.show()

        # if der_P[i] > P_lit_der[i]:
        #     print('critical derP', i, t)
        #     plt.plot(der_P)
        #     plt.plot(P_lit_der)
        #     plt.title('Pressure derivative')
        #     plt.show()
        #     break
    #plt.plot(P_total)
    #plt.plot(P_lit)

    #plt.show()
#plt.plot(P_in_point)
plt.plot(P_lit_reverse)
plt.show()
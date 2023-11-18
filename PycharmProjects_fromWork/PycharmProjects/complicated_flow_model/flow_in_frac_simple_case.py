# рассчитывается давление в трещине по модели PKN с заданным давлением на скважине ( в центре трещины) и с миним. напряжением на кончике трещины.
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    l_fr0 = 0.07
    hx = 0.005
    t_step = 1
    N_fr = int(l_fr0/hx)
    nu = 0.2
    H = 0.07
    E = 3*10**9
    G = E/2/(1+nu)
    k = 4*(1-nu)*H/3.14/G
    mu = 0.1
    perm = 2*10**(-15)
    alpha = 1/12/mu/k
    Pinj = 6.48*10**6
    Sh = 4.52*10**6
    Pres = 1*10**5
    w0 = k*(Pinj - Sh) # 2*10**(-4)
    #coef = -perm/mu*(Pinj/2-Pres)/hy/2
    #q = np.ones((N_fr-1, 1))*coef
    q = np.zeros((N_fr-1, 1))

    T_exp = 100
    w = np.ones((N_fr - 1, 1))*k*(Sh - Sh)

def Pressure_in_frac(N_fr, t_step, alpha, w0, q, w, k, Sh):

    A = np.zeros((N_fr - 1, N_fr - 1))
    B = np.zeros((N_fr - 1, 1))
    for n in range(1, N_fr-2):
        w_right3 = ((w[n+1]+w[n])/2)**3
        w_left3 = ((w[n]+w[n-1])/2)**3
        A[n][n] = 1/t_step + alpha/hx*(w_right3/hx + w_left3/hx)
        A[n][n-1] = -alpha/hx*w_left3/hx
        A[n][n+1] = -alpha/hx*w_right3/hx

    w_right3_0 = ((w[0 + 1] + w[0]) / 2) ** 3
    w_left3_0 = ((w[0] + w0) / 2) ** 3
    A[0][0] = 1/t_step + alpha/hx*(w_right3_0/hx + w_left3_0/hx)
    A[0][1] = -alpha/hx*w_right3/hx

    w_right3_end = ((0 + w[N_fr-2]) / 2) ** 3
    w_left3_end = ((w[N_fr-2] + w[N_fr - 3]) / 2) ** 3
    A[N_fr-2][N_fr-2] = 1/t_step + alpha/hx*(w_right3_end/hx + w_left3_end/hx)
    A[N_fr-2][N_fr-3] = -alpha/hx*w_left3_end/hx

    for n in range(0, N_fr-1):
        B[n] = 1/t_step*w[n] - q[n]

    w_left3_another = ((w[0] + w0) / 2) ** 3
    B[0] = B[0] - (-alpha/hx*w_left3_another/hx)*w0

    w_new = np.linalg.solve(A,B)

    w = w_new.reshape(N_fr-1, 1)
    P_new = w / k + Sh

    return P_new, w



if __name__ == '__main__':

    for t in range(T_exp):
        P_new, w_new = Pressure_in_frac(N_fr, t_step, alpha, w0, q, w, k, Sh)
        w = w_new

    print(P_new)
    print(w_new)
    fig = plt.figure()
    surf = plt.plot(P_new)
    plt.show()



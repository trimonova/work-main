# рассчитывается давление в трещине по модели PKN с заданным давлением на скважине ( в центре трещины) и с миним. напряжением на кончике трещины.
import numpy as np
import matplotlib.pyplot as plt

delta_r = 0.005
delta_r_fine = 0.001
R_for_fine = 0.015
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine/delta_r_fine)
delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*round((R-r_well-R_for_fine)/delta_r)
N_r_full = len(delta_r_list)

if __name__ == '__main__':
    l_fr0 = 0.1
    hx = 0.01
    hy = 0.01
    t_step = 0.1
    N_fr = int(l_fr0/hx)
    N_well = int(N_fr/2-1)+2
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
    w0 = k*(Pinj - Sh) # 2*10**(-4)
    #coef = -perm/mu*(Pinj/2-Pres)/hy/2
    #q = np.ones((N_fr-1, 1))*coef
    q = np.zeros((N_fr-1, 1))

    T_exp = 50
    w = np.ones((N_fr - 1, 1))*k*(Sh - Sh)

def Pressure_in_frac(N_fr, t_step, N_well, alpha, w0, q, w, k, Sh):

    A = np.zeros((N_fr - 1, N_fr - 1))
    B = np.zeros((N_fr - 1, 1))
    for n in range(1, N_fr-2):
        w_right3 = ((w[n+1]+w[n])/2)**3
        w_left3 = ((w[n]+w[n-1])/2)**3
        delta_r_right = delta_r_list[n+1]
        delta_r_left = delta_r_list[n]
        A[n][n] = 1/t_step + alpha/(delta_r_right/2 + delta_r_left/2)*(w_right3/delta_r_right + w_left3/delta_r_left)
        A[n][n-1] = -alpha/(delta_r_right/2 + delta_r_left/2)*w_left3
        A[n][n+1] = -alpha/(delta_r_right/2 + delta_r_left/2)*w_right3

    w_right3_0 = ((w[0 + 1] + w[0]) / 2) ** 3
    w_left3_0 = ((w[0] + w0) / 2) ** 3
    delta_r_right_0 = delta_r_list[0 + 1]
    delta_r_left_0 = delta_r_list[0]
    A[0][0] = 1/t_step + alpha/(delta_r_right_0/2 + delta_r_left_0/2)*(w_right3_0/delta_r_right_0 + w_left3_0/delta_r_left_0)
    A[0][1] = -alpha/(delta_r_right_0/2 + delta_r_left_0/2)*w_right3

    w_right3_end = ((0 + w[N_fr-2]) / 2) ** 3
    w_left3_end = ((w[N_fr-2] + w[N_fr - 3]) / 2) ** 3
    delta_r_right_end = delta_r_list[N_fr-1]
    delta_r_left_end = delta_r_list[N_fr-2]
    A[N_fr-2][N_fr-2] = 1/t_step + alpha/(delta_r_right_end/2 + delta_r_left_end/2)*(w_right3_end/delta_r_right_end + w_left3_end/delta_r_left_end)
    A[N_fr-2][N_fr-3] = -alpha/(delta_r_right_end/2 + delta_r_left_end/2)*w_left3_end

    for n in range(0, N_fr-1):
        B[n] = 1/t_step*w[n] - q[n]

    w_left3_another = ((w[0] + w0) / 2) ** 3
    delta_r_right_another = delta_r_list[0 + 1]
    delta_r_left_another = delta_r_list[0]
    B[0] = B[0] - (-alpha/(delta_r_right_another/2 + delta_r_left_another/2)*w_left3_another)*w0
    # for n in range(0, N_fr-1):
    #     if n+1 == N_well:
    #         A[n][n + 1] = 0
    #         B[n] = B[n] - alpha * (w[n] ** 3 + w[n+1]**3)*w0
    #     if n-1 == N_well:
    #         A[n][n - 1] = 0
    #         B[n] = B[n] - alpha * (w[n] ** 3 + w[n-1]**3)*w0
    #
    # A = np.delete(A, N_well, axis=0)
    # A = np.delete(A, N_well, axis=1)
    # B = np.delete(B, N_well)

    w_new = np.linalg.solve(A,B)
    # w_new = np.insert(w_new, N_well, w0)


    w = w_new.reshape(N_fr-1, 1)
    P_new = w / k + Sh

    return P_new, w



if __name__ == '__main__':

    for t in range(T_exp):
        P_new, w_new = Pressure_in_frac(N_fr, t_step, N_well, alpha, w0, q, w, k, Sh)
        w = w_new

    print(P_new)
    print(w_new)
    fig = plt.figure()
    surf = plt.plot(P_new)
    plt.show()



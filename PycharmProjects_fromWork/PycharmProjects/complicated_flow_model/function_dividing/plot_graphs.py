import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
def plot_graphs(delta_r_list, delta_fi_list, r_well, funcToPlot):
    N_r_full = len(delta_r_list)
    M_fi_full = len(delta_fi_list)
    X = np.zeros((N_r_full, M_fi_full))
    Y = np.zeros((N_r_full, M_fi_full))
    for m in range(M_fi_full):
        for n in range(N_r_full):
            X[n][m] = (r_well + sum(delta_r_list[0:n + 1])) * np.cos(sum(delta_fi_list[0:m]))
            Y[n][m] = (r_well + sum(delta_r_list[0:n + 1])) * np.sin(sum(delta_fi_list[0:m]))

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]
    funcToPlot_list = [k for k in funcToPlot.flat]

    xi = np.linspace(min(X_list), max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    funcToPlot_i = interpolate.griddata((X_list, Y_list), funcToPlot_list, (xig, yig), method='cubic')

    fig = plt.figure()
    surf = plt.contourf(xig, yig, funcToPlot_i, linewidth=0.2, cmap=plt.get_cmap('jet'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


def plot_scatter(bound_coord):
    set_bound_coord = set(bound_coord)
    plt.figure()
    for coord_pair in set_bound_coord:
        plt.scatter(coord_pair[0], coord_pair[1])
    plt.show()
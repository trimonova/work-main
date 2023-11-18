import numpy as np
import shapely.geometry as geom
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm

if __name__ == '__main__':

    delta_r = 0.005
    delta_r_fine = 0.005
    R_for_fine = 0.015
    R = 0.215
    r_well = 0.0075
    N_r_fine = round(R_for_fine/delta_r_fine)

    delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*round((R-r_well-R_for_fine)/delta_r)
    N_r_full = len(delta_r_list)

    delta_fi = np.pi / 60 # угол-шаг в радианах
    delta_fi_fine = np.pi/180
    fi_for_fine = np.pi/18
    M_fi_fine = round(fi_for_fine / delta_fi_fine)
    #delta_fi_list_first = [delta_fi]*round((frac_angle_1-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle_1 - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi))
    #angle_lack = round((2*np.pi - sum(delta_fi_list_first))/delta_fi)
    delta_fi_list = [delta_fi] * round(2 * np.pi / delta_fi)
    M_fi_full = len(delta_fi_list)


def replace_boundary(r_well, delta_fi_list, delta_r_list):
    r_bound = 0.05
    Func_coord_dict = {}

    for j in range(len(delta_fi_list)):
        for i in range(len(delta_r_list)):
            x = (r_well+sum(delta_r_list[0:i])) * np.cos(sum(delta_fi_list[0:j]))
            y = (r_well+sum(delta_r_list[0:i])) * np.sin(sum(delta_fi_list[0:j]))
            if sum(delta_r_list[0:i]) < r_bound:
                Func_coord_dict[(i,j)] = (x,y,1)
            else:
                Func_coord_dict[(i,j)]=(x, y, -1)


    #-------------------------------------------------------------------------------------

    bound_coord = []

    def func(temp):
        return temp[1], temp[0]
    #
    keys_list_sorted = list(sorted(Func_coord_dict.keys(), key=func))
    for i in range(len(keys_list_sorted)-1):
        if Func_coord_dict[keys_list_sorted[i]][2] != Func_coord_dict[keys_list_sorted[i+1]][2] and keys_list_sorted[i][0] != len(delta_r_list)-1:
            bound_coord.append((Func_coord_dict[keys_list_sorted[i]][0], Func_coord_dict[keys_list_sorted[i]][1]))


    line = geom.LineString(bound_coord)
    area = geom.polygon.Polygon(line)

    for i in range(len(delta_r_list)):
        r_const_row = np.zeros((1, len(delta_fi_list)))
        for j in range(len(delta_fi_list)):
            coord_x = Func_coord_dict[(i,j)][0]
            coord_y = Func_coord_dict[(i,j)][1]
            point = geom.Point(coord_x, coord_y)
            min_dist = point.distance(line)
            if area.contains(point):
                r_const_row[0][j] = -min_dist
            else:
                r_const_row[0][j] = min_dist
        if i == 0:
            Func_matrix = r_const_row
        else:
            Func_matrix = np.vstack((Func_matrix, r_const_row))

    return bound_coord, Func_coord_dict, Func_matrix


if __name__ == '__main__':
    bound_coord, Func_coord_dict, Func_matrix = replace_boundary(r_well, delta_fi_list, delta_r_list)
    X_list = []
    Y_list = []
    Func_list = []
    for key in Func_coord_dict:
        X_list.append(Func_coord_dict[key][0])
        Y_list.append(Func_coord_dict[key][1])
        Func_list.append(Func_coord_dict[key][2])


    xi = np.linspace(min(X_list), max(X_list), 500)
    yi = np.linspace(min(Y_list), max(Y_list), 500)
    xig, yig = np.meshgrid(xi, yi)

    coord_x_list = [Func_coord_dict[(i, j)][0] for i in range(len(delta_r_list)) for j in range(len(delta_fi_list))]
    coord_y_list = [Func_coord_dict[(i, j)][1] for i in range(len(delta_r_list)) for j in range(len(delta_fi_list))]

    Func_matrix_i = interpolate.griddata((coord_x_list, coord_y_list), Func_matrix.flat, (xig, yig), method='cubic')
    levels = list(range(0, 5000000, 10000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Func_matrix_i, cmap=cm.jet, antialiased=True, linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    set_bound_coord = set(bound_coord)
    fig1 = plt.figure()
    print(bound_coord)
    for coord_pair in set_bound_coord:
        plt.scatter(coord_pair[0], coord_pair[1])
    plt.show()



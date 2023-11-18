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
    M_fi_full = len(delta_fi_list)
    N_r_full = len(delta_r_list)

    X = np.zeros((N_r_full, M_fi_full))
    Y = np.zeros((N_r_full, M_fi_full))
    for m in range(M_fi_full):
        for n in range(N_r_full):
            X[n][m] = (r_well + sum(delta_r_list[0:n + 1])) * np.cos(sum(delta_fi_list[0:m]))
            Y[n][m] = (r_well + sum(delta_r_list[0:n + 1])) * np.sin(sum(delta_fi_list[0:m]))

    X_list = [i for i in X.flat]
    Y_list = [j for j in Y.flat]

    xi = np.linspace(min(X_list), max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)

    r_bound = 0.005
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

    Func_matrix_list = [l for l in Func_matrix.flat]
    Func_matrix_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
    surf = plt.contour(xig, yig, Func_matrix_i, 0)
    zero_contour = surf.collections[1].get_paths()[0].vertices
    x_zero = zero_contour[:, 0]
    y_zero = zero_contour[:, 1]
    new_bound_coord = [(x_zero[i], y_zero[i]) for i in range(len(x_zero))]

    line_1 = geom.LineString(new_bound_coord)
    area_1 = geom.polygon.Polygon(line_1)

    for i in range(len(delta_r_list)):
        r_const_row = np.zeros((1, len(delta_fi_list)))
        for j in range(len(delta_fi_list)):
            coord_x = Func_coord_dict[(i, j)][0]
            coord_y = Func_coord_dict[(i, j)][1]
            point = geom.Point(coord_x, coord_y)
            min_dist = point.distance(line_1)
            if area_1.contains(point):
                r_const_row[0][j] = -min_dist
            else:
                r_const_row[0][j] = min_dist
        if i == 0:
            Func_matrix_remake = r_const_row
        else:
            Func_matrix_remake = np.vstack((Func_matrix_remake, r_const_row))

    Func_matrix_list = [l for l in Func_matrix_remake.flat]
    Func_matrix_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
    surf = plt.contour(xig, yig, Func_matrix_i, 0)
    zero_contour = surf.collections[1].get_paths()[0].vertices
    x_zero = zero_contour[:, 0]
    y_zero = zero_contour[:, 1]
    bound_coord_remake = [(x_zero[i], y_zero[i]) for i in range(len(x_zero))]

    line_2 = geom.LineString(bound_coord_remake)
    area_2 = geom.polygon.Polygon(line_2)

    for i in range(len(delta_r_list)):
        r_const_row = np.zeros((1, len(delta_fi_list)))
        for j in range(len(delta_fi_list)):
            coord_x = Func_coord_dict[(i, j)][0]
            coord_y = Func_coord_dict[(i, j)][1]
            point = geom.Point(coord_x, coord_y)
            min_dist = point.distance(line_2)
            if area_2.contains(point):
                r_const_row[0][j] = -min_dist
            else:
                r_const_row[0][j] = min_dist
        if i == 0:
            Func_matrix_remake_2 = r_const_row
        else:
            Func_matrix_remake_2 = np.vstack((Func_matrix_remake_2, r_const_row))

    # print(Func_matrix == Func_matrix_remake)
    # print(min(Func_matrix.flat), max(Func_matrix.flat))
    # print(min(Func_matrix_remake.flat), max(Func_matrix_remake.flat))
    # print(new_bound_coord == bound_coord_remake)
    # print(Func_matrix_remake == Func_matrix_remake_2)
    # print(min(Func_matrix_remake_2.flat), max(Func_matrix_remake_2.flat))

    Func_matrix_list = [l for l in Func_matrix_remake_2.flat]
    Func_matrix_i = interpolate.griddata((X_list, Y_list), Func_matrix_list, (xig, yig), method='cubic')
    surf = plt.contour(xig, yig, Func_matrix_i, 0)
    zero_contour = surf.collections[1].get_paths()[0].vertices
    x_zero = zero_contour[:, 0]
    y_zero = zero_contour[:, 1]
    bound_coord_remake_2 = [(x_zero[i], y_zero[i]) for i in range(len(x_zero))]
    #print(bound_coord_remake == bound_coord_remake_2)
    print('done')



    return bound_coord_remake, Func_coord_dict, Func_matrix_remake


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
    surf = plt.contour(xig, yig, Func_matrix_i, 0)
    zero_contour = surf.collections[1].get_paths()[0].vertices
    x_zero = zero_contour[:, 0]
    y_zero = zero_contour[:, 1]
    plt.plot(x_zero, y_zero)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    set_bound_coord = set(bound_coord)
    fig1 = plt.figure()
    print(bound_coord)
    for coord_pair in set_bound_coord:
        plt.scatter(coord_pair[0], coord_pair[1])
    plt.show()



import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm

# frac_angle_1 = np.pi/6
# frac_angle_2 = np.pi/4*3
# frac_angle_list = [frac_angle_1, frac_angle_2]
# delta_r = 0.005
# delta_r_fine = 0.005
# R_for_fine = 0.015
# R = 0.215
# r_well = 0.0075
# N_r_fine = round(R_for_fine/delta_r_fine)
#
# delta_r_list = [delta_r_fine]*N_r_fine + [delta_r]*round((R-r_well-R_for_fine)/delta_r)
# N_r_full = len(delta_r_list)
#
# delta_fi = np.pi / 60 # угол-шаг в радианах
# delta_fi_fine = np.pi/180
# fi_for_fine = np.pi/18
# M_fi_fine = round(fi_for_fine / delta_fi_fine)
# delta_fi_list_first = [delta_fi]*round((frac_angle_1-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle_1 - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi))
# angle_lack = round((2*np.pi - sum(delta_fi_list_first))/delta_fi)
# delta_fi_list = [delta_fi]*round((frac_angle_1-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle_1 - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi)+angle_lack)
# M_fi_full = len(delta_fi_list)


def replace_boundary(frac_angle_1, frac_angle_2, r_well, delta_fi_list, delta_r_list, frac_length):
    frac_angle_list = [frac_angle_1, frac_angle_2]
    rect_length = frac_length + 0.01
    rect_width = 0.02
    Func_coord_dict = {}
    X_list = []
    Y_list = []
    Func_list = []
    for frac_angle in frac_angle_list:
        rect_center_x = (r_well + rect_length / 2) * np.cos(frac_angle)
        rect_center_y = (r_well + rect_length / 2) * np.sin(frac_angle)
        h1 = rect_width / 2 / np.cos(frac_angle)
        h2 = (r_well + rect_length) / np.sin(frac_angle)
        h3 = r_well / np.sin(frac_angle)
        for i in range(len(delta_r_list)):
            for j in range(len(delta_fi_list)):
                x = (r_well+sum(delta_r_list[0:i])) * np.cos(sum(delta_fi_list[0:j]))
                y = (r_well+sum(delta_r_list[0:i])) * np.sin(sum(delta_fi_list[0:j]))
                if frac_angle <= np.pi/2:
                    if y > h3 - x/np.tan(frac_angle) and y < h2 - x/np.tan(frac_angle) and y > -h1 + np.tan(frac_angle)*x and y < h1 + np.tan(frac_angle)*x:
                        Func_coord_dict[(i,j)] = (x,y,1)
                    else:
                        if (i,j) in Func_coord_dict and Func_coord_dict[(i,j)][2] == 1:
                            continue
                        else:
                            Func_coord_dict[(i,j)]=(x, y, -1)
                elif np.pi/2 < frac_angle <= np.pi:
                    if y > h3 - x/np.tan(frac_angle) and y < h2 - x/np.tan(frac_angle) and y > h1 + np.tan(frac_angle)*x and y < -h1 + np.tan(frac_angle)*x:
                        Func_coord_dict[(i, j)] = (x, y, 1)
                    else:
                        if (i,j) in Func_coord_dict and Func_coord_dict[(i,j)][2] == 1:
                            continue
                        else:
                            Func_coord_dict[(i,j)]=(x, y, -1)
                elif np.pi < frac_angle <= np.pi*3/2:
                    if y < h3 - x/np.tan(frac_angle) and y > h2 - x/np.tan(frac_angle) and y > h1 + np.tan(frac_angle)*x and y < -h1 + np.tan(frac_angle)*x:
                        Func_coord_dict[(i, j)] = (x, y, 1)
                    else:
                        if (i,j) in Func_coord_dict and Func_coord_dict[(i,j)][2] == 1:
                            continue
                        else:
                            Func_coord_dict[(i,j)]=(x, y, -1)
                elif np.pi*3/2 < frac_angle <= 2*np.pi:
                    if y < h3 - x/np.tan(frac_angle) and y > h2 - x/np.tan(frac_angle) and y > -h1 + np.tan(frac_angle)*x and y < h1 + np.tan(frac_angle)*x:
                        Func_coord_dict[(i, j)] = (x, y, 1)
                    else:
                        if (i,j) in Func_coord_dict and Func_coord_dict[(i,j)][2] == 1:
                            continue
                        else:
                            Func_coord_dict[(i,j)]=(x, y, -1)

    for key in Func_coord_dict:
        X_list.append(Func_coord_dict[key][0])
        Y_list.append(Func_coord_dict[key][1])
        Func_list.append(Func_coord_dict[key][2])
    #if __name__ == '__main__':
    xi = np.linspace(min(X_list),max(X_list), 500)
    yi = np.linspace(min(Y_list), max(Y_list), 500)
    xig, yig = np.meshgrid(xi, yi)
    Func_i = interpolate.griddata((X_list,Y_list), Func_list, (xig, yig), method='cubic')

    #levels = list(range(0,1,0.1))
    fig, ax = plt.subplots()
    surf = plt.contourf(xig, yig, Func_i, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Func_i), vmax=np.nanmax(Func_i))
    plt.show()
    #-------------------------------------------------------------------------------------
    keys_list_sorted = list(sorted(Func_coord_dict.keys()))
    bound_coord = []
    bound_coord_cell = []
    for i in range(len(keys_list_sorted)-1):
        if Func_coord_dict[keys_list_sorted[i]][2] != Func_coord_dict[keys_list_sorted[i+1]][2]:
            if Func_coord_dict[keys_list_sorted[i]][2] == 1:
                bound_coord.append((Func_coord_dict[keys_list_sorted[i]][0], Func_coord_dict[keys_list_sorted[i]][1]))
                bound_coord_cell.append(keys_list_sorted[i])
            elif Func_coord_dict[keys_list_sorted[i+1]][2] == 1:
                bound_coord.append((Func_coord_dict[keys_list_sorted[i]][0], Func_coord_dict[keys_list_sorted[i]][1]))
                bound_coord_cell.append(keys_list_sorted[i+1])

    print(len(bound_coord))
    def func(temp):
        return temp[1], temp[0]

    keys_list_sorted = list(sorted(Func_coord_dict.keys(), key=func))
    for i in range(len(keys_list_sorted)-1):
        if Func_coord_dict[keys_list_sorted[i]][2] != Func_coord_dict[keys_list_sorted[i+1]][2]:
            if Func_coord_dict[keys_list_sorted[i]][2] == 1:
                bound_coord.append((Func_coord_dict[keys_list_sorted[i]][0], Func_coord_dict[keys_list_sorted[i]][1]))
            elif Func_coord_dict[keys_list_sorted[i+1]][2] == 1:
                bound_coord.append((Func_coord_dict[keys_list_sorted[i]][0], Func_coord_dict[keys_list_sorted[i]][1]))

    print(len(bound_coord))
    print(len(set(bound_coord)))
    set_bound_coord = set(bound_coord)
    #if __name__ == '__main__':
    fig1 = plt.figure()
    print(bound_coord)
    for coord_pair in set_bound_coord:
        plt.scatter(coord_pair[0], coord_pair[1])
    plt.show()

    #-----------------------------------------------------------------------------------------

    coord_x_list = []
    coord_y_list = []
    all_dist_list = []
    for i in range(len(delta_r_list)):
        r_const_row = np.zeros((1, len(delta_fi_list)))
        for j in range(len(delta_fi_list)):
            coord_xy = (Func_coord_dict[(i,j)][0], Func_coord_dict[(i,j)][1])
            coord_x_list.append(coord_xy[0])
            coord_y_list.append(coord_xy[1])
            Func_value = Func_coord_dict[(i,j)][2]
            dist_list = []
            for bound_coord_pair in bound_coord:
                dist = ((bound_coord_pair[0]-coord_xy[0])**2+(bound_coord_pair[1]-coord_xy[1])**2)**0.5
                dist_list.append(dist)
            min_dist = min(dist_list)
            all_dist_list.append(min_dist)
            if Func_value < 0:
                r_const_row[0][j] = -min_dist
            else:
                r_const_row[0][j] = min_dist
        if i == 0:
            Func_matrix = r_const_row
        else:
            Func_matrix = np.vstack((Func_matrix, r_const_row))

    #if __name__ == '__main__':
    xi = np.linspace(min(coord_x_list), max(coord_x_list), 700)
    yi = np.linspace(min(coord_y_list), max(coord_y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((coord_x_list, coord_y_list), all_dist_list, (xig, yig), method='cubic')

    #levels = list(range(0, 5000000, 10000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    plt.show()
    return bound_coord, Func_matrix, bound_coord_cell


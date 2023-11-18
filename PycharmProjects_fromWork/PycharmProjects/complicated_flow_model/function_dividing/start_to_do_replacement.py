import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm


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

    def func(temp):
        return temp[1], temp[0]

    keys_list_sorted = list(sorted(Func_coord_dict.keys(), key=func))
    for i in range(len(keys_list_sorted)-1):
        if Func_coord_dict[keys_list_sorted[i]][2] != Func_coord_dict[keys_list_sorted[i+1]][2]:
            if Func_coord_dict[keys_list_sorted[i]][2] == 1:
                bound_coord.append((Func_coord_dict[keys_list_sorted[i]][0], Func_coord_dict[keys_list_sorted[i]][1]))
            elif Func_coord_dict[keys_list_sorted[i+1]][2] == 1:
                bound_coord.append((Func_coord_dict[keys_list_sorted[i]][0], Func_coord_dict[keys_list_sorted[i]][1]))

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

    return bound_coord, Func_matrix


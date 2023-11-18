import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
import shapely.geometry as geom
from scipy.spatial import ConvexHull

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

import numpy as np
from matplotlib import pyplot as plt
from math import pi, cos, sin

def find_oval(t_rot, b, a):
    u = 0.008*np.cos(t_rot)       #x-position of the center
    v = 0.008*np.sin(t_rot)
    # u = 0
    # v = 0 #y-position of the center
    # a=0.1       #radius on the x-axis
    # b=0.005     #radius on the y-axis
    # t_rot=pi/4 #rotation angle

    t = np.linspace(0, pi/2, 50) + np.linspace(3*pi/2, 2*pi, 50)
    Ell = np.array([a*np.cos(t) , b*np.sin(t)])
         #u,v removed to keep the same center location
    R_rot = np.array([[cos(t_rot) , -sin(t_rot)],[sin(t_rot) , cos(t_rot)]])
         #2-D rotation matrix

    Ell_rot = np.zeros((2,Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])

    plt.plot( u+Ell[0,:] , v+Ell[1,:] )     #initial ellipse

    plt.plot( u+Ell_rot[0,:] , v+Ell_rot[1,:],'darkorange' )    #rotated ellipse
    plt.grid(color='lightgray',linestyle='--')
    plt.show()

    return Ell_rot

def replace_boundary(frac_angle_1, frac_angle_2, r_well, frac_length, M_fi_full, N_r_full, coord_matrix):
    Func_matrix = np.zeros((N_r_full, M_fi_full))
    u_1 = 0.008 * np.cos(frac_angle_1)  # x-position of the center
    v_1 = 0.008 * np.sin(frac_angle_1)
    # u_1 = 0
    # v_1 = 0
    oval_length = frac_length/2 + 0.001
    oval_width = r_well
    half_oval_1 = find_oval(frac_angle_1, oval_width, oval_length)
    x_y_bound_1 = list(zip(half_oval_1[0]+u_1, half_oval_1[1]+v_1))

    u_2 = 0.008 * np.cos(frac_angle_2)  # x-position of the center
    v_2 = 0.008 * np.sin(frac_angle_2)
    # u_2 = 0
    # v_2 = 0
    half_oval_2 = find_oval(frac_angle_2, oval_width, oval_length)
    x_y_bound_2 = list(zip(half_oval_2[0]+u_2, half_oval_2[1]+v_2))

    plt.plot( half_oval_1[0,:]+u_1, half_oval_1[1,:]+v_1)    #rotated ellipse
    plt.plot( half_oval_2[0,:] + u_2, half_oval_2[1,:] + v_2)    #rotated ellipse

    plt.grid(color='lightgray',linestyle='--')
    plt.show()



    hull_1 = ConvexHull(x_y_bound_1).vertices
    sort_bound_coords_1 = [x_y_bound_1[i] for i in hull_1]
    print(sort_bound_coords_1)

    sort_bound_coords_1.append(sort_bound_coords_1[0])
    print(len(x_y_bound_1), len(sort_bound_coords_1))

    hull_2 = ConvexHull(x_y_bound_2).vertices
    sort_bound_coords_2 = [x_y_bound_2[i] for i in hull_2]
    print(sort_bound_coords_2)

    sort_bound_coords_2.append(sort_bound_coords_2[0])
    print(len(x_y_bound_2), len(sort_bound_coords_2))

    # print(bound_coords)
    for elem in sort_bound_coords_1:
        plt.scatter(elem[0], elem[1])
    for elem in sort_bound_coords_2:
        plt.scatter(elem[0], elem[1])
    plt.Circle((0, 0), 0.008, color='r')
    plt.show()

    line_1 = geom.LineString(sort_bound_coords_1)
    polygon_1 = geom.Polygon(sort_bound_coords_1)

    line_2 = geom.LineString(sort_bound_coords_2)
    polygon_2 = geom.Polygon(sort_bound_coords_2)

    for m in range(M_fi_full):
        for n in range(N_r_full):
            point = geom.Point(coord_matrix[n][m])
            if polygon_1.contains(point):
                Func_matrix[n][m] = point.distance(line_1)
            elif polygon_2.contains(point):
                Func_matrix[n][m] = point.distance(line_2)
            else:
                point_dist_line_1 = point.distance(line_1)
                point_dist_line_2 = point.distance(line_2)
                Func_matrix[n][m] = -min([point_dist_line_1, point_dist_line_2])



    # line_2 = geom.LineString(sort_bound_coords_2)
    # polygon_2 = geom.Polygon(sort_bound_coords_2)
    #
    # for m in range(M_fi_full):
    #     for n in range(N_r_full):
    #         point = geom.Point(coord_matrix[n][m])
    #         if polygon_2.contains(point):
    #             Func_matrix[n][m] = point.distance(line_2)
    #         else:
    #             Func_matrix[n][m] = -point.distance(line_2)

    return Func_matrix

delta_r = 0.0001
delta_r_fine = 0.0001
R_for_fine = 0.02
R = 0.215
r_well = 0.0075
N_r_fine = round(R_for_fine / delta_r_fine)
delta_r_list = [delta_r_fine] * N_r_fine + [delta_r] * round((R - r_well - R_for_fine) / delta_r)
N_r_full = len(delta_r_list)

delta_fi = np.pi / 180  # угол-шаг в радианах
delta_fi_fine = np.pi / 180
fi_for_fine = np.pi / 6
M_fi_fine = round(fi_for_fine / delta_fi_fine)

frac_angle = np.pi / 6
frac_angle_2 = np.pi + np.pi/6

delta_fi_list_first = [delta_fi] * round((frac_angle - fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine * 2) + [
    delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine * 2) + [
                          delta_fi] * (round((2 * np.pi - frac_angle_2 - fi_for_fine) / delta_fi))
angle_lack = round((2 * np.pi - sum(delta_fi_list_first)) / delta_fi)
# delta_fi_list = [delta_fi]*round((frac_angle-fi_for_fine)/delta_fi) + [delta_fi_fine]*(M_fi_fine*2) + [delta_fi] * round((frac_angle_2 - frac_angle - 2 * fi_for_fine) / delta_fi) + [delta_fi_fine] * (M_fi_fine*2) + [delta_fi] * (round((2*np.pi - frac_angle_2 - fi_for_fine)/delta_fi)+angle_lack)
delta_fi_list = [delta_fi] * round(2 * np.pi / delta_fi)
M_fi_full = len(delta_fi_list)

coord_matrix_rad = []
coord_matrix_cart = []
Y_list = []
X_list = []
for n in range(len(delta_r_list)):
    coord_line_rad = []
    coord_line_cart = []
    r = sum(delta_r_list[0:n]) + r_well
    for m in range(len(delta_fi_list)):
        fi = sum(delta_fi_list[0:m])
        coord_line_rad.append((r, fi))
        coord_line_cart.append((r * np.cos(fi), r * np.sin(fi)))
        X_list.append(r * np.cos(fi))
        Y_list.append(r * np.sin(fi))
        # coord_matrix_rad[n][m] = (r, fi)
        # coord_matrix_cart[n][m] = (r*np.cos(fi), r*np.sin(fi))
    coord_matrix_rad.append(coord_line_rad)
    coord_matrix_cart.append(coord_line_cart)
    #print(min(coord_matrix_cart), max(coord_matrix_cart))


fig = plt.figure()
func_matrix = replace_boundary(frac_angle, frac_angle_2, 0.003, 0.2, M_fi_full, N_r_full, coord_matrix_cart)

xi = np.linspace(min(X_list),max(X_list), 500)
yi = np.linspace(min(Y_list), max(Y_list), 500)
xig, yig = np.meshgrid(xi, yi)
Func_i = interpolate.griddata((X_list,Y_list), func_matrix.flat, (xig, yig), method='cubic')
surf = plt.contourf(xig, yig, Func_i, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Func_i), vmax=np.nanmax(Func_i))
plt.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

# surf = plt.contourf(Y_list, X_list, func_matrix, cmap=cm.jet, antialiased=True, vmin=np.nanmin(func_matrix.flat), vmax=np.nanmax(func_matrix.flat),linewidth=0.2)
#
# plt.show()


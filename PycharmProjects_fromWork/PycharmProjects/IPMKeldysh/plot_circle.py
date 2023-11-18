from numpy import *
from matplotlib import pyplot as plt
r = 0.215
t = arange(0, 2*pi, 0.01)
x = r*sin(t)
y = r*cos(t)
fig = plt.plot(x, y)
x_sensor = [-0.121, -0.057, 0.057, 0.121, 0, 0.065, 0.127, -0.127, -0.121, 0.057, 0, 0.07, 0, 0, 0, -0.185, -0.057]
y_sensor = [0.121, 0.127, 0.127, 0.121, 0.07, 0.065, 0, 0, -0.121, -0.127, -0.185, 0, 0.127, -0.07, 0, 0, -0.127]
axis_font = {'size':'20'}
plt.xlabel('X, m', fontsize=16 )
plt.ylabel('Y, m', fontsize=16)
plt.scatter(x_sensor, y_sensor, c='g', marker='o')
# plt.scatter((-0.121, 0.121), (-0.057, 0.127), (0.057, 0.127),
#             (0.121, 0.121))
# plt.scatter((0, 0.07), (0.065, 0.065), (0.127, 0), (-0.127, 0))
# plt.scatter((-0.121, -0.121), (0.057, -0.127), (0, -0.185), (0.07, 0))
#
# plt.scatter((0, 0.127), (0, -0.07), (0, 0), (-0.185, 0))
# plt.scatter((-0.057, -0.127))
plt.tick_params(axis='both', labelsize=16)
plt.show()

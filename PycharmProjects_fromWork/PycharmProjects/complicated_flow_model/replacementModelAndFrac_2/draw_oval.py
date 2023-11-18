import numpy as np
from matplotlib import pyplot as plt
from math import pi, cos, sin

u=0.       #x-position of the center
v=0.      #y-position of the center
a=0.1       #radius on the x-axis
b=0.005     #radius on the y-axis
t_rot=pi/4 #rotation angle

t = np.linspace(0, pi/2, 50) + np.linspace(3*pi/2, 2*pi, 50)
Ell = np.array([a*np.cos(t) , b*np.sin(t)])
     #u,v removed to keep the same center location
R_rot = np.array([[cos(t_rot) , -sin(t_rot)],[sin(t_rot) , cos(t_rot)]])
     #2-D rotation matrix

Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])

plt.plot( u+Ell[0,:] , v+Ell[1,:] )     #initial ellipse
plt.plot( u+Ell_rot[0,:] , v+Ell_rot[1,:],'darkorange' ) #rotated ellipse
print(u+Ell_rot[0,:])
print(v+Ell_rot[1,:])
print(Ell_rot)
plt.grid(color='lightgray',linestyle='--')
plt.show()

plt.plot()
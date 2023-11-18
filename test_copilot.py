
# draw plot with parabola using matplotlib
import matplotlib.pyplot as plt
import numpy as np

# draw ellipsoid with given length and width
def draw_ellipsoid(length, width):
    x = np.linspace(0, length, 100)
    y = width * np.sqrt(1 - (x/length)**2)
    plt.plot(x, y)
    plt.show()


def draw_graphic():
    x = np.linspace(0, 10, 100)
    y = x**2
    plt.plot(x, y)
    plt.show()

#draw_graphic()

draw_ellipsoid(10, 5)
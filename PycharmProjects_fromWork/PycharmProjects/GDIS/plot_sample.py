import matplotlib.pyplot as plt

circle1 = plt.Circle((0, 0), 0.215, color='r', fill=False)
circle2 = plt.Circle((0.121, 0.121), 0.0075, color='blue')
circle3 = plt.Circle((-0.121, -0.121), 0.0075, color='blue')

circle4 = plt.Circle((-0.127, 0), 0.003, color='g', fill=True)
circle5 = plt.Circle((0.127, 0), 0.003, color='g', fill=True)
circle6 = plt.Circle((0.057, -0.127), 0.003, color='g', fill=True)
circle7 = plt.Circle((-0.057, 0.127), 0.003, color='g', fill=True)
circle8 = plt.Circle((0, 0.07), 0.003, color='g', fill=True)
circle9 = plt.Circle((0, -0.185), 0.003, color='g', fill=True)
circle10 = plt.Circle((0.065, 0.065), 0.003, color='g', fill=True)
circle11 = plt.Circle((-0.121, 0.121), 0.003, color='g', fill=True)
circle12= plt.Circle((0.057, 0.127), 0.003, color='g', fill=True)

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
ax.add_artist(circle4)
ax.add_artist(circle5)
ax.add_artist(circle6)
ax.add_artist(circle7)
ax.add_artist(circle8)
ax.add_artist(circle9)
ax.add_artist(circle10)
ax.add_artist(circle11)
ax.add_artist(circle12)

plt.grid()

ax.set_xlim((-0.3, 0.3))
ax.set_ylim((-0.3, 0.3))

plt.show()

from Visualize.mpl_visualize import plot_balls, plot_circles, plot_edges
import matplotlib.pyplot as plt
from System.sys_funcs.calcs.circle import calc_circ
from System.Network.edges.build_edge import build_edge
import numpy as np


locs = np.array([[-5, 0, 0], [5, 0, 0], [0, 0.5, 0]])
rads = [3, 3, 1.0]

my_circ = calc_circ(*locs, *rads, return_both=True)

# print(np.linalg.norm([]))


my_edge = build_edge(locs, rads, [np.array([0.0, -0.12018873,  5.12002788]), np.array([0.0, -0.12007157, -5.12007157])], res=0.01)
my_edge1 = build_edge(locs, rads, [np.array([0.0, -0.12018873,  5.12002788]), my_circ[2]], res=0.01)
my_edge2 = build_edge(locs, rads, [my_circ[2], np.array([0.0, -0.12007157, -5.12007157])], res=0.01)


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

plot_circles([my_circ[0], my_circ[2]], [my_circ[1], my_circ[3]], fig=fig, ax=ax, colors=['k', 'k'])

plot_balls(locs, rads, fig=fig, ax=ax, colors=['b', 'b', 'b'], res=10, alpha=0.5)

plot_edges([my_edge[0], my_edge1[0], my_edge2[0]], fig=fig, ax=ax, thickness=1, colors=['r', 'r', 'r'])
ax.set_xlim([-10, 10])
ax.set_ylim([-10, 10])
ax.set_zlim([-10, 10])
plt.show()

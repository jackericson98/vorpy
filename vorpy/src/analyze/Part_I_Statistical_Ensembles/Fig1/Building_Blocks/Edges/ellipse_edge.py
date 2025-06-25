import matplotlib.pyplot as plt

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.visualize.mpl_visualize import plot_balls, plot_circles, plot_edges, plot_surfs
from vorpy.src.calculations.edge import calc_circ
from vorpy.src.network.build_edge import build_edge
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

plot_balls([my_circ[0], my_circ[2]], [my_circ[1], my_circ[3]], fig=fig, ax=ax, colors=['k', 'k'], alpha=0.1)
plot_circles([my_circ[0], my_circ[2]], [my_circ[1], my_circ[3]], fig=fig, ax=ax, colors=['k', 'k'])
ax.scatter([my_circ[0][0], my_circ[2][0], locs[0][0], locs[1][0], locs[2][0]],
           [my_circ[0][1], my_circ[2][1], locs[0][1], locs[1][1], locs[2][1]],
           [my_circ[0][2], my_circ[2][2], locs[0][2], locs[1][2], locs[2][2]], c='k', alpha=0.5)

plot_balls(locs, rads, fig=fig, ax=ax, colors=['b', 'b', 'b'], res=10, alpha=0.5)

plot_edges([my_edge[0], my_edge1[0], my_edge2[0]], fig=fig, ax=ax, thickness=1, colors=['r', 'r', 'r'])


ax.set_xlim([-10, 10])
ax.set_ylim([-10, 10])
ax.set_zlim([-10, 10])
plt.show()

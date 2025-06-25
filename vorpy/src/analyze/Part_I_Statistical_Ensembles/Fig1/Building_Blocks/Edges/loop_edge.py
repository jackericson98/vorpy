import numpy as np
import matplotlib.pyplot as plt

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.visualize.mpl_visualize import plot_balls, plot_edges, plot_verts
from vorpy.src.network.build_edge import build_edge

atoms = [([0, 0, 10.25], 10), ([0, 0, -10.25], 10), ([0, 0, 0], 0.25)]
# verts = [[0.5, 0.5, 0],
#          [-2 * (1 + np.sqrt(2)) / 8, 2 * (np.sqrt(2) - 1) / 8, 0],
#          [32 * (np.sqrt(2) - 1) / 8, - 2 * (1 + np.sqrt(2)) / 8, 0]]
verts1 = [[-0.50637631, -0.08104263,  0.        ], [ 0.38077303, -0.34350659,  0.        ], [-0.10890493,  0.50112333,  0.        ]]
edges = [build_edge([_[0] for _ in atoms], [_[1] for _ in atoms], [verts1[0], verts1[1]], res=0.02, straight=False),
         build_edge([_[0] for _ in atoms], [_[1] for _ in atoms], [verts1[1], verts1[2]], res=0.02)]
         # build_edge([_[0] for _ in atoms], [_[1] for _ in atoms], [verts[2], verts[0]], res=0.02)]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')


plot_balls([_[0] for _ in atoms], [_[1] for _ in atoms], fig=fig, ax=ax, alpha=0.2, res=8)
# plot_verts(verts1, [0.1, 0.1, 0.1], spheres=True, fig=fig, ax=ax)
plot_edges([_[0] for _ in edges], fig=fig, ax=ax, thickness=2)

# Set the scales for the figure
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)

plt.show()

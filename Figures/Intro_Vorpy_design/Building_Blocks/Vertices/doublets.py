import numpy as np
from numpy import array as ar
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_atoms, plot_verts
from System.Network.verts.calc_vert import calc_vert


"""
Doublet plotting code. Choose doublet type from below

    1. Doublet type 1 - 2 Edges, 1 Surface, All Balls Equal
    2. Doublet type 1 - 2 Edges, 1 Surface, All Balls Different
    3. Doublet type 2 - 3 Edges, 3 Surfaces, All Balls Equal
    4. Doublet type 2 - 3 Edges, 3 Surfaces, All Balls Different

"""

# Choose Doublet type here V
doublet_type = 4


# Choose other settings here
show_edges = False
show_surfs = False
atom_alpha = 0.5


# Set my vert to None, for fake generation
my_vert = None

# Type 1 Doublet - Perfect
if doublet_type == 1:
    rads = [1.1, 1.1, 1.0, 1.0]
    locs = [2.1, 0.0, 0.0], [-2.1, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, -2.0, 0.0]
    my_vert = ([0.0, 0.0, 2.5], 2.2, [0.0, 0.0, -2.5], 2.2)
    title = 'Doublet Type 1'

# No Overlap Vertex 1
elif doublet_type == 2:
    rads = [1.6, 1.5, 1.5, 0.72]
    locs = [1.31, -0.02, -0.9], [-1.03, 1.74, 0.58], [2.17, 0.21, 1.72], [0.35, 0.95, 1.63]
    title = 'Doublet Type 1, Verts on Either Side'

# Type 1 Doublet - Both on the same size
elif doublet_type == 3:
    rads = [1.8, 1.8, 1.3, 1.3]
    locs = [0.23, 0.39, 0.5], [-1.1, -0.37, 0.77], [0.3, 0.67, -0.56], [-1.82, 0.41, -1.48]
    title = 'Doublet Type 1, Verts on Same side'

# Type 2 Doublet -
elif doublet_type == 4:
    rads = [1.0, 1.0, 1.0, 0.5]
    locs = [1.5, np.sqrt(3)/2, 0.0], [-1.5, np.sqrt(3)/2, 0.0], [0.0, -np.sqrt(3), 0.0], [0.0, 0.0, 0.0]
    my_vert = ([0.0, 0.0, 2.75], 2.25, [0.0, 0.0, -2.75], 2.25)
    title = 'Doublet Type 2'


# Calculate the vertex
if my_vert is None:
    my_vert = calc_vert(locs=ar([ar(_) for _ in locs]), rads=ar(rads))


# Make the plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')


# Plot the atoms
plot_atoms(locs, rads, fig=fig, ax=ax, res=10, alpha=atom_alpha)
# Plot the vertices
plot_verts([my_vert[0]], [abs(my_vert[1])], fig=fig, ax=ax, spheres=True, res=10)
plot_verts([my_vert[2]], [abs(my_vert[3])], fig=fig, ax=ax, spheres=True, res=10)
# Plot the edges


# Set the axes lines
# ax.plot([-3, -2], [-3, -3], [0, 0])
# ax.plot([-3, -3], [-2, -3], [0, 0])
# ax.plot([-3, -3], [-3, -3], [0, 1])
#
# # Set the axes labels
# ax.text(x=-2, y=-3, z=-0.25, s='x')
# ax.text(x=-3, y=-2, z=-0.25, s='y')
# ax.text(x=-3, y=-3, z=1, s='z')

# Set the scales for the figure
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)

# Set the title
ax.set_title(title, font=dict(size=20, family='serif'))

# Show the plot
plt.show()

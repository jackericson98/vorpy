import numpy as np
from numpy import array as ar
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_balls, plot_verts, plot_edges, plot_surfs, plot_circles
from System.Network.verts.calc_vert import calc_vert
from System.Network.edges.build_edge import build_edge
from System.sys_funcs.calcs.circle import calc_circ


"""
Edge plotting code. Choose an edge type below.

    1. Flat Edge: All balls equal
    2. Semi Curved Edge 1: No overlap, 2 balls the same size, relatively close
    3. Semi Curved Edge 2: No overlap, all balls different sizes, relatively close
    4. Semi Curved Edge 3: Some overlap, 2 balls the same size, very close
    5. Semi curved Edge 4. Some overlap, all balls different sizes, very close
    ...
"""

# Choose here
edge_choice = 8

# Additional Settings
atom_alpha = 0.4
vert_alpha = 0.1
edge_thickness = 2

show_edge = True
show_verts = True
show_vert_spheres = False
show_atoms = True
show_circs = True

edge_color = 'red'
atom_color = 'blue'

# Declare variables
title = ''


# Flat edge
if edge_choice == 1:
    locs = [1.5, np.sqrt(3)/2, 0.0], [-1.5, np.sqrt(3)/2, 0.0], [0.0, -np.sqrt(3), 0.0]
    rads = 1.0, 1.0, 1.0
    title = 'Flat Edge - All Balls Equal'

# Semi Curved edge 1
elif edge_choice == 2:
    locs = [1.5, np.sqrt(3)/2, 0.0], [-1.5, np.sqrt(3)/2, 0.0], [0.0, -np.sqrt(3), 0.0]
    rads = 1.0, 0.75, 0.75
    title = 'Curved Edge - No Overlap'

# Semi Curved Edge 2
elif edge_choice == 3:
    locs = [1.5, np.sqrt(3)/2, 0.0], [-1.5, np.sqrt(3)/2, 0.0], [0.0, -np.sqrt(3), 0.0]
    rads = 1.0, 0.75, 0.5

# Semi Curved Edge 3
elif edge_choice == 4:
    locs = [1.5, np.sqrt(3)/2, 0.0], [-1.5, np.sqrt(3)/2, 0.0], [0.0, -np.sqrt(3), 0.0]
    rads = 1.0, 2.0, 2.0
    title = 'Curved Edge - Two Overlapping Atoms'

elif edge_choice == 5:
    locs = [-5.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 1.0, 0.0]
    rads = 4.0, 4.0, 2.0
    title = ''

# Semi Curved Edge 4
elif edge_choice == 6:
    locs = [1.5, np.sqrt(3)/2, 0.0], [-1.5, np.sqrt(3)/2, 0.0], [0.0, -np.sqrt(3), 0.0]
    rads = 1.0, 1.0, 2.5
    title = 'Curved Edge - All Atoms Overlap'


# Semi Curved Edge 6
elif edge_choice == 7:
    locs = [1.5, np.sqrt(3)/2, 0.0], [-1.5, np.sqrt(3)/2, 0.0], [0.0, -np.sqrt(3), 0.0]
    rads = 1.65, 1.65, 3.2
    title = 'Curved Edge - All Atoms Overlap'

elif edge_choice == 8:
    locs = [0.0, 0.0, 0.0], [5.0, 0.1, 0.0], [-5.0, 0.2, 0.0]
    rads = 3.0, 1.1, 0.9

# Set the radii and distances for the surrounding atoms
r, d = 0.5, 15.0
my_vert_atoms = [(ar([0.0, 0.0, d]), r), (ar([0.0, 0.0, -d]), r)]

# Create the edge atoms
edge_atoms = [(ar(locs[i]), rads[i]) for i in range(3)]


# Calculate the vertices
vert_atoms = [edge_atoms + [_] for _ in my_vert_atoms]
my_verts = [calc_vert(ar([_[0] for _ in my_atoms]), ar([_[1] for _ in my_atoms])) for my_atoms in vert_atoms]

# Calculate the Edge
if edge_choice < 8:
    my_edge = build_edge(locs=ar([_[0] for _ in edge_atoms]), rads=ar([_[1] for _ in edge_atoms]),
                         vlocs=ar([_[0] for _ in my_verts]), res=0.5)


# Make the plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Plot the edge
if show_edge and edge_choice < 8:
    plot_edges([my_edge[0]], fig=fig, ax=ax, colors=[edge_color], thickness=edge_thickness)
# Plot the vertices
if show_verts and edge_choice < 8:
    plot_verts([_[0] for _ in my_verts], [_[1] for _ in my_verts], fig=fig, ax=ax, spheres=show_vert_spheres, alpha=vert_alpha)
# Plot the atoms
if show_atoms:
    plot_balls(alocs=[_[0] for _ in edge_atoms], arads=[_[1] for _ in edge_atoms], colors=[atom_color for _ in edge_atoms], alpha=atom_alpha, fig=fig, ax=ax)

# Plot the circs
if show_circs and edge_choice < 8:
    my_circ = calc_circ(*locs, *rads, return_both=True)
    plot_circles([my_circ[0]], [my_circ[1]], fig=fig, ax=ax, center_point=True)
    if len(my_circ) == 4:
        plot_circles([my_circ[2]], [my_circ[3]], fig=fig, ax=ax, center_point=True)
# Set the scales for the figure
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)

# Set the title for the plot
ax.set_title(title, font=dict(size=20, family='serif'))
plt.show()

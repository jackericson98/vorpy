from numpy import array as ar
import numpy as np
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_atoms, plot_verts, plot_edges, plot_surfs
from System.Network.verts.calc_vert import calc_vert, calc_flat_vert
from System.Network.edges.build_edge import build_edge
from System.Network.surfs.build_surf import build_surf


"""
Plotting Surfaces. Choose a surface type from below
    
    1. Flat Surface - Equally Sized Balls
    2. Curved Surface - Ball Size Ratio 1:2, Far Apart
    3. Curved Surface - Ball Size Ratio 1:2, Close but not Touching
    4. Curved Surface - Ball Size Ratio 9:10, Balls Kissing
    5. Curved Surface - Ball Size Ratio 1:2, Balls Kissing
    6. Curved Surface - Ball Size Ratio 9:10, Balls Slightly Overlapping
    7. Curved Surface - Ball Size Ratio 1:2, Balls Slightly Overlapping
    8. Curved Surface - Ball Size Ratio 9:10, Balls Overlapping
    9. Curved Surface - Ball Size Ratio 1:2, Balls Overlapping

"""

"""Set the Surface number from above here"""
surface_type = 7


""" Other Settings """

# Atoms
show_atoms = True
atom_alpha = 0.2

# Vertices
show_verts = True
show_vert_balls = False
vert_alpha = 0.1

# Edges
show_edges = True
edge_thickness = 2

# Surfaces
show_surf = True
surf_res = 0.5
surf_alpha = 0.2

# Plot
show_title = True
show_axes = False
xyz_min, xyz_max = -5, 5


# Different surface types
net_type = 'vor'
xlocs, rads, title = None, None, None


# 1. Flat Surface
if surface_type == 1:
    rads = 1.0, 1.0
    xlocs = -5.0, 5.0
    net_type = 'pow'
    title = 'Flat Surface - Equally Sized Balls'

# 2. Curved Surface 1
elif surface_type == 2:
    rads = 0.5, 1.0
    xlocs = -3.0, 3.0
    title = 'Curved Surface - Ball Size Ratio 1:2, Far Apart'

# 3. Curved Surface 2
elif surface_type == 3:
    rads = 0.5, 1.0
    xlocs = -1.0, 1.0
    title = 'Curved Surface - Ball Size Ratio 1:2, Close'

# 4. Curved Surface 3
elif surface_type == 4:
    rads = 0.9, 1.0
    xlocs = -0.95, 0.95
    title = 'Curved Surface - Ball Size Ratio 9:10, Balls Kissing'

# 5. Curved Surface 4
elif surface_type == 5:
    rads = 0.5, 1.0
    xlocs = -0.75, 0.75
    title = 'Curved Surface - Ball Size Ratio 1:2, Balls Kissing'

# 6. Curved Surface 5
elif surface_type == 6:
    rads = 0.9, 1.0
    xlocs = -0.92, 0.92
    title = 'Curved Surface - Ball Size Ratio 9:10, Balls Slightly Overlapping'

# 7. Curved Surface 6
elif surface_type == 7:
    rads = 0.5, 1.0
    xlocs = -0.72, 0.72
    title = 'Curved Surface -  Ball Size Ratio 1:2, Slight Overlap'

# 8. Curved Surface 7
elif surface_type == 8:
    rads = 0.9, 1.0
    xlocs = -0.75, 0.75
    title = 'Curved Surface - Ball Size Ratio 9:10, Balls Overlapping'

# 9. Curved Surface 8
elif surface_type == 9:
    rads = 0.5, 1.0
    xlocs = -0.55, 0.55
    title = 'Curved Surface - Ball Size Ratio 1:2, Balls Overlapping'

# Set the radii and distances for the surrounding atoms
r, d = 0.5, 20.0
my_vert_atoms = [(ar([0.0, d, d]), r), (ar([0.0, -d, d]), r), (ar([0.0, -d, -d]), r), (ar([0.0, d, -d]), r)]
# Create the Surface Atoms
my_surf_atoms = [(ar([xlocs[0], 0.0, 0.0]), rads[0]), (ar([xlocs[1], 0.0, 0.0]), rads[1])]

# Calculate the vertices
vert_atoms = [[_ for _ in my_surf_atoms] + [my_vert_atoms[i], my_vert_atoms[(i+1) % 4]] for i in range(4)]
if net_type == 'vor':
    my_verts = [calc_vert(ar([_[0] for _ in my_atoms]), ar([_[1] for _ in my_atoms])) for my_atoms in vert_atoms]
else:
    my_verts = [calc_flat_vert(ar([_[0] for _ in my_atoms]), ar([_[1] for _ in my_atoms]),
                               power=True if net_type == 'pow' else False) for my_atoms in vert_atoms]

# Calculate the Edges
edge_atoms = [[_ for _ in my_surf_atoms] + [my_vert_atoms[(i+1)%4]] for i in range(4)]
my_edge_verts = [(my_verts[j], my_verts[(j+1)%4]) for j in range(4)]
my_edges = [build_edge(alocs=ar([_[0] for _ in edge_atoms[i]]), arads=ar([_[1] for _ in edge_atoms[i]]),
                       vlocs=ar([_[0] for _ in my_edge_verts[i]]), res=0.5,
                       straight=False if net_type == 'vor' else True) for i in range(4)]

# Calculate the surfaces
my_surf = build_surf(alocs=[_[0] for _ in my_surf_atoms], arads=[_[1] for _ in my_surf_atoms],
                     epnts=[_[0] for _ in my_edges], res=surf_res, net_type=net_type)

# Create the figure
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Plot the surfaces
if show_surf:
    plot_surfs([my_surf[0]], [my_surf[1]], fig=fig, ax=ax, alpha=surf_alpha)

# Plot the edges
if show_edges:
    plot_edges([_[0] for _ in my_edges], fig=fig, ax=ax, colors=['k'] * 4, thickness=edge_thickness)

# Plot the vertices
if show_verts:
    plot_verts([_[0] for _ in my_verts], [_[1] for _ in my_verts], fig=fig, ax=ax, colors=['r']*4,
               spheres=show_vert_balls, alpha=vert_alpha)
# Plot the atoms
if show_atoms:
    plot_atoms(alocs=[_[0] for _ in my_surf_atoms], arads=[_[1] for _ in my_surf_atoms], alpha=atom_alpha, fig=fig, ax=ax,
               res=10)

# Set the scales for the figure
ax.set_xlim(xyz_min, xyz_max)
ax.set_ylim(xyz_min, xyz_max)
ax.set_zlim(xyz_min, xyz_max)

# Set the axes lines
if show_axes:
    ax.plot([-3, -2], [-3, -3], [0, 0])
    ax.plot([-3, -3], [-2, -3], [0, 0])
    ax.plot([-3, -3], [-3, -3], [0, 1])

    # Set the axes labels
    ax.text(x=-2, y=-3, z=-0.25, s='x')
    ax.text(x=-3, y=-2, z=-0.25, s='y')
    ax.text(x=-3, y=-3, z=1, s='z')

# Set the title for the figure
if show_title:
    ax.set_title(title, font=dict(size=20, family='serif'))

# Show the plot
plt.show()

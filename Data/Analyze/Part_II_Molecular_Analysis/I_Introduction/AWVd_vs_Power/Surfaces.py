from numpy import array as ar
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_balls, plot_verts, plot_edges, plot_surfs
from System.sys_funcs.calcs.vert import calc_vert
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
    7. Curved Surface -  Ball Size Ratio 1:2, Balls Slightly Overlapping
    8. Curved Surface - Ball Size Ratio 9:10, Balls Overlapping
    9. Curved Surface - Ball Size Ratio 1:2, Balls Overlapping
    10. 
"""

# Set the surface type here:
surface_type = 9

# Other settings
res = 0.5
show_edges = True
show_verts = False

# Different surface types
net_type = 'vor'

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
    title = 'Curved Surface - Ball Size Ratio 1:2, Close but not Touching'

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
    title = 'Curved Surface -  Ball Size Ratio 1:2, Balls Slightly Overlapping'

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
vert_atoms = [[_ for _ in my_surf_atoms] + [my_vert_atoms[i], my_vert_atoms[(i + 1) % 4]] for i in range(4)]
my_verts = [calc_vert(ar([_[0] for _ in my_atoms]), ar([_[1] for _ in my_atoms])) for my_atoms in vert_atoms]

# ## Plot the verts and their atoms
# for i in range(4):
#     fig = plt.figure()
#     ax = fig.add_subplot(projection='3d')
#     plot_atoms([_[0] for _ in vert_atoms[i]], [_[1] for _ in vert_atoms[i]], fig=fig, ax=ax)
#     plot_verts([my_verts[i][0]], [my_verts[i][1]], fig=fig, ax=ax, Show=True, spheres=True)

# Calculate the Edges
edge_atoms = [[_ for _ in my_surf_atoms] + [my_vert_atoms[(i + 1) % 4]] for i in range(4)]
my_edge_verts = [(my_verts[j], my_verts[(j + 1) % 4]) for j in range(4)]
my_edges = [build_edge(locs=ar([_[0] for _ in edge_atoms[i]]), rads=ar([_[1] for _ in edge_atoms[i]]),
                       vlocs=ar([_[0] for _ in my_edge_verts[i]]), res=0.5) for i in range(4)]

# Plot the edges and their atoms
# for i in range(4):
#     fig = plt.figure()
#     ax = fig.add_subplot(projection='3d')
#     plot_atoms([_[0] for _ in edge_atoms[i]], [_[1] for _ in edge_atoms[i]], fig=fig, ax=ax)
#     plot_edges([my_edges[i][0]], fig=fig, ax=ax)
#     plot_verts([_[0] for _ in my_edge_verts[i]], [_[1] for _ in my_edge_verts[i]], fig=fig, ax=ax, Show=True)

# Calculate the surfaces
my_surf = build_surf(locs=[_[0] for _ in my_surf_atoms], rads=[_[1] for _ in my_surf_atoms],
                     epnts=[_[0] for _ in my_edges], res=res, net_type=net_type)

# Plot everything
# Create the figure
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Plot the surfaces
plot_surfs([my_surf[0]], [my_surf[1]], fig=fig, ax=ax, alpha=0.5)
# Plot the edges
if show_edges:
    plot_edges([_[0] for _ in my_edges], fig=fig, ax=ax, colors=['k'] * 4, thickness=2)
# Plot the vertices
if show_verts:
    plot_verts([_[0] for _ in my_verts], [_[1] for _ in my_verts], fig=fig, ax=ax, colors=['r'] * 4)
# Plot the atoms
plot_balls(alocs=[_[0] for _ in my_surf_atoms], arads=[_[1] for _ in my_surf_atoms], alpha=0.2, fig=fig, ax=ax, res=10)

# Set the scales for the figure
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)

# Set the title for the figure
ax.set_title(title, font=dict(size=20, family='serif'))
plt.show()



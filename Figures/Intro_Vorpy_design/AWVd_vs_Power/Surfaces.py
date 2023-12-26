from numpy import array as ar
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_atoms, plot_verts, plot_edges, plot_surfs
from System.Network.verts.calc_vert import calc_vert
from System.Network.edges.build_edge import build_edge
from System.Network.surfs.build_surf import build_surf

# Set the radii and distances for the surrounding atoms
r, d = 0.5, 25.0
my_vert_atoms = [(ar([0.0, d, d]), r), (ar([0.0, -d, d]), r), (ar([0.0, -d, -d]), r), (ar([0.0, d, -d]), r)]
# Create the Surface Atoms
my_surf_atoms = [(ar([2.5, 0.0, 0.0]), 0.5), (ar([-2.5, 0.0, 0.0]), 2.5)]

# Calculate the vertices
vert_atoms = [[_ for _ in my_surf_atoms] + [my_vert_atoms[i], my_vert_atoms[(i+1) % 4]] for i in range(4)]
my_verts = [calc_vert(ar([_[0] for _ in my_atoms]), ar([_[1] for _ in my_atoms])) for my_atoms in vert_atoms]

# ## Plot the verts and their atoms
# for i in range(4):
#     fig = plt.figure()
#     ax = fig.add_subplot(projection='3d')
#     plot_atoms([_[0] for _ in vert_atoms[i]], [_[1] for _ in vert_atoms[i]], fig=fig, ax=ax)
#     plot_verts([my_verts[i][0]], [my_verts[i][1]], fig=fig, ax=ax, Show=True, spheres=True)

# Calculate the Edges
edge_atoms = [[_ for _ in my_surf_atoms] + [my_vert_atoms[(i+1)%4]] for i in range(4)]
my_edge_verts = [(my_verts[j], my_verts[(j+1)%4]) for j in range(4)]
my_edges = [build_edge(alocs=ar([_[0] for _ in edge_atoms[i]]), arads=ar([_[1] for _ in edge_atoms[i]]),
                       vlocs=ar([_[0] for _ in my_edge_verts[i]]), res=0.5) for i in range(4)]
# Plot the edges and their atoms
# for i in range(4):
#     fig = plt.figure()
#     ax = fig.add_subplot(projection='3d')
#     plot_atoms([_[0] for _ in edge_atoms[i]], [_[1] for _ in edge_atoms[i]], fig=fig, ax=ax)
#     plot_edges([my_edges[i][0]], fig=fig, ax=ax)
#     plot_verts([_[0] for _ in my_edge_verts[i]], [_[1] for _ in my_edge_verts[i]], fig=fig, ax=ax, Show=True)

# Calculate the surfaces
my_surf = build_surf([_[0] for _ in my_surf_atoms], [_[1] for _ in my_surf_atoms], [_[0] for _ in my_edges], 0.1, 'vor')

# Plot everything
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_surfs([my_surf[0]], [my_surf[1]], fig=fig, ax=ax, alpha=0.2)
plot_edges([_[0] for _ in my_edges], fig=fig, ax=ax)
plot_verts([_[0] for _ in my_verts], [_[1] for _ in my_verts], fig=fig, ax=ax)
plot_atoms(alocs=[_[0] for _ in my_surf_atoms], arads=[_[1] for _ in my_surf_atoms], fig=fig, ax=ax, Show=True)



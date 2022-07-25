import os
from System.system import System
from System.Network.network import Network
from Presentation.Visualize.visualize import plot_edges
import matplotlib.pyplot as plt
from System.Network.net_funcs import find_vertices
os.chdir("../..")


# Files
m_file = "./Data/test_data/Na_W_cluster5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

# Get the System
sys = System(m_file)


find_vertices(sys.net)
#
# sys.net.verts = []
# sys.net.edges = []
# # Add voronota data
# sys.add_vta_data(b_file, v_file)
#
# # Get the new locations for the vertices
# for i in range(len(sys.net.verts)):
#     sys.net.verts[i] = calc_vert(sys.net.verts[i].atoms)
#
#
verts = []
b = sys.calc_box(1)
for vert in sys.net.verts:
    if b[0][0] < vert.loc[0] < b[1][0] and b[0][1] < vert.loc[1] < b[1][1] and b[0][2] < vert.loc[2] < b[1][2]:
        verts.append(vert)
sys.net.verts = verts

# # Connect the network
sys.net.connect()
# # Build the meshes
sys.net.build_meshes(min_dist=0.5)
# Print the data
# print("\n\nAtom Data: ###################################################################\n")
# for i in range(len(sys.net.atoms)):
#     print("\nAtom: {}".format(i))
#     print("Vertices:", sys.atoms[i].verts)
#     print("Edges:", sys.atoms[i].edges)
#     print("Surfaces:", sys.atoms[i].surfs)
#
# print("\n\nVertex Data: ###################################################################\n")
# for i in range(len(sys.net.verts)):
#     print("\nVertex: {}".format(i))
#     print("Atoms:", sys.net.verts[i].atoms)
#     print("Edges:", sys.net.verts[i].edges)
#     print("Surfaces:", sys.net.verts[i].surfs)
#
# print("\n\nEdge Data: ###################################################################\n")
# for i in range(len(sys.net.edges)):
#     print("\nEdge: {}".format(i))
#     print("Atoms:", sys.net.edges[i].atoms)
#     print("Vertices:", sys.net.edges[i].verts)
#     print("Surfaces:", sys.net.edges[i].surfs)
#
# print("\n\nSurface Data: ###################################################################\n")
# for i in range(len(sys.net.surfs)):
#     print("\nSurface: {}".format(i))
#     print("Atoms:", sys.net.surfs[i].atoms)
#     print("Vertices:", sys.net.surfs[i].verts)
#     print("Edges:", sys.net.surfs[i].edges)

# myVert = find_vertices(sys)
# # Plot the System
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# plot_atoms(sys.atoms[:15], fig=fig, ax=ax, alpha=.1, colors=['w' for i in range(len(sys.atoms))])
# plot_verts(sys.atoms[0].verts, fig=fig, ax=ax, colors=['r' for i in range(len(sys.net.verts))])
plot_edges(sys.net.edges, fig=fig, ax=ax, Show=True)
# plot_surfs(sys.atoms[0].surfs, fig=fig, ax=ax, alpha=1, Show=True)

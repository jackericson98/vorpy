import os
from system.objects import System
from build_network.calculators import calc_dist
from Presentation.Visualize.visualize import plot_atoms, plot_verts
import matplotlib.pyplot as plt
from build_network.find_vertices import find_network
os.chdir("../..")


# Files
m_file = "./Data/test_data/Na_W_cluster5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

# Get the system
sys = System(m_file)
find_network(sys)
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
# # Connect the network
# sys = connect_network(sys)
#
# # Build the meshes
# build_meshes(sys, min_dist=0.5)
for atom in sys.net.verts[0].atoms:
    print(calc_dist(atom.loc, sys.net.verts[0].loc) - (atom.rad + sys.net.verts[0].rad))
# myVert = find_vertices(sys)
# # Plot the system
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_atoms(sys.atoms, fig=fig, ax=ax, alpha=0.1)
plot_verts(sys.net.verts, plot_spheres=True, fig=fig, ax=ax, Show=True, colors=['w'] + ['r' for i in range(len(sys.net.verts))])
# plot_edges(sys.net.atoms[0].edges, fig=fig, ax=ax)
# plot_surfs(sys.net.atoms[0].surfs, fig=fig, ax=ax, alpha=1)

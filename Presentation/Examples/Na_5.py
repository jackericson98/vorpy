import os
from System.system import System
from Network.calculators import calc_dist
from Presentation.Visualize.visualize import plot_atoms, plot_verts, plot_edges, plot_surfs
import matplotlib.pyplot as plt
from Network.find_vertices import find_network
from Network.build_network import build_network
from Network.connect_network import connect_network
from Cells.build_mesh import build_meshes
os.chdir("../..")


# Files
m_file = "./Data/test_data/Na_W_cluster5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

# Get the System
sys = System(m_file)

for atom in sys.atoms:
    atom.rad = atom.rad/10

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
sys = connect_network(sys)
# # Build the meshes
build_meshes(sys, min_dist=0.5)

# myVert = find_vertices(sys)
# # Plot the System
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_atoms(sys.atoms[:10], fig=fig, ax=ax, alpha=1, colors=['w' for i in range(len(sys.atoms))])
plot_verts(sys.atoms[0].verts, fig=fig, ax=ax, colors=['r' for i in range(len(sys.net.verts))])
plot_edges(sys.net.atoms[0].edges, fig=fig, ax=ax)
plot_surfs(sys.net.atoms[0].surfs, fig=fig, ax=ax, alpha=1, Show=True)

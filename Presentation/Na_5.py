import os
from objects import System
from build_network import build_network
from build_mesh import build_meshes
from visualize import plot_atoms, plot_verts, plot_surfs, plot_edges
from connect_network import connect_network
import matplotlib.pyplot as plt
os.chdir("..")

# Files
m_file = "./Data/test_data/Na_W_cluster5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

# Get the system
sys = System(m_file)
sys.add_vta_data(b_file, v_file)
sys = connect_network(sys)

# Build the meshes
build_meshes(sys, min_dist=0.5)

# Plot the system
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_atoms(sys.net.atoms[:1], fig=fig, ax=ax, alpha=0.1)
plot_verts(sys.net.atoms[0].verts, fig=fig, ax=ax)
plot_edges(sys.net.atoms[0].edges, fig=fig, ax=ax)
plot_surfs(sys.net.atoms[0].surfs, fig=fig, ax=ax, alpha=1, Show=True)
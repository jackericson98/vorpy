import os
from System.system import System
from System.sys_funcs import build_vrnta_meshes
from Visualize.visualize import *

os.chdir("../..")
# Files
m_file = "./Data/test_data/Na_W_cluster5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

# Get the System
sys = System(m_file)

# Add the voronota data
sys.add_vta_data(b_file, v_file)
# Connect the network
sys.net.connect()

build_vrnta_meshes(sys, 20)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
# plot_atoms(sys.net.atoms, fig=fig, ax=ax, dfo=20)
plot_edges(sys.net.edges, fig=fig, ax=ax, dfo=20)
plot_verts(sys.net.verts, fig=fig, ax=ax)
plot_surfs(sys.net.surfs, simps=True, fig=fig, ax=ax, Show=True)

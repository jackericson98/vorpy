import os
from System.system import System
from Visualize.visualize import *
from System.input import add_vta_data

os.chdir("../..")
# Files
m_file = "./Data/test_data/Na_5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

# Get the System
sys = System(m_file)

# Add the voronota data
add_vta_data(sys, b_file, v_file)

# Build the network
sys.net.build_surfs()

atoms, edges, surfs, verts = [],[], [], []
for atom in sys.net.atoms:
    cell = True
    for vert in atom.verts:
        if len(vert.edges) < 3:
            cell = False
            break
    if cell:
        atoms.append(atom)
        edges += atom.edges
        verts += atom.verts
        surfs += atom.surfs

# Plot everything
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
# plot_atoms(sys.net.atoms, fig=fig, ax=ax, dfo=20)
plot_edges(edges, fig=fig, ax=ax)
plot_verts(verts, fig=fig, ax=ax)
plot_surfs(surfs, simps=True, fig=fig, ax=ax, Show=True)

import os
from System.system import System
from Visualize.visualize import plot_edges, plot_verts, plot_atoms, plot_surfs
import matplotlib.pyplot as plt
os.chdir("../..")


# Files
m_file = "./Data/test_data/Na_W_cluster5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

# Get the System
sys = System(m_file)
sys.build_network(1)

# Plot the System
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_atoms(sys.atoms[:1], fig=fig, ax=ax, alpha=.1, colors=['w' for i in range(len(sys.atoms))])
plot_verts(sys.atoms[0].verts, fig=fig, ax=ax, colors=['r' for i in range(len(sys.net.verts))])
plot_edges(sys.atoms[0].edges, fig=fig, ax=ax, Show=True)
# plot_surfs(sys.atoms[0].surfs, fig=fig, ax=ax, alpha=1, simps=True, Show=True)

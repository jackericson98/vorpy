import os
from System.system import System
from Visualize.visualize import plot_edges, plot_verts
import matplotlib.pyplot as plt
os.chdir("../..")


# Files
m_file = os.getcwd() + "./Data/test_data/Na_W_cluster5.pdb"
v_file = os.getcwd() + "./Data/test_data/Na5_verts1.txt"

# Get the System
sys = System()
sys.load_sys(file="C:/Users/jacke/PycharmProjects/vorpy/Data/test_data/EDTA_Mg.pdb")
sys.build_network(sol_verts=True, min_dist=0.5, max_vert=20)


# # Plot the System
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# # plot_atoms(sys.atoms[:1], fig=fig, ax=ax, alpha=.1, colors=['w' for i in range(len(sys.atoms))])
# plot_verts(sys.net.verts, fig=fig, ax=ax, colors=['r' for i in range(len(sys.net.verts))])
# plot_edges(sys.net.edges, fig=fig, ax=ax)
# # plot_surfs(sys.atoms[0].surfs, fig=fig, ax=ax, alpha=1, simps=True)
# plt.show()

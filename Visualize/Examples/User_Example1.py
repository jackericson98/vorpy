from Visualize.mpl_visualize import *
import os
from System.system import System
os.chdir("../..")

# Set up a list of files to analyze
files = ["DB1976.pdb"]
my_dir = os.getcwd()
# Analyze each of the files
for file in files:

    os.chdir(my_dir)
    # Get the System
    sys = System(os.getcwd() + "/Data/test_data/" + file)
    sys.build_network(sol_verts=False, surf_res=0.2, max_vert=5, box_size=1.5)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    plot_atoms(sys.atoms[0:4], fig=fig, ax=ax)

    plot_verts(sys.net.verts, fig=fig, ax=ax)
    plot_edges(sys.net.edges, fig=fig, ax=ax, Show=True)


    print(sys.net.doublets)

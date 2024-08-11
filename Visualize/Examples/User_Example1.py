from Visualize.mpl_visualize import *
import os
from System.system import System
os.chdir("../..")

# Set up a list of files to analyze
files = ["Complex1_frame1.pdb"]
my_dir = os.getcwd()
# Analyze each of the files
for file in files:

    os.chdir(my_dir)
    # Get the System
    sys = System(os.getcwd() + "/Data/test_data/" + file)
    sys.build_network(surf_res=0.1, max_vert=20, box_size=2.5)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    plot_balls(sys.atoms[0:4], fig=fig, ax=ax)

    plot_verts(sys.net.verts, fig=fig, ax=ax)
    plot_edges(sys.net.edges, fig=fig, ax=ax, Show=True)


    print(sys.net.doublets)

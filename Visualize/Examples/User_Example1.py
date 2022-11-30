from Visualize.mpl_visualize import *
import os
from System.system import System
os.chdir("../..")

# Set up a list of files to analyze
files = ["Na5.pdb"]
my_dir = os.getcwd()
# Analyze each of the files
for file in files:
    os.chdir(my_dir)
    # Get the System
    sys = System()

    sys.load_sys(file=os.getcwd() + "/Data/test_data/" + file)
    sys.build_network(sol_verts=True, surf_res=0.2, max_vert=10, box_size=2)
    plot_edges(sys.net.edges, Show=True)


    print(sys.name + " Completed")

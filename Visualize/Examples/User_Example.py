import os
from System.system import System
from Visualize.visualize import plot_edges, plot_verts
import matplotlib.pyplot as plt
os.chdir("../..")


# Files
m_file = os.getcwd() + "./Data/test_data/Na_W_cluster5.pdb"
v_file = os.getcwd() + "./Data/test_data/Na5_verts1.txt"




files = ["DB1976.pdb", "Complex1_frame1.pdb"]
my_dir = os.getcwd()
for file in files:
    os.chdir(my_dir)
    # Get the System
    sys = System()

    sys.load_sys(file=os.getcwd() + "/Data/test_data/" + file)
    sys.build_network(sol_verts=True, min_dist=0.5, max_vert=10, box_size=2)



    print(sys.name + " Completed")
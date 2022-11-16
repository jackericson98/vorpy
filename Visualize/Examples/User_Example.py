import os
from System.system import System, Network, Vertex
from Visualize.visualize import plot_edges, plot_verts, plot_atoms
import matplotlib.pyplot as plt
os.chdir("../..")


# Files
m_file = "C:/Users/jacke/PycharmProjects/vorpy/Data/test_data/Na5.pdb"
v_file = "C:/Users/jacke/PycharmProjects/vorpy/Data/test_data/Na_W_cluster5_vertices.txt"
b_file = "C:/Users/jacke/PycharmProjects/vorpy/Data/test_data/Na_W_cluster5_balls.txt"


# Get the System
loaded_sys = System(base_file=m_file)
loaded_sys.load_verts(vert_file=v_file, vta_ball_file=b_file)

# built_sys = System(base_file=m_file)
# built_sys.build_network(min_dist=0.5, max_vert=10, box_size=2, flat_faces=True)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
plot_verts([loaded_sys.net.verts[2]], True, fig=fig, ax=ax)

new_vert = Vertex(loaded_sys.net.verts[2].atoms)
new_vert.calc_ff_vert()
plot_verts([new_vert], fig=fig, ax=ax, colors=['r'], spheres=True)
plot_atoms(loaded_sys.net.verts[2].atoms, fig=fig, ax=ax, Show=True)
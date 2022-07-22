import os
from System.system import System, Surface
from Network.connect_network import connect_network
from Meshes.build_meshes import make_mesh
from Presentation.Visualize.visualize import *
import numpy as np

os.chdir("../..")
# Files
m_file = "./Data/test_data/Complex1.pdb"
b_file = "./Data/test_data/Complex1_balls.txt"
v_file = "./Data/test_data/Complex1_vertices.txt"

# Get the System
sys = System(m_file)

# Add the voronota data
sys.add_vta_data(b_file, v_file)
# Connect the network
connect_network(sys)


def build_vrnta_meshes(sys, min_dist):
    # Add the edge points for Voronota's vertices
    for edge in sys.net.edges:
        # Get the vertices of the edge
        v0, v1 = edge.verts
        # Define the vector from v0 to v1
        r01 = np.array(v1.loc) - np.array(v0.loc)
        # Get the magnitude of that vector
        mag = np.linalg.norm(r01)
        # Find the unit vector in the direction from v0 to v1
        r0_hat = r01/mag
        # Figure out the number of steps to take along the edge based on the minimum distance
        num_steps = max(int(mag/min_dist), 10)
        # Define the singular vector step
        step = r0_hat*(mag/num_steps)
        # Add the first vertex to the edge's points
        edge.points = [v0.loc]
        # # Go through step by step adding points until we reach the other vertex
        for i in range(num_steps):
            p_1 = edge.points[-1]
            edge.points.append([p_1[0] + step[0], p_1[1] + step[1], p_1[2] + step[2]])
        # Add v1's location
        edge.points.append(v1.loc)
    # Build each mesh
    for i in range(len(sys.net.surfs)):
        # Pull the surf out and name it better
        surf = sys.net.surfs[i]
        # Add the edge points to the surfaces edge points
        for edge in surf.edges:
            surf.edge_points += edge.points
        # Get the atoms and make sure the smaller one is a0
        a0, a1 = sys.net.surfs[i].atoms
        if a0.rad > a1.rad:
            a0, a1 = a1, a0
        # Get the normal vector along the direction from a1 to a0
        r10 = np.array(a0.loc) - np.array(a1.loc)
        r10_hat = r10/np.linalg.norm(r10)
        # Get the difference in the radii
        r_diff = a1.rad - a0.rad
        # Change the radius and location of the smaller atom to represent a flat surface
        a0.loc = [a0.loc[0] + r_diff * r10_hat[0], a0.loc[1] + r_diff * r10_hat[1], a0.loc[2] + r_diff * r10_hat[2]]
        a0.rad = a1.rad
        surf.atoms = [a0, a1]
        # Make the mesh
        make_mesh(surf, 0.1, vta=True)
        # Return a0 to its radius and location
        a0.rad = a0.rad - r_diff
        a0.loc = [a0.loc[0] - r_diff * r10_hat[0], a0.loc[1] - r_diff * r10_hat[1], a0.loc[2] - r_diff * r10_hat[2]]
    return sys


build_vrnta_meshes(sys, 20)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
# plot_atoms(sys.net.atoms, fig=fig, ax=ax, dfo=20)
plot_edges(sys.net.edges, fig=fig, ax=ax, dfo=20)
plot_verts(sys.net.verts, fig=fig, ax=ax)
plot_surfs(sys.net.surfs, simps=True, fig=fig, ax=ax, Show=True)

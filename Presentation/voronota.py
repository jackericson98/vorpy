from objects import Vertex, System, Edge, Surface
import os
from visualize import plot_atoms, plot_verts, plot_surfs, plot_edges
import matplotlib.pyplot as plt
from build_mesh import build_meshes, calc_surf
os.chdir("..")


# Check surf function. Takes in a set of atoms and a list of surfs and returns the corresponding surf or None if no surf
def check_surf(s_atoms, surf_list):
    # Go through each surf in the surf list
    for surf in surf_list:
        # Check if the given atoms correspond to the atoms in the surf
        if s_atoms.issubset(surf.atoms):
            # Return the surf
            return surf
    return


# Check edge function. Takes in a set of atoms and a list of edges and returns the corresponding edge or None if no edge
def check_edge(e_atoms, edge_list):
    # Go through each edge in the edge list
    for edge in edge_list:
        # Check if the given atoms correspond to the atoms in the edge
        if e_atoms.issubset(edge.atoms):
            # Return the edge
            return edge
    return


# Check vert function. Takes in a set of atoms and a list of verts and returns the corresponding edge or None if no vert
def check_vert(v_atoms, vert_list):
    # Go through each edge in the edge list
    for vert in vert_list:
        # Check if the given atoms correspond to the atoms in the edge
        if v_atoms.issubset(vert.atoms):
            # Return the edge
            return vert
    return


# Build voronota system function. Takes in voronota data and returns a system
def build_vta_sys(mol_file, ball_file, vert_file):
    # Create the system and load the files
    sys = System(mol_file)
    vert_file = open(vert_file).readlines()
    ball_file = open(ball_file).readlines()
    # Interpret the balls
    balls = []
    for i in range(len(ball_file)):
        # Split the data
        data = ball_file[i].split(" ")
        # Grab the data reference for the atoms
        balls.append(sys.atoms[int(data[5])])
    # Interpret the vertices
    for i in range(len(vert_file)):
        # Split the data
        data = vert_file[i].split(" ")
        # Add the vertex data
        loc, rad = [float(data[4]), float(data[5]), float(data[6])], float(data[7])
        atoms = [balls[int(data[0])], balls[int(data[1])], balls[int(data[2])], balls[int(data[3])]]
        myVert = Vertex(loc, rad, atoms=atoms)
        sys.net.verts.append(myVert)

    # Create edges and add connections between verts and edges
    # Go through each vertex and find its edges
    for vert1 in sys.net.verts:
        # Check every combination of vert atoms as an edge
        for i in range(4):
            # Grab the atoms
            atoms = {vert1.atoms[i], vert1.atoms[(i + 1) % 4], vert1.atoms[(i + 2) % 4]}
            verts = []
            # Find the possible verts
            for vert2 in sys.net.verts:
                if atoms.issubset(vert2.atoms):
                    verts.append(vert2)
            # Find which edge, if any, go nowhere
            if len(verts) < 2:
                continue
            # Check to see if the edge has been found
            my_edge = check_edge(atoms, sys.net.edges)
            if my_edge is None:
                # Create the edge
                my_edge = Edge(list(atoms), verts)
                # Add the edge to the system
                sys.net.edges.append(my_edge)
                # Add the edge to the verts
                verts[0].edges.append(my_edge)
                verts[1].edges.append(my_edge)

    # Create surfaces and add connections for edges and verts
    for vert1 in sys.net.verts:
        # Go through each combination of sets atom in the vertices' atom list
        t_ndxs = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        for ndxs in t_ndxs:
            # Grab the atoms
            t_atoms = {vert1.atoms[ndxs[0]], vert1.atoms[ndxs[1]]}
            # Check to see if we have recorded this surface before
            if check_surf(t_atoms, sys.net.surfs):
                continue
            # Put together a list of edges that have our atoms
            edges = []
            for edge in sys.net.edges:
                if t_atoms.issubset(edge.atoms):
                    edges.append(edge)
            # Put together a list of verts that have our atoms
            verts = []
            for vert2 in sys.net.verts:
                if t_atoms.issubset(vert2.atoms):
                    verts.append(vert2)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(verts) == len(edges):
                my_surf = Surface(list(t_atoms), calc_surf(list(t_atoms)), verts=verts, edges=edges)
                sys.net.surfs.append(my_surf)
                list(t_atoms)[0].surfs.append(my_surf)
                list(t_atoms)[1].surfs.append(my_surf)
                list(t_atoms)[0].edges += edges
                list(t_atoms)[1].edges += edges
                list(t_atoms)[0].verts += verts
                list(t_atoms)[1].verts += verts
    # Return the system we have created
    return sys


m_file = "./Data/test_data/Na_W_cluster5.pdb"
b_file = "./Data/test_data/Na_W_cluster5_balls.txt"
v_file = "./Data/test_data/Na_W_cluster5_vertices.txt"

sys = build_vta_sys(m_file, b_file, v_file)
build_meshes(sys, min_dist=0.1)
net = sys.net
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
print(net.atoms[0].verts)
# plot_atoms(net.atoms[:20], fig=fig, ax=ax, alpha=0.1)
plot_verts(net.atoms[0].verts, fig=fig, ax=ax)
plot_edges(net.atoms[0].edges[:1], fig=fig, ax=ax, Show=True)
# plot_surfs(net.atoms[0].surfs, fig=fig, ax=ax, alpha=1, Show=True)

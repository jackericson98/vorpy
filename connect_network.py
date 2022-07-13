from objects import Vertex, System, Edge, Surface
from visualize import plot_atoms, plot_verts, plot_surfs, plot_edges
import matplotlib.pyplot as plt
from build_mesh import build_meshes, calc_surf
from build_network import build_network


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


# Connect network function. Takes in a network with
def connect_network(sys):
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

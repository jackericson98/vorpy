from System.system import Edge, Surface
from Network.find_vertices import *


# Connect network function. Takes in a network with vertices defined and returns a filled out network of atoms,
# surfaces, edges and vertices
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
            if len(verts) == 1:
                continue
            # Check to see if the edge has been found
            my_edge = check_edge(atoms, sys.net.edges)
            if my_edge is None:
                # Create the edge
                my_edge = Edge(list(atoms), verts)
                # Add the edge to the System
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
                my_surf = Surface(list(t_atoms), verts=verts, edges=edges)
                sys.net.surfs.append(my_surf)
                list(t_atoms)[0].surfs.append(my_surf)
                list(t_atoms)[1].surfs.append(my_surf)
                list(t_atoms)[0].edges += edges
                list(t_atoms)[1].edges += edges
                list(t_atoms)[0].verts += verts
                list(t_atoms)[1].verts += verts

    # Add the surfaces to the edges
    for edge in sys.net.edges:
        edge.surfs = []
        for surf in sys.net.surfs:
            if set(surf.atoms).issubset(edge.atoms):
                edge.surfs.append(surf)
    # Add the surfaces to the vertices
    for vert in sys.net.verts:
        vert.surfs = []
        for surf in sys.net.surfs:
            if set(surf.atoms).issubset(vert.atoms):
                vert.surfs.append(surf)

    # Return the System we have created
    return sys

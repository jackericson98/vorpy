from System.calcs import *
from System.Network.edge import Edge
from System.Network.surface import Surface


# Connect network method.
def connect(net):
    # Create the edges
    for vert1 in net.verts:
        # Check every combination of vert atoms as an edge
        for i in range(4):
            # Grab the atoms
            atoms = {vert1.atoms[i], vert1.atoms[(i + 1) % 4], vert1.atoms[(i + 2) % 4]}
            # If the edge has been found before, continue
            if check_edge(atoms, net.edges):
                continue
            verts = []
            # Find the possible verts (the original vert and the new vert)
            for vert2 in net.verts:
                if atoms.issubset(vert2.atoms):
                    verts.append(vert2)
            # If the number of valid vertices for the edge is 1
            if len(verts) == 1:
                continue
            # Create the edge
            my_edge = Edge(list(atoms), verts, net)
            # Add the edge to the System
            net.edges.append(my_edge)

    # Create the surfaces
    net.surfs = []
    for edge1 in net.edges:
        # Go through the edge's atoms combinations
        for i in range(3):
            atoms = {edge1.atoms[i], edge1.atoms[(i + 1) % 3]}
            # If the surface has been found before continue
            if check_surf(atoms, net.surfs):
                continue
            # Put together a list of edges that have our atoms
            edges = []
            for edge2 in net.edges:
                if atoms.issubset(edge2.atoms):
                    edges.append(edge2)
            # Put together a list of verts that have our atoms
            verts = []
            for vert2 in net.verts:
                if atoms.issubset(vert2.atoms):
                    verts.append(vert2)
            # In order to be a true surface the number of edges need to be equal to the number of verts
            if len(verts) == len(edges):
                my_surf = Surface(list(atoms), verts=verts, net=net, edges=edges)
                net.surfs.append(my_surf)
            else:
                pass

    # Add the vertices, edges and surfs to the atoms
    for atom in net.atoms:
        # Reset the atom's vert list
        atom.verts = []
        # Go through the verts in the network
        for vert in net.verts:
            # If the atom is in the vertices atoms add the vertex to the atom's list of vertices
            if {atom}.issubset(vert.atoms):
                atom.verts.append(vert)
        # Reset the atom's edge list
        atom.edges = []
        # Go through the edges in the network
        for edge in net.edges:
            # If the atom is in the edge's list of atoms add the edge to the atoms list of edges
            if {atom}.issubset(edge.atoms):
                atom.edges.append(edge)
        # Reset the atom's surf list
        atom.surfs = []
        # Go through the surfs in the network
        for surf in net.surfs:
            # If the atom is in the surfs list of atoms add the surf to the atoms list of surfs
            if {atom}.issubset(surf.atoms):
                atom.surfs.append(surf)

    # Add the edges and surfs to the vertices
    for vert in net.verts:
        # Reset the vertexes edge list
        vert.edges = []
        # Go through the edges in the network
        for edge in net.edges:
            # If the edges atoms are in the vertices atoms add it to the vertex
            if set(edge.atoms).issubset(vert.atoms):
                vert.edges.append(edge)
        # Reset the vertexes surf list
        vert.surfs = []
        # Go through the surfaces in the network
        for surf in net.surfs:
            # If the surfaces atoms are in the vertexes atoms add it to the vertex
            if set(surf.atoms).issubset(vert.atoms):
                vert.surfs.append(surf)

    # Add the surfs to the edges
    for edge in net.edges:
        # Reset the edges surf list
        edge.surfs = []
        # Go through the surfaces in the network
        for surf in net.surfs:
            # If the surfaces atoms are in the edges atoms add it to the edge
            if set(surf.atoms).issubset(edge.atoms):
                edge.surfs.append(surf)

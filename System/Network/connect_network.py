from System.calcs import *
from System.Network.edge import Edge
from System.Network.surface import Surface


# Connect network function. Takes in a broken network of vertices and spits out a connected network
def connect(net):

    ################################################# Create the edges #################################################

    # Reset the network's list of edges
    net.edges = []
    net.surfs = []
    # Go through the vertices in the network searching for potential edges
    for vert1 in net.verts:
        # Set up the edges and doublet indices
        my_edges = []
        dub_edges = []
        # Check to see if the vertex is a doublet
        if vert1.doublet:
            # If it is, we need to find what type of doublet it is (i.e. the # of edges)
            for i in range(4):
                # Get each combination of atoms
                atoms = [vert1.atoms[i], vert1.atoms[(i + 1) % 4], vert1.atoms[(i + 2) % 4]]
                # Calculate the inscribed circle for the current set of test atoms
                circ = calc_circ(atoms)
                # If this circle doesn't overlap with the other atom this is doublet edge
                if calc_dist(circ[0], vert1.atoms[(i + 3) % 4].loc) > circ[1] + vert1.atoms[(i + 3) % 4].rad:
                    # Add the
                    my_edges.append(Edge(atoms, [vert1], net, doublet=True))
                    dub_edges.append(i)
            # Mark the vertex as a 1 or 3 surface doublet type
            if len(dub_edges) == 2:
                # Set the vertex doublet type
                vert1.d_type = "1"
                # Get the atoms in the surface
                surf_atoms = [atom for atom in my_edges[0].atoms if atom in my_edges[1].atoms]
                # Create the surface
                mySurf = Surface(surf_atoms, edges=my_edges)
                net.surfs.append(mySurf)
            # If the vertex is a 3 surface type create those surfaces
            elif len(dub_edges) == 3:
                # Set the vertex doublet type
                vert1.d_type = "3"
                # Get the center atom
                min_rad, myAtom = np.inf, None
                # Find the smallest atom in the bunch
                for atom1 in vert1.atoms:
                    if atom1.rad < min_rad:
                        myAtom = atom1
                        min_rad = myAtom.rad
                # Get a list of the other atoms in the surfaces
                other_atoms = [atom for atom in vert1.atoms if atom != myAtom]
                # Create the surfaces
                vert_surfs = []
                edges = my_edges
                for k in range(3):
                    surf_edges = [edges[k], edges[(k + 1) % 3]]
                    surf_atoms = [atom for atom in vert1.atoms if atom in surf_edges[0].atoms and atom in surf_edges[1].atoms]
                    vert_surfs.append(Surface(surf_atoms, net, surf_edges))
            elif len(dub_edges) == 3:
                vert1.d_type = "3"
        # Check every combination of vert atoms as an edge
        for i in range(4):
            # Check for previously made doublet edges
            if i in dub_edges:
                continue
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
            my_edges.append(my_edge)
        # Add the edges to the network's list of edges
        net.edges += my_edges

    ################################################### Create the surfaces ############################################

    # Reset the network's list of surfaces
    net.surfs = []
    # Go through the edges in the network
    for edge1 in net.edges:
        # If the edge is a doublet, skip it
        if edge1.doublet:
            continue
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

    ################################################# Connect the atoms ################################################

    # Go through the atoms in the network adding vertices, edges and surfaces
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

    ################################################# Connect the vertices #############################################

    # Go through the vertices in the network adding the edges and surfaces that share atoms
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

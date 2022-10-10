from System.calcs import *
from System.Network.edge import Edge
from System.Network.surface import Surface


############################################ Filter System #############################################################


# Filter vertices function. Filters out vertices that are repeated, outside the box or larger than the max vertex value
def filter_verts(net):

    # Set up a list of vertex ndxs and vertices
    vert_ndxs = []
    verts = []
    # Go through the vertices
    for i in range(len(net.verts)):

        # Boolean for whether the vertex is inside the box or not
        loc_in_box = True
        vert = net.verts[i]
        # Check for None vertices
        if vert.loc is None:
            continue

        # Doublet verts
        if vert.doublet:
            loc2_in_box = True
            # Check if the vertex is inside the box
            for j in range(3):
                if vert.loc[j] < net.box[0][j] or vert.loc[j] > net.box[1][j]:
                    loc_in_box = False
                if net.verts[i].loc2[j] < net.box[0][j] or net.verts[i].loc2[j] > net.box[1][j]:
                    loc2_in_box = False
            # If loc is outside we replace the first vertex and make the vertex a non doublet
            if not loc_in_box and loc2_in_box:
                vert.loc, vert.rad = vert.loc2, vert.rad2
                vert.loc2, vert.rad2 = None, None
                vert.doublet = False
            # If just the second vertex is outside make the loc2, rad2 values None and make the vertex non-doublet
            elif not loc2_in_box and loc_in_box:
                vert.loc2, vert.rad2 = None, None
                vert.doublet = False
            # If both vertices are outside make both locations None and make it a non doublet
            elif not loc_in_box and not loc2_in_box:
                vert.loc, vert.rad = None, None
                vert.loc2, vert.rad2 = None, None
                vert.doublet = False
        else:
            # Check if the vertex is inside the box
            for j in range(3):
                if vert.loc[j] < net.box[0][j] or vert.loc[j] > net.box[1][j]:
                    loc_in_box = False

        # Search the list of vertices for the vertex
        if vert.ndx not in vert_ndxs and loc_in_box:
            vert_ndxs.append(net.verts[i].ndx)
            verts.append(net.verts[i])

    # Set the networks vertices
    net.verts = verts


############################################## Doublets ################################################################

# Doublet creating function. Fills in the surfaces and edges inside the doublet vertex
def doublet(vert, net):

    # Set up a variable for doublet edges
    dub_edges = []
    my_edges = []

    # Find what type of doublet it is (i.e. the # of edges) by counting the number of "free" inscribed circles
    for i in range(4):
        # Get the current combination of atoms to test doubletness of the edge
        atoms = [vert.atoms[i], vert.atoms[(i + 1) % 4], vert.atoms[(i + 2) % 4]]
        # Calculate the inscribed circle for the current set of test atoms
        circ = calc_circ(atoms)
        # If this circle doesn't overlap with the other atom this is doublet edge
        if calc_dist(circ[0], vert.atoms[(i + 3) % 4].loc) > circ[1] + vert.atoms[(i + 3) % 4].rad:
            # Add the doublet edge to the list of edges
            my_edges.append(Edge(atoms, net, [vert], doublet=True))
            dub_edges.append(i)

    # Add the edges to the network
    net.edges += my_edges

    # If there are 2 edges involved in the doublet it is a type 1 doublet and has 1 surface
    if len(dub_edges) == 2:
        # Set the vertex doublet type
        vert.d_type = "1"
        # Get the atoms in the surface
        surf_atoms = [atom for atom in my_edges[0].atoms if atom in my_edges[1].atoms]
        # Create the surface
        mySurf = Surface(surf_atoms, net, edges=my_edges, doublet=True)
        # Add the surface to the network
        net.surfs.append(mySurf)

    # If there are 3 edges involved in the doublet it is a type 3 doublet and has 3 surfaces
    elif len(dub_edges) == 3:
        # Set the vertex doublet type
        vert.d_type = "3"
        # Create the surfaces for the vertex surfaces
        for k in range(3):
            edges = [my_edges[k], my_edges[(k + 1) % 3]]
            atoms = [atom for atom in edges[0].atoms if atom in edges[1].atoms]
            net.surfs.append(Surface(atoms, net, edges, doublet=True))


# Make objects function. Checks the atoms of the vertices for patterns and creates edges and surfaces
def make_objects(net):

    # Reset the network's list of edges and surfaces for a clean slate
    net.edges = []
    net.surfs = []

    ################################################# Create the edges #################################################

    # Go through the vertices in the network searching for potential edges
    for vert1 in net.verts:
        # Set up the edges and doublet indices
        my_edges = []
        # If the vertex is a doublet, fill in the insides of the doublet
        if vert1.doublet:
            doublet(vert1, net)
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
            my_edge = Edge(list(atoms), net, verts)
            # Add the edge to the System
            my_edges.append(my_edge)
        # Add the edges to the network's list of edges
        net.edges += my_edges

    ################################################### Create the surfaces ############################################

    # Go through the edges in the network
    for edge1 in net.edges:
        # Go through the edge's atoms combinations
        for i in range(3):
            atoms = {edge1.atoms[i], edge1.atoms[(i + 1) % 3]}
            # If the surface has been found before continue
            if check_surf(atoms, net.surfs):
                continue
            # Put together a list of edges that have our atoms
            edges = []
            # Create a variable for the edges that are doublets
            num_dubs = 0
            # Go through the edges in the system
            for edge2 in net.edges:
                # If the surface's atoms are in the edge add it
                if atoms.issubset(edge2.atoms):
                    edges.append(edge2)
                    # If the edge is a doublet, the vertex will only count once, so we need to add another to the count
                    if edge2.doublet:
                        num_dubs += 1
            # Put together a list of verts that have our atoms
            verts = []
            for vert2 in net.verts:
                # If the surface's atoms are shared with the vertex, add it to the list
                if atoms.issubset(vert2.atoms):
                    verts.append(vert2)
            # In order to be a true surface the number of edges need to be equal to the number of verts (+ doublets)
            if len(verts) + num_dubs == len(edges):
                my_surf = Surface(list(atoms), verts=verts, net=net, edges=edges)
                net.surfs.append(my_surf)


# Connect network function. Takes in a disconnected network of atoms, vertices, surfaces and edges and connects it
def connect(net):

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

    ############################################# Connect the vertices #################################################

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

    ########################################## Connect the Edges #######################################################

    # Add the surfs to the edges
    for edge in net.edges:
        # Reset the edges surf list
        edge.surfs = []
        # Go through the surfaces in the network
        for surf in net.surfs:
            # If the surfaces atoms are in the edges atoms add it to the edge
            if set(surf.atoms).issubset(edge.atoms):
                edge.surfs.append(surf)


# Create network function. Takes in a broken network of vertices and spits out a connected network
def build(net):

    # Filter the vertices
    filter_verts(net)
    # Make the surface and edge objects
    make_objects(net)
    # Connect the network
    connect(net)


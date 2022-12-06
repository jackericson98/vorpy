from System.calcs import *
from System.Network.vertex import Vertex
from System.Network.edge import Edge
from Visualize.mpl_visualize import *


# Find v0 function. Uses the atom finding functions to find a real verified site in the network
def find_v0(net, a0=None):
    # Find the middle sub_box of the set of boxes and
    mid = len(net.sub_boxes) // 2
    if a0 is None:
        a0s = []
        inc = 0
        # Keep grabbing atoms until we have enough to get the current a0 increment
        while len(a0s) < 5:
            a0s = net.get_atoms([[mid, mid, mid]], inc)
            inc += 1
        # Pull an atom from the atoms list
        a0 = a0s[0]
    a1s = []
    inc = 0
    # Get the 5 closest atoms to a0
    while len(a1s) < 5:
        a1s = net.get_atoms([a0.box], inc)
        inc += 1
    # Set up the a2s lists
    a2s = []
    # Check the a1s for verifiable
    for j in range(len(a1s)):
        # Add the circle check
        a2s.append([])
        inc = 0
        # Get the 20 closest atoms to a0 and the current a1
        if len(net.atoms) < 20:
            a2s[j] = net.atoms
        else:
            while len(a2s[j]) < 20:
                a2s[j] = net.get_atoms([[mid, mid, mid]], inc)
                inc += 1

        # Filter out the circles that don't work

        # Set up verified circles list for this a1
        verified_circles = []
        # Check each of the combinations for this a1
        for a2 in a2s[j]:
            # Use an edge object as a vehicle for calculating and verifying the inscribed circle
            edge = Edge(atoms=[a0, a1s[j], a2])
            edge.get_loc()
            # If a circle can be made and the site does not overlap with any other atoms, add it to the list
            if edge.loc is not None and edge.rad < net.max_vert and verify_site(edge, net):
                verified_circles.append(edge.atoms)

        # Test for verified sites

        # Try to make a verified v0 site with the verified circles
        for circle in verified_circles:
            # Try to create a vertex
            myVert = find_site(net, circle)
            # Check for a real site
            if myVert is not None and myVert[0].loc is not None:
                return myVert[0]


# Verify site function. Compares a vertex to the atoms around to see if they overlap
def verify_site(vert, net):
    # Grad the location and radius of the check vertex
    loc, rad = vert.loc, vert.rad
    # Find the indices of the sub-box for the vertex
    vi = int((loc[0] - net.box[0][0]) / net.sub_box_size[0])
    vj = int((loc[1] - net.box[0][1]) / net.sub_box_size[1])
    vk = int((loc[2] - net.box[0][2]) / net.sub_box_size[2])
    # Get the number of boxes that an overlapping atom could possibly be away from the vertex sub-box
    atom_range = int(rad / min(net.sub_box_size) + net.max_atom_rad / min(net.sub_box_size)) + 2
    # Get the atoms in that range
    overlap_test_atoms = net.get_atoms([[vi, vj, vk]], atom_range)
    # Set up the overlap tracker
    overlap = False
    # Go through the atoms in the overlap test atom list
    for atom2 in overlap_test_atoms:
        # If the atom is one of the vertex atoms move on
        if atom2 in vert.atoms:
            continue
        # If the distance between the vertex and the test atom is less than the sum of their radii, they overlap
        if calc_dist(atom2.loc, loc) < atom2.rad + rad:
            overlap = True
            break
    # If we make it all the way through the list of close atoms without overlapping it is a viable vertex
    if not overlap:
        return True
    return False


# Find site function. Used a vertex and a combination of it's edge atoms to find the connecting vertex
def find_site(net, edge_atoms, vn_1=None):
    # Get the atoms that should not ba a part of the new vertex
    if vn_1 is None:
        vert_atoms = edge_atoms
    else:
        vert_atoms = vn_1.atoms
    # Set up check sol variable
    check_sol = False
    # Check the edge atoms to see if they are sol atoms
    if (not net.sol_verts) and edge_atoms[0].res.lower() == 'sol' and \
            edge_atoms[1].res.lower() == 'sol' and edge_atoms[2].res.lower() == 'sol':
        check_sol = True
    # Set up a list of atoms to test our edge atoms with
    test_atoms = []
    inc = 0
    # Grab the atoms we want to test against
    while (len(test_atoms) < 10 or len(test_atoms) < len(net.atoms)) and inc < 5:
        test_atoms += net.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], inc)
        inc += 1
    # Instantiate the vertex list and the size limit for vertices found
    verts = []
    vert_ndx_list_locs = []
    # Go through each atom in the given test atoms. Extremely optimized
    for atom in test_atoms:
        # Reset the doublet variable
        doublet = None
        # If the atom is in the previous vertex move on
        if atom in vert_atoms or check_sol and atom.res.lower() == 'sol':
            continue
        # If we have found the vertex before it is not the previous vertex return
        atom_ndxs = [net.atoms.index(atom1) for atom1 in edge_atoms + [atom]]
        atom_ndxs.sort()
        # Get the vertex's index/insert index
        vert_ndx = search_verts(net.vert_ndxs, atom_ndxs)
        # If the vertex has been found before connect it to the previous one and return
        if vert_ndx < len(net.vert_ndxs) and net.vert_ndxs[vert_ndx] == atom_ndxs:
            return
        # Create the vertex and calculate its value
        vert = Vertex(edge_atoms + [atom], net=net)
        vert.calc_vert()
        # Catch the none location case
        if vert.loc is None:
            continue

        # Create the vertex's doublet if it exists
        if vert.loc2 is not None and abs(vert.rad2) < net.max_vert:
            # Create the alternate vertex for the doublet site
            doublet = Vertex(location=vert.loc2, radius=vert.rad2, atoms=vert.atoms, net=net, doublet=vert,
                             loc2=vert.loc, rad2=vert.rad, ndx=vert.ndx)

        # Filter the vertex out if it is too large or not able to be made
        if abs(vert.rad) < net.max_vert and verify_site(vert, net):
            verts.append(vert)
            vert_ndx_list_locs.append(vert_ndx)
            # If the first vertex site is a valid site add it to the list of check vertices and add its index
            if doublet is not None and verify_site(doublet, net):
                vert.doublet = doublet
        # Check to see if the doublet's site is verified
        elif doublet is not None and verify_site(doublet, net):
            doublet.doublet = None
            verts.append(doublet)
            vert_ndx_list_locs.append(vert_ndx)


    # If no verts have been found return
    if len(verts) == 0:
        return
    # If we find only 1 vertex, return it
    elif len(verts) == 1:
        return verts[0], vert_ndx_list_locs[0]
    # If there are multiple vertices find the smallest one
    elif len(verts) >= 2:
        # Instantiate the return vertex, its relative location in the vertex list and the comparison radius
        myVert = None
        min_rad = np.inf
        myVert_ndx = None
        # Go through the list of vertices
        for j in range(len(verts)):
            vert = verts[j]
            if vert.rad < min_rad:
                myVert = vert
                min_rad = vert.rad
                myVert_ndx = vert_ndx_list_locs[j]
        # Return the smallest vertex and where it belongs in the network's list of sorted vertex indices
        return myVert, myVert_ndx


# Find network function. Keeps searching the network until all verts are found
def find_verts(net, a0=None):
    # Calculate the total number of vertices
    if net.sol_verts:
        tot_verts = 6 * len(net.atoms)
    else:
        sol_len = 0
        if net.sys.sol is not None:
            sol_len = len(net.sys.sol)
        tot_verts = 16 * (len(net.atoms) - sol_len)
    # Find the first verified vertex
    if len(net.atoms) == 4:
        v0 = Vertex(net.atoms, net)
        v0.calc_vert()
    else:
        v0 = find_v0(net, a0)

    # If no v0 is possible (e.g., a lone atom) return
    if v0 is None:
        return
    # Check if this is the first go around
    if net.verts is None:
        net.verts = [v0]
        net.vert_ndxs = [v0.ndx]
        net.edges = []
    else:
        my_ndx = search_verts(net.vert_ndxs, v0.ndx)
        net.verts.insert(my_ndx, v0)
        net.vert_ndxs.insert(my_ndx, v0.ndx)
    # Set up the vertex stack
    vert_stack = [v0]
    # While the verts stack is not empty
    while vert_stack:
        # Get the vertex from the top of the stack
        vert = vert_stack.pop()
        # Set up the edge stack
        e_stack = [[[vert.atoms[i], vert.atoms[(i + 1) % 4], vert.atoms[(i + 2) % 4]], vert] for i in range(4)]
        # While the edge stack is not empty
        while e_stack:
            # Get the percentage and print it
            percentage = min((len(net.verts) / tot_verts) * 100, 100)
            print("\rFinding vertices: {:.2f} %".format(percentage), end="")
            # Get the edge from the top of the stack
            edge_atoms, vert = e_stack.pop()
            # Find the next site in the network
            myVert = find_site(net, edge_atoms, vert)
            # If the vertex is none continue
            if myVert is None:
                continue

            myVert, myVert_ndx = myVert

            # Add the vertex to the stack and the network
            vert_stack.append(myVert)
            # Insert the vertices in order of increasing atom indices
            net.verts.insert(myVert_ndx, myVert)
            net.vert_ndxs.insert(myVert_ndx, myVert.ndx)
            # Remove the atoms from the
            for atom in myVert.atoms:
                atom_ndx = net.atoms.index(atom)
                if net.atoms[atom_ndx].res.lower() != 'sol' and atom_ndx in net.atom_ndxs:
                    net.atom_ndxs.remove(atom_ndx)

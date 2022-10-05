from System.calcs import *
from System.Network.vertex import Vertex


# Find v0 function. Finds the first vertex in the network
def find_v0(net, a0=None, n=0):
    # Check for if the network does not have enough atoms
    if len(net.atoms) < 4:
        return
    # If no a0 is given
    if a0 is None:
        # Find the middle sub_box of the set of boxes and
        mid = len(net.sub_boxes) // 2
        atoms = []
        inc = 1
        while not atoms:
            atoms = net.get_atoms([[mid, mid, mid]], inc)
            inc += 1
        a0 = atoms[-1 - n]
    # Find the set of atoms with the minimum distance between surfaces
    min_dist = np.inf
    a1 = None
    inc = 0
    while a1 is None and inc <= len(net.sub_boxes) + 1:
        # Get the atoms in the system
        atoms = net.get_atoms([a0.box], inc)
        # Go through each atom determining the atom with the minimum distance between it and a0's surfaces
        for atom in atoms:
            # Skip a0
            if atom == a0:
                continue
            # Set the new atom distances
            a_dist = calc_dist(a0.loc, atom.loc) - (a0.rad + atom.rad)
            # If the new atom distance is less than the previous minimum distance update the variables
            if a_dist < min_dist:
                min_dist = a_dist
                a1 = atom
        inc += 1
    # Find the set of atoms with the minimum inscribed circle
    min_rad = np.inf
    a2 = None
    inc = 0
    while a2 is None and inc <= len(net.sub_boxes) + 1:
        # Get the atoms from the network
        atoms = net.get_atoms([a0.box, a1.box], inc + 1)
        # Go through each other atom to determine the smallest circle that can be made with our 2 atoms and a third
        for atom in atoms:
            # Skip a0, a1
            if atom == a0 or atom == a1:
                continue
            # Calculate the circle made with the 3 atoms
            circ = calc_circ([a0, a1, atom])
            # If the radius of the inscribed circle is smaller than the previous smallest found circle's rad replace
            if circ and abs(circ[1]) < min_rad:
                min_rad = abs(circ[1])
                a2 = atom
        inc += 1
    if a2 is None:
        return
    # Find the set of atoms with the minimum inscribed sphere
    myVert = find_site(net, [a0, a1, a2])
    # Keep recursively calling the find_v0 function until it finds a valid site that isn't a doublet
    if myVert is None:
        myVert = find_v0(net, a0=net.atoms[n], n=n+1)
    else:
        # If we find v0 return the vertex that was found
        myVert = myVert[0]
    # Return the vertex
    return myVert


# Verify site function. Compares a vertex to the atoms around to see if they overlap
def verify_site(vert, net, doublet_check=False):
    # Grad the location and radius of the check vertex
    loc, rad = vert.loc, vert.rad
    if doublet_check:
        loc, rad = vert.loc2, vert.rad2
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


# Find site function. Currently, overkill, searching through all atoms for overlap and
def find_site(net, edge_atoms, vn_1=None):
    # Get the atoms that should not ba a part of the new vertex
    if vn_1 is None:
        vert_atoms = edge_atoms
    else:
        vert_atoms = vn_1.atoms
    # Set up a list of atoms to test our edge atoms with
    test_atoms = []
    inc = 0
    # Grab the atoms we want to test against
    while (len(test_atoms) < 10 or len(test_atoms) < len(net.atoms)) and inc < 5:
        test_atoms += net.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], inc, exclusive=True)
        inc += 1
    # Instantiate the vertex list and the size limit for vertices found
    verts = []
    vert_ndx_list_locs = []
    # Go through each atom in the network --> This can easily be improved
    for atom in test_atoms:
        # If the atom is in the previous vertex move on
        if atom in vert_atoms:
            continue
        # If we have found the vertex before it is not the previous vertex return
        atom_ndxs = [net.atoms.index(atom1) for atom1 in edge_atoms + [atom]]
        atom_ndxs.sort()
        # Get the vertex's index/insert index
        vert_ndx = search_verts(net.vert_ndxs, atom_ndxs)
        # If the found vertex index is less than the # of vertices and the found vertex index list matches ours, return
        if vert_ndx < len(net.vert_ndxs) and net.vert_ndxs[vert_ndx] == atom_ndxs:
            return
        # Create the vertex
        vert = Vertex(edge_atoms + [atom], net=net)
        # Filter the vertex out if it is too large or not able to be made
        if vert.loc is None or vert.rad > net.beta_val:
            continue
        # For doublet cases verify differently
        if vert.loc2 is not None:
            # If the first vertex site is a valid site add it to the list of check vertices and add its index
            if verify_site(vert, net):
                verts.append(vert)
                vert_ndx_list_locs.append(vert_ndx)
                # If the second vertex's radius is less than the min_rad, and it is a verified site mark it as a doublet
                if vert.rad2 < net.beta_val and verify_site(vert, net, doublet_check=True):
                    vert.doublet = True
            # If the first vertex is not a verified site, test the second location
            elif vert.rad2 < net.beta_val and verify_site(vert, net, doublet_check=True):
                # Replace the location and radius for the vertex and add it to the list
                vert.loc, vert.rad = vert.loc2, vert.rad2
                verts.append(vert)
                vert_ndx_list_locs.append(vert_ndx)
            # If neither sites are verifiable continue
            else:
                continue
        # If the site is verified add it to the list of potential vertices
        elif verify_site(vert, net):
            verts.append(vert)
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
def find_vertices(net, v0=None, print_loc=None):
    # If print_loc is given, set the
    if print_loc is not None:
        tot_verts = print_loc[1]
        benchmark = print_loc[0]
    else:
        tot_verts = len(net.atoms) * 6
        benchmark = 0
    # If no vert is given get one
    if v0 is None:
        # Find the first vertex in the System
        v0 = find_v0(net)
    # Add v0 to the System
    net.verts.append(v0)
    # Set up the vertex stack
    vert_stack = [v0]
    # While the verts stack is not empty
    while vert_stack:
        # Running print statement giving an estimate for percentage of the network that has been created
        percentage = float(np.round(100 * ((benchmark + len(net.verts)) / tot_verts), 2))
        print("\rBuilding Network:  ",
              '#' * (int(percentage) // 10) + ' ' * (10 - (int(percentage) // 10)), percentage, "%", end='')
        # Get the vertex from the top of the stack
        vert = vert_stack.pop()
        # Set up the edge stack
        e_stack = [[[vert.atoms[i], vert.atoms[(i + 1) % 4], vert.atoms[(i + 2) % 4]], vert] for i in range(4)]
        # While the edge stack is not empty
        while e_stack:
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
            # Add the vertex to the atoms
            for atom in myVert.atoms:
                atom.verts.append(myVert)

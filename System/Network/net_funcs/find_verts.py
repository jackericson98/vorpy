from System.sys_funcs.calcs import ndx_search, calc_dist
from System.Network.net_objs.vertex import Vertex
from System.Network.net_objs.edge import Edge
from numpy import sqrt, array, square


# Find v0 function. Uses the atom finding functions to find a real verified site in the network
def find_v0(net, a0=None, group_atoms=None):
    # Check to see if we need a group atom's box
    if a0 is not None:
        my_box = a0.box
    elif group_atoms is not None:
        my_box = net.atoms[group_atoms[0]].box
    else:
        # Find the middle sub_box of the set of boxes and
        mid = len(net.sub_boxes) // 2
        my_box = [mid, mid, mid]
    if a0 is None:
        a0s = []
        inc = 0
        # Keep grabbing atoms until we have enough to get the current a0 increment
        while len(a0s) < 5:
            a0s = net.get_atoms([my_box], inc)
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
    a2s, j = [], 0
    # Check the a1s for verifiable
    while len(a1s) > 0:
        # Get the a1
        a1 = a1s.pop()
        # Add the circle check
        a2s.append([])
        inc = 0
        # Get the 20 closest atoms to a0 and the current a1
        if len(net.atoms) < 20:
            a2s[j] = net.atoms.copy()
        else:
            while len(a2s[j]) < 20:
                a2s[j] = net.get_atoms([my_box], inc)
                inc += 1
        # Set up verified circles list for this a1
        verified_circles = []
        # Check each of the combinations for this a1
        for a2 in a2s[j]:
            # Use an edge object as a vehicle for calculating and verifying the inscribed circle
            edge = Edge(atoms=[a0, a1, a2], net=net)
            edge.get_loc()
            # If a circle can be made and the site does not overlap with any other atoms, add it to the list
            if edge.loc is not None and edge.rad < net.max_vert and verify_site(edge, net):
                verified_circles.append(edge.atoms)
        # Try to make a verified v0 site with the verified circles
        for circle in verified_circles:
            # Try to create a vertex
            myVert = find_site(net, circle, group_atoms=group_atoms)
            # Check for a real site
            if myVert is not None and myVert[0].loc is not None:
                return myVert[0]
        j += 1


# Verify site function. Compares a vertex to the atoms around to see if they overlap
def verify_site(vert, net):
    # Grad the location and radius of the check vertex
    loc, rad = vert.loc, vert.rad
    # Find the indices of the sub-box for the vertex
    vi, vj, vk = [int((loc[i] - net.box[0][i]) / net.sub_box_size[i]) for i in range(3)]
    # Check to see if the sub box even exists
    if vi > net.box_max[0] or vj > net.box_max[1] or vk > net.box_max[2] or vi < 0 or vj < 0 or vk < 0:
        return False
    # Checked atoms list
    checked_atoms = [_ for _ in vert.ndx]
    # Quick check to see if any atoms exist inside the vertex's box
    quick_atoms = net.sub_boxes[vi][vj][vk]
    for atom in quick_atoms:
        # If the atom is one of the vertex atoms move on
        if atom.num in checked_atoms:
            continue
        my_radius = atom.rad
        if net.type == 'del':
            if sqrt(sum(square(array(atom.loc) - array(loc)))) < rad:
                return False
        # I don't know how to verify power yet
        elif net.type == 'pow':
            if sqrt(sum(square(array(atom.loc) - array(loc)))) ** 2 - my_radius ** 2 < rad:
                return False
        # Verification for a voronoi network
        elif net.type == 'vor':
            if sqrt(sum(square(array(atom.loc) - array(loc)))) < my_radius + rad:
                return False
        # Add the atom to the checked atoms list
        checked_atoms.append(atom.num)
    inc = 1
    # Get the number of boxes that an overlapping atom could possibly be away from the vertex sub-box
    min_sub_box_size = min(net.sub_box_size)
    atom_range = int(rad / min_sub_box_size + net.sys.max_atom_rad / min_sub_box_size) + 2
    while inc <= atom_range:
        # Get the atoms in that range
        overlap_test_atoms = net.get_atoms([[vi, vj, vk]], inc)
        # Go through the atoms in the overlap test atom list
        for atom in overlap_test_atoms:
            # If the atom is one of the vertex atoms move on
            if atom.num in checked_atoms:
                continue
            my_radius = atom.rad
            if net.type == 'del':
                if sqrt(sum(square(array(atom.loc) - array(loc)))) < rad:
                    return False
            # I don't know how to verify power yet
            elif net.type == 'pow':
                if sqrt(sum(square(array(atom.loc) - array(loc)))) ** 2 - my_radius ** 2 < rad:
                    return False
            # Verification for a voronoi network
            elif net.type == 'vor':
                if sqrt(sum(square(array(atom.loc) - array(loc)))) < my_radius + rad:
                    return False
            # If the distance between the vertex and the test atom is less than the sum of their radii, they overlap

            checked_atoms.append(atom.num)
        inc += 1
    # If we make it all the way through the list of close atoms without overlapping it is a viable vertex
    return True


# Find site function. Used a vertex and a combination of it's edge atoms to find the connecting vertex
def find_site(net, edge_atoms, vn_1=None, first=False, group_atoms=None):
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = [_.num for _ in edge_atoms]
    # Check if the edge contains a group atom or not
    check_atoms = True
    for ndx in edge_ndxs:
        if group_atoms is not None and ndx in group_atoms:
            check_atoms = False
            break
    # If the previous vertex has been provided, add the other atom to the not allowed atoms
    if vn_1 is None:
        vert_atom_ndxs = edge_ndxs
    else:
        vert_atom_ndxs = vn_1.ndx
    # Set up a list of atoms to test our edge atoms with
    test_atoms = []
    inc = 0
    max_inc = int(net.max_vert / min(net.sub_box_size) - net.sys.max_atom_rad) + 1
    if first:
        max_inc = 5
    # Grab the atoms we want to test against
    while (len(test_atoms) < 10 or len(test_atoms) < len(net.atoms)) and inc < max_inc:
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
        if atom.num in vert_atom_ndxs:
            continue
        # Check if we need to check and if so check for the atom in the list
        if check_atoms and atom.num not in group_atoms:
            continue
        # If we have found the vertex before it is not the previous vertex return
        atom_ndxs = edge_ndxs + [atom.num]
        atom_ndxs.sort()
        # Get the vertex's index/insert index
        vert_ndx = ndx_search(net.vert_ndxs, atom_ndxs)
        # If the vertex has been found before connect it to the previous one and return
        if vert_ndx < len(net.vert_ndxs) and net.vert_ndxs[vert_ndx] == atom_ndxs:
            return
        # Create the vertex and calculate its value
        vert = Vertex(edge_atoms + [atom], net=net)
        # Calculate the correct vertex values
        if net.type == 'pow':
            vert.calc_flat_vert(power=True)
        elif net.type == 'del':
            vert.calc_flat_vert(power=False)
        else:
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
            if len(verts) > 0 and verts[0].rad < vert.rad:
                return verts[0], vert_ndx_list_locs[0]
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
    elif len(verts) == 1 or verts[0].rad < verts[1].rad:
        return verts[0], vert_ndx_list_locs[0]
    return verts[1], vert_ndx_list_locs[1]


# Find network function. Keeps searching the network until all verts are found
def find_verts(net, a0=None, my_group=None):
    # Get the group atoms from which to check vertices against
    if my_group is None or (my_group is not None and len(my_group.atoms) == len(net.atoms)):
        group_atoms = [i for i in range(len(net.atoms))]
        # Calculate the rough number of vertices
        tot_verts = 7 * len(net.atoms)
    # If a group was provided make sure to get its indices
    elif my_group is not None:
        group_atoms = my_group.atom_ndxs
        # Calculate the number of vertices
        tot_verts = 7 * len(group_atoms) + int(60 * sqrt(len(group_atoms)))
    else:
        return
    # Find the first verified vertex
    if len(group_atoms) == 4:
        v0 = Vertex(net.atoms, net)
        v0.calc_vert()
    else:
        v0 = find_v0(net, a0, group_atoms)
    # If no v0 is possible (e.g., a lone atom) return
    if v0 is None:
        return
    # Check if this is the first go around
    if net.verts is None:
        net.verts = [v0]
        net.vert_ndxs = [v0.ndx]
        net.edges = []
    else:
        my_ndx = ndx_search(net.vert_ndxs, v0.ndx)
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
            print("\rfinding vertices: {:.2f} %".format(percentage), end="")
            # Get the edge from the top of the stack
            edge_atoms, vert = e_stack.pop()
            # Find the next site in the network
            myVert = find_site(net=net, edge_atoms=edge_atoms, vn_1=vert, group_atoms=group_atoms)
            # If the vertex is none continue
            if myVert is None:
                continue
            # Set the vertex and its index
            myVert, myVert_ndx = myVert
            # Add the vertex to the stack and the network
            vert_stack.append(myVert)
            # Insert the vertices in order of increasing atom indices
            net.verts.insert(myVert_ndx, myVert)
            net.vert_ndxs.insert(myVert_ndx, myVert.ndx)
            # Remove the atoms from the
            for atom in myVert.atoms:
                if atom.num in net.atom_ndxs:
                    net.atom_ndxs.remove(atom.num)

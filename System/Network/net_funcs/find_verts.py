from System.sys_funcs.calcs import ndx_search
from System.Network.net_objs.vertex import Vertex
from System.Network.net_objs.edge import Edge
import numpy as np


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
                a2s[j] = net.get_atoms([[mid, mid, mid]], inc)
                inc += 1
        # Filter out the circles that don't work

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
            myVert = find_site(net, circle)
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
    # Set up the numpy checks for faster reference
    sqrt, array, square = np.sqrt, np.array, np.square
    # Quick check to see if any atoms exist inside the vertex's box
    quick_atoms = net.sub_boxes[vi][vj][vk]
    for atom in quick_atoms:
        # If the atom is one of the vertex atoms move on
        if atom.num in checked_atoms:
            continue
        # If the distance between the vertex and the test atom is less than the sum of their radii, they overlap
        if sqrt(sum(square(array(atom.loc) - array(loc)))) < atom.rad + rad:
            return False
        # Add the atom to the checked atoms list
        checked_atoms.append(atom.num)
    inc = 1
    # Get the number of boxes that an overlapping atom could possibly be away from the vertex sub-box
    min_sub_box_size = min(net.sub_box_size)
    atom_range = int(rad / min_sub_box_size + net.max_atom_rad / min_sub_box_size) + 2
    while inc <= atom_range:
        # Get the atoms in that range
        overlap_test_atoms = net.get_atoms([[vi, vj, vk]], inc)
        # Go through the atoms in the overlap test atom list
        for atom in overlap_test_atoms:
            # If the atom is one of the vertex atoms move on
            if atom.num in checked_atoms:
                continue
            # If the distance between the vertex and the test atom is less than the sum of their radii, they overlap
            if sqrt(sum(square(array(atom.loc) - array(loc)))) < atom.rad + rad:
                return False
            checked_atoms.append(atom.num)
        inc += 1
    # If we make it all the way through the list of close atoms without overlapping it is a viable vertex
    return True


# Find site function. Used a vertex and a combination of it's edge atoms to find the connecting vertex
def find_site(net, edge_atoms, vn_1=None, first=False):
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = [_.num for _ in edge_atoms]
    if vn_1 is None:
        vert_atom_ndxs = edge_ndxs
    else:
        vert_atom_ndxs = vn_1.ndx
    # Set up a list of atoms to test our edge atoms with
    test_atoms = []
    inc = 0
    max_inc = int(net.max_vert / min(net.sub_box_size) - net.max_atom_rad) + 1
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
def find_verts(net, a0=None, group=None):
    # Calculate the total number of vertices
    tot_verts = 7 * len(net.atoms)
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
            myVert = find_site(net, edge_atoms, vert)
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

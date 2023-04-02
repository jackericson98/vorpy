from System.sys_funcs.calcs import ndx_search, get_time
from System.Network.net_objs.vertex import Vertex
from System.Network.net_objs.edge import Edge
from System.Network.net_funcs.build_edge import get_edge_loc
from numpy import sqrt, array, square
import numpy as np
import time


# Calculate vertex function. Takes in 4 atoms, calculates the loc and rad of the inscribed sphere and adds the
def calc_vert(locs, rads):
    # The real location and radius of the base sphere
    if type(locs[0]) == list:
        locs = [np.array(_) for _ in locs]
    r0, r1, r2, r3 = rads
    r0_2 = r0 ** 2
    # Find the recalculated location of the atoms
    l0, l1, l2, l3 = locs[0], locs[1] - locs[0], locs[2] - locs[0], locs[3] - locs[0]
    # Calculate our System of linear equations coefficients
    a1, b1, c1, d1, f1 = 2 * l1[0], 2 * l1[1], 2 * l1[2], 2 * (r1 - r0), r0_2 - r1 ** 2 + l1[0] ** 2 + l1[
        1] ** 2 + l1[2] ** 2
    a2, b2, c2, d2, f2 = 2 * l2[0], 2 * l2[1], 2 * l2[2], 2 * (r2 - r0), r0_2 - r2 ** 2 + l2[0] ** 2 + l2[
        1] ** 2 + l2[2] ** 2
    a3, b3, c3, d3, f3 = 2 * l3[0], 2 * l3[1], 2 * l3[2], 2 * (r3 - r0), r0_2 - r3 ** 2 + l3[0] ** 2 + l3[
        1] ** 2 + l3[2] ** 2
    # Calculate the F values
    F = a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1
    F_2 = F ** 2
    F10 = b1 * c2 * f3 - b1 * c3 * f2 - b2 * c1 * f3 + b2 * c3 * f1 + b3 * c1 * f2 - b3 * c2 * f1
    F11 = -b1 * c2 * d3 + b1 * c3 * d2 + b2 * c1 * d3 - b2 * c3 * d1 - b3 * c1 * d2 + b3 * c2 * d1
    F20 = -a1 * c2 * f3 + a1 * c3 * f2 + a2 * c1 * f3 - a2 * c3 * f1 - a3 * c1 * f2 + a3 * c2 * f1
    F21 = a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1
    F30 = a1 * b2 * f3 - a1 * b3 * f2 - a2 * b1 * f3 + a2 * b3 * f1 + a3 * b1 * f2 - a3 * b2 * f1
    F31 = -a1 * b2 * d3 + a1 * b3 * d2 + a2 * b1 * d3 - a2 * b3 * d1 - a3 * b1 * d2 + a3 * b2 * d1
    # Calculate the ranks of the matrices
    m_rank, f_rank = 3, 3
    if F == 0:
        my_mtx = [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f3]]
        m_rank = np.linalg.matrix_rank(np.array(my_mtx[:-1]))
        if m_rank != 3:
            f_rank = np.linalg.matrix_rank(np.array(my_mtx))
    verts = []
    xs, ys, zs, Rs = [], [], [], []
    # Case 1:
    if F != 0:
        # Calculate the radius polynomial coefficients
        a = ((F11 ** 2 + F21 ** 2 + F31 ** 2) / F_2) - 1
        b = 2 * (((F10 * F11 + F20 * F21 + F30 * F31) / F_2) - r0)
        c = ((F10 ** 2 + F20 ** 2 + F30 ** 2) / F_2) - r0 ** 2
        # If the discriminant is positive, find the real positive roots of the quadratic
        if -4 * a * c + b ** 2 >= 0:
            Rs = [R for R in np.roots([a, b, c]) if np.isreal(R)]
        # Instantiate the verts array
        verts = []
        # Go through each radius and calculate the vertex
        for R in Rs:
            x, y, z = F10 / F + R * F11 / F, F20 / F + R * F21 / F, F30 / F + R * F31 / F
            # Move the vertex back to the actual location of the atoms
            verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
    # Case 2:
    elif a1 * b2 - a2 * b1 != 0 and m_rank == 3 and f_rank == 3 and F > 0:
        # Calculate the _ polynomial coefficients
        a = F_2 + F11 ** 2 + F21 ** 2 - F31 ** 2
        b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * r0)
        c = F10 ** 2 + F20 ** 2 - (F30 + F * r0)
        roots, verts = [], []
        # Check the discriminant
        disc = -4 * a * c + b ** 2
        if disc > 0:
            roots = [root for root in np.roots([a, b, c]) if np.isreal(root)]
        # Case 2 subcases:
        # Case 2.1
        if F31 != 0:
            # Go through each radius and calculate the vertex
            for z in roots:
                x, y, R = F10 / F + z * F11 / F, F20 / F + z * F21 / F, F30 / F + z * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
        # Case 2.2
        elif F21 != 0:
            # Go through each radius and calculate the vertex
            for y in roots:
                x, R, z = F10 / F + y * F11 / F, F20 / F + y * F21 / F, F30 / F + y * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
        # Case 2.3
        elif F11 != 0:
            # Go through each radius and calculate the vertex
            for x in roots:
                R, y, z = F10 / F + x * F11 / F, F20 / F + x * F21 / F, F30 / F + x * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append([[x + l0[0], y + l0[1], z + l0[2]], R])
    loc, rad, loc2, rad2 = None, None, None, None
    # If one root exists return it
    if len(verts) == 1:
        loc, rad = verts[0][0], verts[0][1]
    # If two roots exist:
    elif len(verts) == 2:
        # Get the largest atom's radius
        max_atom_rad = max([r0, r1, r2, r3])
        # Set the locations and radii, so that the smaller vertex is first
        if abs(verts[0][1]) > abs(verts[1][1]):
            verts[0], verts[1] = verts[1], verts[0]
        # Set the locations and radii variables
        locs, rads = [verts[0][0], verts[1][0]], [verts[0][1], verts[1][1]]
        # If either radii are negative (I'm not sure if this is possible, but let's catch it anyway)
        if rads[0] < 0 or rads[1] < 0:
            if rads[0] > 0 or abs(rads[0]) < max_atom_rad:
                loc, rad = locs[0], rads[0]
                if rads[1] > 0 or abs(rads[1]) < max_atom_rad:
                    loc2, rad2 = locs[1], rads[1]
            elif rads[1] > 0 or abs(rads[1]) < max_atom_rad:
                loc, rad = locs[1], rads[1]
        # If both radii are positive we have a doublet. Choose the smaller vertex to be the lead vertex and set loc2
        else:
            loc, loc2, rad, rad2 = locs[0], locs[1], rads[0], rads[1]
    return loc, rad, loc2, rad2


def calc_flat_vert(locs, rads, power=False):
    """
    Calculates the flat vertex between 4 atoms by finding the intersection of the mid-point planes between the first
    atom and the others
    :param rads:
    :param locs:
    :param power:
    :return:
    """
    atom_rads = [(x, _) for _, x in sorted(zip(rads, locs), key=lambda pair: pair[0])]
    # Get the plane equations
    coeffs = []
    # Go through the atoms to make the planes
    for an in atom_rads[1:]:
        # Get the point between the atoms
        r = np.array(an[0]) - np.array(atom_rads[0][0])
        norm = np.linalg.norm(r)
        rn = r / norm
        if power:
            d0 = 0.5 * (norm ** 2 + atom_rads[0][1] ** 2 - an[1] ** 2) / norm
            center = atom_rads[0][0] + d0 * rn
        else:
            center = 0.5 * r + np.array(atom_rads[0][0])
        coeffs.append(rn.tolist() + [np.dot(rn, center)])

    x1, y1, z1, c1 = coeffs[0]
    x2, y2, z2, c2 = coeffs[1]
    x3, y3, z3, c3 = coeffs[2]

    disc = z1 * y2 * x3 - y1 * z2 * x3 - z1 * x2 * y3 + x1 * z2 * y3 + y1 * x2 * z3 - x1 * y2 * z3
    x_numerator = c1 * z2 * y3 - z1 * c2 * y3 - c1 * y2 * z3 + y1 * c2 * z3 + z1 * y2 * c3 - y1 * z2 * c3
    y_numerator = - c1 * z2 * x3 + z1 * c2 * x3 + c1 * x2 * z3 - x1 * c2 * z3 - z1 * x2 * c3 + x1 * z2 * c3
    z_numerator = c1 * y2 * x3 - y1 * c2 * x3 - c1 * x2 * y3 + x1 * c2 * y3 + y1 * x2 * c3 - x1 * y2 * c3
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
    # Get the radius
    if power:
        rad = np.sqrt(sum(np.square(np.array([x, y, z]) - np.array(atom_rads[0][0])))) ** 2 - atom_rads[0][1] ** 2
    else:
        rad = np.sqrt(sum(np.square(np.array([x, y, z]) - np.array(atom_rads[0][0]))))
    return [x, y, z], rad


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
            get_edge_loc(edge)
            # If a circle can be made and the site does not overlap with any other atoms, add it to the list
            if edge.loc is not None and edge.rad < net.max_vert and verify_site(edge.loc, edge.rad, edge.ndx, net):
                verified_circles.append(edge.atoms)
        # Try to make a verified v0 site with the verified circles
        for circle in verified_circles:
            # Try to create a vertex
            my_vert = find_site(net, circle, group_atoms=group_atoms)
            # Check for a real site
            if my_vert is not None and my_vert[0].loc is not None:
                return my_vert[0]
        j += 1


# Verify site function. Compares a vertex to the atoms around to see if they overlap
def verify_site(loc, rad, ndx, net):
    # Find the indices of the sub-box for the vertex
    vi, vj, vk = [int((loc[i] - net.box[0][i]) / net.sub_box_size[i]) for i in range(3)]
    # Check to see if the sub box even exists
    if vi > net.box_max[0] or vj > net.box_max[1] or vk > net.box_max[2] or vi < 0 or vj < 0 or vk < 0:
        return False
    # Checked atoms list
    checked_atoms = [_ for _ in ndx]
    # Quick check to see if any atoms exist inside the vertex's box
    quick_atoms = net.sub_boxes[vi][vj][vk]
    for atom in quick_atoms:
        # If the atom is one of the vertex atoms move on
        if atom.num in checked_atoms:
            continue
        arad = atom.rad
        if net.type == 'del':
            if sqrt(sum(square(array(atom.loc) - array(loc)))) < rad:
                return False
        # Verify power
        elif net.type == 'pow':
            if sqrt(sum(square(array(atom.loc) - array(loc)))) ** 2 - arad ** 2 < rad:
                return False
        # Verification for a voronoi network
        elif net.type == 'vor':
            if sqrt(sum(square(array(atom.loc) - array(loc)))) < arad + rad:
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
            arad = atom.rad
            if net.type == 'del':
                if sqrt(sum(square(array(atom.loc) - array(loc)))) < rad:
                    return False
            # I don't know how to verify power yet
            elif net.type == 'pow':
                if sqrt(sum(square(array(atom.loc) - array(loc)))) ** 2 - arad ** 2 < rad:
                    return False
            # Verification for a voronoi network
            elif net.type == 'vor':
                if sqrt(sum(square(array(atom.loc) - array(loc)))) - arad < rad:
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
            vert.loc, vert.rad = calc_flat_vert(locs=[_.loc for _ in vert.atoms], rads=[_.rad for _ in vert.atoms], power=True)
        elif net.type == 'del':
            vert.loc, vert.rad = calc_flat_vert(locs=[_.loc for _ in vert.atoms], rads=[_.rad for _ in vert.atoms], power=False)
        else:
            vert.loc, vert.rad, vert.loc2, vert.rad2 = calc_vert(locs=[_.loc for _ in vert.atoms], rads=[_.rad for _ in vert.atoms])
        # Catch the none location case
        if vert.loc is None:
            continue
        # Create the vertex's doublet if it exists
        if vert.loc2 is not None and abs(vert.rad2) < net.max_vert:
            # Create the alternate vertex for the doublet site
            doublet = Vertex(location=vert.loc2, radius=vert.rad2, atoms=vert.atoms, net=net, doublet=vert,
                             loc2=vert.loc, rad2=vert.rad, ndx=vert.ndx)
        # Filter the vertex out if it is too large or not able to be made
        if abs(vert.rad) < net.max_vert and verify_site(vert.loc, vert.rad, vert.ndx, net):
            if len(verts) > 0 and verts[0].rad < vert.rad:
                return verts[0], vert_ndx_list_locs[0]
            verts.append(vert)
            vert_ndx_list_locs.append(vert_ndx)
            # If the first vertex site is a valid site add it to the list of check vertices and add its index
            if doublet is not None and verify_site(doublet.loc, doublet.rad, doublet.ndx, net):
                vert.doublet = doublet
        # Check to see if the doublet's site is verified
        elif doublet is not None and verify_site(doublet.loc, doublet.rad, doublet.ndx, net):
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
        v0.loc, v0.rad, v0.loc2, v0.rad2 = calc_vert(locs=[_.loc for _ in v0.atoms], rads=[_.rad for _ in v0.atoms])
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
            my_time = time.perf_counter() - net.my_time
            h, m, s = get_time(my_time)
            print("\rRun Time = {}:{}:{:.2f} - Process: finding vertices: {:.2f} %".format(int(h), int(m), round(s, 2), percentage), end="")
            # Get the edge from the top of the stack
            edge_atoms, vert = e_stack.pop()
            # Find the next site in the network
            vert_ndx_pr = find_site(net=net, edge_atoms=edge_atoms, vn_1=vert, group_atoms=group_atoms)
            # If the vertex is none continue
            if vert_ndx_pr is None:
                continue
            # Set the vertex and its index
            my_vert, my_vert_ndx = vert_ndx_pr
            # Add the vertex to the stack and the network
            vert_stack.append(my_vert)
            # Insert the vertices in order of increasing atom indices
            net.verts.insert(my_vert_ndx, my_vert)
            net.vert_ndxs.insert(my_vert_ndx, my_vert.ndx)
            # Remove the atoms from the
            for atom in my_vert.atoms:
                if atom.num in net.atom_ndxs:
                    net.atom_ndxs.remove(atom.num)

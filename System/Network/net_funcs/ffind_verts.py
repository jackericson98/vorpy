from System.sys_funcs.calcs import calc_circ, np, calc_dist, ndx_search, calc_angle
from System.Network.net_objs.vertex import Vertex
from System.Network.net_objs.edge import Edge
from itertools import combinations

"""
Throw an f in front of each function to differentiate when importing. Functions needed:
1. ffind_v0 - Finds the start vertex by closeness
2. ffind_verts - Stack of verts, each with a sub_stack of edges. Goes from site to edge to site till net is complete
3. ffind_site - Takes in an edge and a vert and finds the closest atom
4.
"""


def calc_flat_vert_pow(atoms):
    """
    Calculates the flat vertex between 4 atoms by finding the intersection of the mid-point planes between the first
    atom and the others
    :param atoms:
    :return:
    """
    rads = [_.rad for _ in atoms]
    atom_rads = [x for _, x in sorted(zip(rads, atoms), key=lambda pair: pair[0])]
    # Get the plane equations
    coeffs = []
    # Go through the atoms to make the planes
    for an in atom_rads[1:]:
        # Get the point between the atoms
        r = np.array(an.loc) - np.array(atoms[0].loc)
        norm = np.linalg.norm(r)
        rn = r / norm
        atom_surf_dist = norm - an.rad - atoms[0].rad
        center = (0.5 * atom_surf_dist + atoms[0].rad) * rn + np.array(atoms[0].loc)
        coeffs.append(rn.tolist() + [np.dot(rn, center)])

    a, b, c, d = coeffs[0]
    e, f, g, h = coeffs[1]
    i, j, k, m = coeffs[2]

    disc = c * f * i - b * g * i - c * e * j + a * g * j + b * e * k - a * f * k
    x_numerator = d * g * j - c * h * j - d * f * k + b * h * k + c * f * m - b * g * m
    y_numerator = - d * g * i + c * h * i + d * e * k - a * h * k - c * e * m + a * g * m
    z_numerator = d * f * i - b * h * i - d * e * j + a * h * j + b * e * m - a * f * m
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
    # Get the radius
    rad = calc_dist([x, y, z], atom_rads[0].loc) - atom_rads[0].rad
    return [x, y, z], rad


def calc_flat_vert(atoms):
    """
    Calculates the flat vertex between 4 atoms by finding the intersection of the mid-point planes between the first
    atom and the others
    :param atoms:
    :return:
    """
    rads = [_.rad for _ in atoms]
    atom_rads = [x for _, x in sorted(zip(rads, atoms), key=lambda pair: pair[0])]
    # Get the plane equations
    coeffs = []
    # Go through the atoms to make the planes
    for an in atom_rads[1:]:
        # Get the point between the atoms
        r = np.array(an.loc) - np.array(atom_rads[0].loc)
        norm = np.linalg.norm(r)
        rn = r / norm
        center = 0.5 * r + np.array(atom_rads[0].loc)
        coeffs.append(rn.tolist() + [np.dot(rn, center)])

    a, b, c, d = coeffs[0]
    e, f, g, h = coeffs[1]
    i, j, k, m = coeffs[2]

    disc = c*f*i - b*g*i - c*e*j + a*g*j + b*e*k - a*f*k
    x_numerator = d*g*j - c*h*j - d*f*k + b*h*k + c*f*m - b*g*m
    y_numerator = - d*g*i + c*h*i + d*e*k - a*h*k - c*e*m + a*g*m
    z_numerator = d*f*i - b*h*i - d*e*j + a*h*j + b*e*m - a*f*m
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
    # Get the radius
    rad = calc_dist([x, y, z], atom_rads[0].loc)
    return [x, y, z], rad


def ffind_v0(net, a0=None):
    # Find the middle sub_box of the set of boxes and
    mid = len(net.sub_boxes) // 2
    mid_atoms = []
    inc = 0
    while len(mid_atoms) < len(net.atoms) and len(mid_atoms) < 20:
        mid_atoms = net.get_atoms([[mid, mid, mid]], inc)
        inc += 1
    if a0 is None:
        a0 = mid_atoms.pop(0)
    min_dist = np.inf
    a1 = None
    for atom in mid_atoms:
        my_dist = calc_dist(atom.loc, a0.loc)
        if my_dist < min_dist:
            a1 = atom
            min_dist = my_dist
    min_rad = np.inf
    a2 = None
    for atom in mid_atoms:
        if atom.num == a1.num:
            continue
        my_circ_rad = calc_circ([a0, a1, atom])
        if my_circ_rad is not None and abs(my_circ_rad[1]) < min_rad:
            a2 = atom
            min_rad = abs(my_circ_rad[1])
    v0_site = ffind_site(net=net, edge_atoms=[a0, a1, a2])
    return v0_site[0]


def verify_site_power(vert, net):
    # Grad the location and radius of the check vertex
    loc, rad = vert.loc, vert.rad
    # Find the indices of the sub-box for the vertex
    vi = int((vert.loc[0] - net.box[0][0]) / net.sub_box_size[0])
    vj = int((vert.loc[1] - net.box[0][1]) / net.sub_box_size[1])
    vk = int((vert.loc[2] - net.box[0][2]) / net.sub_box_size[2])
    # Get the number of boxes that an overlapping atom could possibly be away from the vertex sub-box
    atom_range = int(rad / min(net.sub_box_size) + net.max_atom_rad / min(net.sub_box_size)) + 2
    # Get the atoms in that range
    overlap_test_atoms = net.get_atoms([[vi, vj, vk]], atom_range)
    # Go through the atoms in the overlap test atom list
    for atom2 in overlap_test_atoms:
        # If the atom is one of the vertex atoms move on
        if atom2.num in vert.ndx:
            continue
        # Get the power distance between the atom2 and the vertex
        pow_dist = np.dot(atom2.loc, loc) - atom2.rad ** 2
        # Check that the power distance between the vert and the atom is not less than the radius of the vertex
        if pow_dist < rad:
            return False
    return True


def verify_site(vert, net):
    # Grad the location and radius of the check vertex
    loc, rad = vert.loc, vert.rad
    # Find the indices of the sub-box for the vertex
    vi = int((vert.loc[0] - net.box[0][0]) / net.sub_box_size[0])
    vj = int((vert.loc[1] - net.box[0][1]) / net.sub_box_size[1])
    vk = int((vert.loc[2] - net.box[0][2]) / net.sub_box_size[2])
    # Get the number of boxes that an overlapping atom could possibly be away from the vertex sub-box
    atom_range = int(rad / min(net.sub_box_size) + net.max_atom_rad / min(net.sub_box_size)) + 2
    # Get the atoms in that range
    overlap_test_atoms = net.get_atoms([[vi, vj, vk]], atom_range)
    # Go through the atoms in the overlap test atom list
    for atom2 in overlap_test_atoms:
        # If the atom is one of the vertex atoms move on
        if atom2.num in vert.ndx:
            continue
        # If the distance between the vertex and the test atom is less than the sum of their radii, they overlap
        if np.sqrt(sum(np.square(np.array(atom2.loc) - np.array(loc)))) < rad:
            return False
    # If we make it all the way through the list of close atoms without overlapping it is a viable vertex
    return True


def ffind_site(edge_atoms, net, vn_1=None):
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = [_.num for _ in edge_atoms]
    if vn_1 is None:
        vert_atom_ndxs = edge_ndxs
    else:
        vert_atom_ndxs = vn_1.ndx
    # Get the closest atoms around
    test_atoms = net.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], 5)
    new_test_atoms = []
    new_vert_ndxs = []
    for an in test_atoms:
        # If the atom is in the previous vertex move on
        if an.num in vert_atom_ndxs:
            continue
        new_test_atoms.append(an)
        # If we have found the vertex before it is not the previous vertex return
        atom_ndxs = edge_ndxs + [an.num]
        atom_ndxs.sort()
        # Get the vertex's index/insert index
        vert_ndx = ndx_search(net.vert_ndxs, atom_ndxs)
        # If the vertex has been found before connect it to the previous one and return
        if vert_ndx < len(net.vert_ndxs) and net.vert_ndxs[vert_ndx] == atom_ndxs:
            return
        new_vert_ndxs.append(vert_ndx)
    # Get the planes made by the test atoms and one of the edge atoms
    for i in range(len(new_test_atoms)):
        an, vert_ndx = new_test_atoms[i], new_vert_ndxs[i]
        # Get the intersecting point of the edge and the plane made by a0 and an
        vert_loc, vert_rad = calc_flat_vert(edge_atoms + [an])
        # Check that the vertex is inside the network's box
        outside = False
        for j in range(3):
            if not (net.box[0][j] <= vert_loc[j] <= net.box[1][j]):
                outside = True
        if outside:
            continue
        # Find the distance between the new vert and the old ver
        my_vert = Vertex(atoms=edge_atoms + [an], net=net, location=vert_loc, radius=vert_rad, flat_faces=True)
        # Verify that this site is real
        if verify_site(my_vert, net):
            return my_vert, vert_ndx


# Find network function. Keeps searching the network until all verts are found
def ffind_verts(net, a0=None, group=None):
    # Calculate the total number of vertices
    tot_verts = 5 * len(net.atoms)
    # Find the first verified vertex
    if len(net.atoms) == 4:
        loc, rad = calc_flat_vert(net.atoms)
        v0 = Vertex(atoms=net.atoms, net=net, location=loc, radius=rad, flat_faces=True)
        v0.calc_vert()
    else:
        v0 = ffind_v0(net, a0)
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
            myVert = ffind_site(edge_atoms=edge_atoms, net=net, vn_1=vert)
            # If the vertex is none continue
            if myVert is None or (len(myVert) > 1 and myVert[1] is None):
                continue
            # Set the vertex and its index
            myVert, myVert_ndx = myVert
            # Add the vertex to the stack and the network
            vert_stack.append(myVert)
            # Insert the vertices in order of increasing atom indices
            net.verts.insert(myVert_ndx, myVert)
            net.vert_ndxs.insert(myVert_ndx, myVert.ndx)
            # Add the minimum vertex size
            # Remove the atoms from the
            for atom in myVert.atoms:
                atom_ndx = net.atoms.index(atom)
                if atom_ndx in net.atom_ndxs:
                    net.atom_ndxs.remove(atom_ndx)

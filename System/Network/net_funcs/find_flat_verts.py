from System.sys_funcs.calcs import calc_circ, np, calc_dist, ndx_search
from System.Network.net_objs.vertex import Vertex
from System.Network.net_objs.edge import Edge

"""
Throw an f in front of each function to differentiate when importing. Functions needed:
1. ffind_v0 - Finds the start vertex by closeness
2. ffind_verts - Stack of verts, each with a sub_stack of edges. Goes from site to edge to site till net is complete
3. ffind_site - Takes in an edge and a vert and finds the closest atom
4.
"""


def calc_flat_vert(atoms):
    # Get the plane equations
    coeffs = []
    # Go through the atoms to make the planes
    for an in atoms[1:]:
        # Get the point between the atoms
        r = np.array(an.loc) - np.array(atoms[0].loc)
        rn = r / np.linalg.norm(r)
        center = 0.5 * r + np.array(atoms[0].loc)
        coeffs.append(rn.tolist() + [np.dot(rn, center)])

    a, b, c, d = coeffs[0]
    e, f, g, h = coeffs[1]
    i, j, k, m = coeffs[2]

    disc = c*f*i - b*g*i - c*e*j + a*g*j + b*e*k - a*f*k
    x_numerator = d*g*j - c*h*j - d*f*k + b*h*k + c*f*m - b*g*m
    y_numerator = - d*g*i + c*h*i + d*e*k - a*h*k - c*e*m + a*g*m
    z_numerator = d*f*i - b*h*i - d*e*j + a*h*j + b*e*m - a*f*m
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
    return [x, y, z]


def ffind_v0(net, a0=None):
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
        if a1 == a0:
            continue
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
        # Set up verified circles list for this a1
        min_rad, my_circ, my_a2 = np.inf, None, None
        # Check each of the combinations for this a1
        for a2 in a2s[j]:
            # Calculate the circle between the a2 atoms
            test_circ = calc_circ([a0, a1, a2])
            # Get the smallest non-None circle
            if test_circ is not None and test_circ[1] < min_rad:
                my_circ, min_rad, my_a2 = test_circ, test_circ[1], a2
        # Get the nearest atoms
        test_atoms = net.get_atoms([a0.box, a1.box, my_a2.box], 5)
        # Set up the tracker variables
        min_dist, my_vert, my_vert_loc = np.inf, None, []
        # Go through the test_atoms
        for a3 in test_atoms:
            if a3 in [a0, a1, my_a2]:
                continue
            # Calculate the vertex's
            my_vert_loc = calc_flat_vert([a0, a1, my_a2, a3])
            my_dist = calc_dist(my_vert_loc, my_circ[0])
            # Check the locations distance from the center of the atoms
            if my_dist < min_dist and [_.num for _ in [a0, a1, my_a2, a3]] not in net.vert_ndxs:
                my_vert = Vertex(atoms=[a0, a1, my_a2, a3], net=net, location=my_vert_loc, radius=0, flat_faces=True)
                min_dist = my_dist
        j += 1
        return my_vert


def ffind_site(edge_atoms, net, vn_1):
    # Get the atoms that should not ba a part of the new vertex
    edge_ndxs = [_.num for _ in edge_atoms]
    if vn_1 is None:
        vert_atom_ndxs = edge_ndxs
    else:
        vert_atom_ndxs = vn_1.ndx
    # Get the closest atoms around
    test_atoms = net.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], 5)
    # Set up the list of distances
    min_dist, my_vert, my_vert_ndx = np.inf, None, None
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
            # net.edges.append(Edge(net=net, atoms=edge_atoms, verts=[vn_1, net.verts[vert_ndx]], doublet=False, straight=True))
            return
        new_vert_ndxs.append(vert_ndx)
    # Get the planes made by the test atoms and one of the edge atoms
    for i in range(len(new_test_atoms)):
        an, vert_ndx = new_test_atoms[i], new_vert_ndxs[i]
        # Get the intersecting point of the edge and the plane made by a0 and an
        vert_loc = calc_flat_vert(edge_atoms + [an])
        # Check that the vertex is inside the network's box
        outside = False
        for j in range(3):
            if not (net.box[0][j] <= vert_loc[j] <= net.box[1][j]):
                outside = True
        if outside:
            continue
        # Find the distance between the new vert and the old vert
        vert_dist = calc_dist(vn_1.loc, vert_loc)
        if vert_dist < min_dist:
            my_vert = Vertex(atoms=edge_atoms + [an], net=net, location=vert_loc, radius=0, flat_faces=True)
            my_vert_ndx = vert_ndx
    # Check to see if the vertex is None
    if my_vert is not None:
        return my_vert, my_vert_ndx


# Find network function. Keeps searching the network until all verts are found
def ffind_verts(net, a0=None):
    # Calculate the total number of vertices
    tot_verts = 5 * len(net.atoms)
    # Find the first verified vertex
    if len(net.atoms) == 4:
        v0 = Vertex(atoms=net.atoms, net=net, location=calc_flat_vert(net.atoms), radius=0, flat_faces=True)
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
            if myVert is None:
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

    print(*[_.ndx for _ in net.verts])

from System.calcs import *
from System.Network.vertex import Vertex


# Find v0 function. Finds the first vertex in the network
def find_v0(net, a0=None):
    # If no a0 is given
    if a0 is None:
        # Find the middle sub_box of the set of boxes and
        mid = len(net.sub_boxes) // 2
        atoms = []
        inc = 1
        while not atoms:
            atoms = net.get_atoms([[mid, mid, mid]], inc)
            inc += 1
        a0 = atoms[-1]
    # Find the set of atoms with the minimum distance between surfaces
    min_dist = np.inf
    a1 = None
    inc = 0
    while not a1 and inc <= len(net.sub_boxes) + 1:
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
    while not a2 and inc <= len(net.sub_boxes) + 1:
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
    # Return the vertex
    return myVert


# Find site function. Currently, overkill, searching through all atoms for overlap and
def find_site(net, edge_atoms, vn_1=None):
    # Get the atoms that should not ba a part of the new vertex
    if vn_1 is None:
        vert_atoms = edge_atoms
    else:
        vert_atoms = vn_1.atoms
    verts = []
    min_rad = np.linalg.norm(np.array(net.box[0]) - np.array(net.box[1])) / 4
    for atom in net.atoms:
        if atom in vert_atoms:
            continue
        # If we have found the vertex before it is not the previous vertex return
        vert_found = False
        atom_ndxs = [net.atoms.index(atom1) for atom1 in edge_atoms + [atom]]
        for vert1 in net.verts:
            atom_ndxs.sort()
            if vert1.ndx == atom_ndxs:
                vert_found = True
        if vert_found:
            return
        vert = Vertex(edge_atoms + [atom], net=net)
        if vert.loc is None or vert.rad > min_rad:
            continue

        # Otherwise, find the indices of the sub-box for the vertex
        vi = int((vert.loc[0] - net.box[0][0]) / net.sub_box_size[0])
        vj = int((vert.loc[1] - net.box[0][1]) / net.sub_box_size[1])
        vk = int((vert.loc[2] - net.box[0][2]) / net.sub_box_size[2])
        # Get the number of boxes that an overlapping atom could possibly be away from the vertex sub-box
        atom_range = int(vert.rad / min(net.sub_box_size)) + int(5 / min(net.sub_box_size))
        overlap_test_atoms = net.get_atoms([[vi, vj, vk]], atom_range)
        overlap = False
        for atom2 in overlap_test_atoms:
            if atom2 in edge_atoms + [atom]:
                continue
            if calc_dist(atom2.loc, vert.loc) < atom2.rad + vert.rad:
                overlap = True
                break
        if not overlap:
            verts.append(vert)
    if len(verts) == 0:
        return
    elif len(verts) == 1:
        return verts[0]
    else:
        myVert = None
        min_rad = np.inf
        for vert in verts:
            if vert.rad < min_rad:
                myVert = vert
                min_rad = vert.rad
        return myVert


# Find network function. Keeps searching the network until all verts are found
def find_vertices(net, v0=None, i=0):
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
        tot_verts = len(net.verts) + (len(net.atoms) - i)
        percentage = int(len(net.verts) / tot_verts * 10000)/100
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
            # Add the vertex to the stack and the network
            vert_stack.append(myVert)
            net.verts.append(myVert)
            for atom in myVert.atoms:
                atom.verts.append(myVert)
    net.filter_verts()

from System.sys_calcs import *
from System.Network.vertex import Vertex

########################################################################################################################
"""Finding functions"""


# Doublet making function. Takes in a doublet and adds the two vertices to their respective places.
def doublet(verts, net):
    return verts[0]


# Find v0 function. Finds the first vertex in the network
def find_v0(net):
    # Find the center of mass of the atoms
    com = calc_atoms_com(net.atoms)
    # First choose an appropriate initial atom based of com proximity
    min_dist = np.inf
    a0 = None
    # Go through each atom determining if it is closer to the com
    for atom in net.atoms:
        # Set the new com distance
        com_dist = calc_dist(atom.loc, com)
        # If is less than the current closest atom's distance to the center of mass update the variables
        if com_dist < min_dist:
            min_dist = com_dist
            a0 = atom
    # Find the set of atoms with the minimum distance between surfaces
    min_dist = np.inf
    a1 = None
    # Go through each atom determining the atom with the minimum distance between it and a0's surfaces
    for atom in net.atoms:
        # Skip a0
        if atom == a0:
            continue
        # Set the new atom distances
        a_dist = calc_dist(a0.loc, atom.loc) - (a0.rad + atom.rad)
        # If the new atom distance is less than the previous minimum distance update the variables
        if a_dist < min_dist:
            min_dist = a_dist
            a1 = atom
    # Find the set of atoms with the minimum inscribed circle
    min_rad = np.inf
    a2 = None
    # Go through each other atom to determine the smallest circle that can be made with our 2 atoms and a third
    for atom in net.atoms:
        # Skip a0, a1
        if atom == a0 or atom == a1:
            continue
        # Calculate the circle made with the 3 atoms
        circ = calc_circ([a0, a1, atom])
        # If the radius of the circle that is made with the 3 atoms is smaller than the previous smallest radius replace
        if circ and abs(circ[0][1]) < min_rad:
            min_rad = abs(circ[0][1])
            a2 = atom
    # Find the set of atoms with the minimum inscribed sphere
    min_rad = np.inf
    myVert = None
    # Go through each other atom to determine the smallest possible inscribed sphere
    for atom in net.atoms:
        # Skip a0, a1, a2
        if atom == a0 or atom == a1 or atom == a2:
            continue
        # Get the vertex made from the atoms
        vert = Vertex(atoms=[a0, a1, a2] + [atom])
        print()
        # If the radius of the inscribed
        if vert.loc and vert.rad < min_rad:
            min_rad = vert.rad
            myVert = vert
    # Return the vertex
    return myVert


# Find site function. Takes in an edge and finds the only other vertex that does not overlap with other atoms
def find_site(edge_atoms, vn_1, net):
    # Instantiate the vertex
    myVert = None
    # Loop through the atoms to see if they create a vertex that doesn't overlap with any other atoms
    for atom in net.atoms:
        # This filters out any of the atoms in the edge or the remaining atom from the previous vertex
        if {atom}.issubset(vn_1.atoms):
            continue
        # Calculate the vertex with atom
        vert = Vertex(atoms=edge_atoms + [atom])
        if vert.loc is None:
            continue
        # Check if the vertex overlaps with any of the networks atoms
        overlap = False
        for a_test in net.atoms:
            if {a_test}.issubset(edge_atoms + [atom]):
                continue
            if round(calc_dist(a_test.loc, vert.loc) - (a_test.rad + vert.rad), 7) < 0:
                overlap = True
                break
        if not overlap:
            myVert = vert
            break
    return myVert


# Find network function. Keeps searching the network until all verts are found
def find_vertices(net):
    # Find the first vertex in the System
    v0 = find_v0(net)
    # Add v0 to the System
    net.verts.append(v0)
    # Set up the vertex stack
    vert_stack = [v0]
    # While the verts stack is not empty
    while vert_stack:
        # Running print statement giving an estimate for percentage of the network that has been created
        tot_verts = max(len(net.verts) + int(3*len(vert_stack)/4), 4 * len(net.atoms))
        percentage = int(len(net.verts) / tot_verts * 100)
        print("\rBuilding Network: ", percentage, "%", end='')
        # Get the vertex from the top of the stack
        vert = vert_stack.pop()
        # Set up the edge stack
        e_stack = [[[vert.atoms[i], vert.atoms[(i + 1) % 4], vert.atoms[(i + 2) % 4]], vert] for i in range(4)]
        # While the edge stack is not empty
        while e_stack:
            # Get the edge from the top of the stack
            edge = e_stack.pop()
            # Find the next site in the network
            myVert = find_site(edge[0], edge[1], net)
            # If the vertex is none continue
            if myVert is None:
                continue
            # If the vertex exists in the network add the vertex to the edge and move on to the next edge in the stack
            found_vert = check_vert(set(myVert.atoms), net.verts)
            if not found_vert:
                vert_stack.append(myVert)
                net.verts.append(myVert)


# Calculate direction function. Takes in a vertex and an edge and returns True if it is facing the center
def calc_dir(edge):
    # Grab the previous vertex
    vn_1 = edge.verts[0]
    # Find ak and copy it
    ak = None
    for atom in vn_1.atoms:
        if not {atom}.issubset(edge.atoms):
            ak = atom
    akp = ak
    # Find the direction toward the center of the edge
    r0 = [edge.loc[0] - vn_1.loc[0], edge.loc[1] - vn_1.loc[1], edge.loc[2] - vn_1.loc[2]]
    r0_mag = np.sqrt(r0[0]**2 + r0[1]**2 + r0[2]**2)
    r0_hat = [r0[0]/r0_mag, r0[1]/r0_mag, r0[2]/r0_mag]
    # Move the copy toward the center of the edge.
    akp.loc = [akp.loc[0] + r0_hat[0]*0.1, akp.loc[1] + r0_hat[1]*0.1, akp.loc[2] + r0_hat[2]*0.1]
    # Calculate the new vertex made by akp
    vkp = Vertex(edge.atoms + [akp])
    while not vkp:
        akp.loc = [akp.loc[0] - r0_hat[0]*0.01, akp.loc[1] - r0_hat[1]*0.1, akp.loc[2] - r0_hat[2]*0.1]
        vkp = Vertex(edge.atoms + [akp])
    # If the new inscribed sphere overlaps with ak, flip the direction of tang_hat
    if calc_dist(ak.loc, vkp.loc) - (ak.rad + vkp.rad) < 0:
        return False
    return True


# Calculate relative length function. Takes in 3 points and returns a float value for the relative distance
def calc_rel_dist(v0, v1, edge):
    # Grab the center
    c = np.array(edge.loc)
    # Find the distances between the 3 points
    r0, r1, r2 = np.linalg.norm(c - np.array(v0.loc)), np.linalg.norm(c - np.array(v1.loc)), \
                 np.linalg.norm(np.array(v0.loc) - np.array(v1.loc))
    # Cases 1 and 2: r0 > r1 > r2 and r0 > r2 > r1
    if r0 >= r1 and r0 > r2 and edge.dir:
        rel_dist = r2
    # Cases 3 and 4: r1 > r0 > r2 and r1 > r2 > r0
    elif r1 > r0 and r1 > r2 and not edge.dir:
        rel_dist = r2
    # Cases 5 and 6: r2 > r0 > r1 and r2 > r1 > r0
    elif r2 > r0 and r2 > r1 and edge.dir:
        rel_dist = r0 + r1
    # All other cases should not give a distance
    else:
        rel_dist = np.inf
    # Return the relative distance
    return rel_dist


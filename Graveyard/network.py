"""Imports"""
from System.system import *


########################################################################################################################
"""Gather information functions"""


# Find site function. Takes in an edge and finds the only other vertex that does not overlap with other atoms
def find_site1(self, edge_atoms, vn_1):
    # Instantiate the vertex
    myVert = None
    inc = 0
    # Loop through the atoms to see if they create a vertex that doesn't overlap with any other atoms
    while inc < len(self.sub_boxes):
        atoms = self.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], inc)
        for atom in atoms:
            # This filters out any of the atoms in the edge or the remaining atom from the previous vertex
            if {atom}.issubset(vn_1.atoms):
                continue
            # Calculate the vertex with atom and pass if the vertex location is None
            vert = Vertex(atoms=edge_atoms + [atom], net=self)
            if vert.loc is None:
                continue  ## This is where I implement Hu's Method
            # Find the box that the vertex would be in
            vi = int((vert.loc[0] - self.box[0][0]) / self.sub_box_size[0])
            vj = int((vert.loc[1] - self.box[0][1]) / self.sub_box_size[1])
            vk = int((vert.loc[2] - self.box[0][2]) / self.sub_box_size[2])
            # Get the set of atoms needed to check for the next vertex
            atoms = self.get_atoms([vert.atoms[0].box, vert.atoms[1].box, vert.atoms[2].box, vert.atoms[3].box,
                                    [vi, vj, vk]], int(vert.rad / min(self.sub_box_size)) + 3)
            # Check if the vertex overlaps with any of the networks atoms
            overlap = False
            # Loop over all atoms
            for a_test in atoms:
                if {a_test}.issubset(edge_atoms + [atom]):
                    continue
                if calc_dist(a_test.loc, vert.loc) - (a_test.rad + vert.rad) < 0:
                    overlap = True
                    break
            if not overlap:
                myVert = vert
                break
        inc += 1
    return myVert



# Find site function. When given an edge, this function returns the closest next vertex in the System
def find_site1(edge, net):
    # Get the edges location and radius
    circs = calc_circ(edge.atoms)
    if circs is None:
        return
    edge.loc, edge.rad = circs[0][0], abs(circs[0][1])
    # Get the edge's direction
    edge.dir = calc_dir(edge)
    # Grab the edge's previous vertex
    vn_1 = edge.verts[0]
    min_val = np.inf
    vn = None
    # Go through each atom in the System
    for atom in net.atoms:
        # Check for edge atoms and backwards atom by making sure atom is not in vn_1's list of atoms
        if {atom}.issubset(vn_1.atoms):
            continue
        # Set up an atoms list
        atoms = [atom] + edge.atoms
        # Calculate the vertex
        vert = calc_vert(atoms)
        if vert is None:
            continue
        # Find the relative distance between the verts along the edge
        dist = calc_rel_dist(vn_1, vert, edge)
        if dist < min_val:
            min_val = dist
            vn = vert
    return vn


# Find circle function. Finds the smallest circle between the two given atoms and every other atom and return that atom
def find_circle(a0, a1, net, num_checks=12):
    # Instantiate variables
    a2 = None
    neighbors = sortbyDist([a0, a1], net)
    rad = np.inf
    # Go through the num_checks closest atoms and find the smallest circle
    for atom in neighbors[1:num_checks]:
        # Check if any are the same as the new atom
        if a0 == atom or a1 == atom:
            continue
        # Calculate the radius of the circle made by the three atoms
        new_rad = calc_circ([a0, a1, atom])
        if new_rad:
            new_rad = new_rad[0][1]
        else:
            continue
        # Check the new radius against the smallest found and make it the smallest if it is
        if new_rad < rad:
            rad = new_rad
            a2 = atom
    # Return the atom found to have the smallest circle
    return a2


# Get initial vertex function. Finds an optimal starting vertex for the network.
def find_v0(net):
    # Grab the first atom in the network. This will be replaced later with an optimized
    a0 = net.atoms[0]
    # Find it's closest neighbor
    neighbors0 = sortbyDist([a0], net)
    a1 = neighbors0[0]
    # Find the smallest circle you can make with a0, a1 and a third atom
    a2 = find_circle(a0, a1, net)
    # Find the first site by choosing the smallest interstitial sphere that can be made with the 3 atoms and the 50
    # closest atoms
    neighbors2 = sortbyDist([a0, a1, a2], net)
    r = np.inf
    myVert, my_an = None, None
    # Go through the closest atoms to the triplet and find the smallest vertex that can be made with a neighbor
    for an in neighbors2[:15]:
        if an.loc == a0.loc or an.loc == a1.loc or an.loc == a2.loc:
            continue
        # Calculate the vertex and check if None
        vert = calc_vert([a0, a1, a2, an])
        # Don't worry about None vert types
        if vert is None:
            continue
        # Check if the vn has a smaller radius than the current smallest radius
        if vert.rad < r:
            r = vert.rad
            myVert = vert
            my_an = an
    # Add connections to the network
    net.verts.append(myVert)
    return myVert


########################################################################################################################
"""Recursive network finding function"""


# Find edges function. Recursively traces out the network and records vertex locations,
def find_edges(vertex, net):
    # Create 4 edges with the 4 combinations of atoms that can be created
    for i in range(4):
        myEdge = Edge([vertex.atoms[i], vertex.atoms[(i + 1) % 4], vertex.atoms[(i + 2) % 4]], [vertex])
        etest = check_edge(myEdge, net.edges)
        if etest:
            vertex.edges.append(etest)
        else:
            vertex.edges.append(myEdge)
    # Create the edge objects or grab them from the network and connect them
    for edge in vertex.edges:
        # Check to see if the edge exists. If it does move on the next edge in the vertex
        net_edge = check_edge(set(edge.atoms), net.edges)
        if net_edge:
            continue
        # Create a vertex
        vn = calc_edge1(edge, net)
        # If the vertex is None give the edge a None vertex and continue to the next edge
        if vn is None:
            edge.verts.append(Vertex([np.inf, np.inf, np.inf], np.inf))
            net.edges.append(edge)
            continue
        # Check the vertex to see if it exists in the network
        net_vert = check_vert(vn, net)
        # If it does, add the vertex to the edge and the edge to the network
        if net_vert:
            edge.verts.append(net_vert)
            net.edges.append(edge)
            return
        # If both the edge and the vertex do not exist in the network, we have a true new site
        else:
            edge.verts.append(vn)
            net.edges.append(edge)
            net.verts.append(vn)
            find_edges(vn, net)
    return


########################################################################################################################


# Build Network function. Takes in a System, runs as a shell for the recursive next_site function and returns a Network
def build_network(mySys):
    # Find the first vertex
    v0 = find_v01(mySys.myNet)
    # Initiate the recursive network finding algorithm on the network and the first vertex
    find_edges(v0, mySys.myNet)
    return mySys



# Calculate vertex function. Takes in 4 atoms, calculates the center and radius of the inscribed sphere and returns them
def calc_vert(atoms):
    # The real location and radius of the base sphere
    l1, R1 = atoms[0].loc, atoms[0].rad
    # Set the radii and x, y, z values for the 3 spheres
    R2, R3, R4 = atoms[1].rad, atoms[2].rad, atoms[3].rad
    x2, y2, z2 = atoms[1].loc[0] - l1[0], atoms[1].loc[1] - l1[1], atoms[1].loc[2] - l1[2]
    x3, y3, z3 = atoms[2].loc[0] - l1[0], atoms[2].loc[1] - l1[1], atoms[2].loc[2] - l1[2]
    x4, y4, z4 = atoms[3].loc[0] - l1[0], atoms[3].loc[1] - l1[1], atoms[3].loc[2] - l1[2]

    # Calculate our System of linear equations coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (R2 - R1), R1 ** 2 - R2 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (R3 - R1), R1 ** 2 - R3 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    a3, b3, c3, d3, f3 = 2 * x4, 2 * y4, 2 * z4, 2 * (R4 - R1), R1 ** 2 - R4 ** 2 + x4 ** 2 + y4 ** 2 + z4 ** 2

    A, B, C, d, f = [a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f2]

    # Calculate the ranks of the matrices
    ABC_rank = np.linalg.matrix_rank([A, B, C])
    m_rank = np.linalg.matrix_rank([A, B, C, d])
    f_rank = np.linalg.matrix_rank([A, B, C, d, f])

    # Calculate the F values
    F, F10, F11, F20, F21, F30, F31 = a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1, \
                                      b1*c2*f3 - b1*c3*f2 - b2*c1*f3 + b2*c3*f1 + b3*c1*f2 - b3*c2*f1, \
                                      -b1*c2*d3 + b1*c3*d2 + b2*c1*d3 - b2*c3*d1 - b3*c1*d2 + b3*c2*d1, \
                                      -a1*c2*f3 + a1*c3*f2 + a2*c1*f3 - a2*c3*f1 - a3*c1*f2 + a3*c2*f1, \
                                      a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1, \
                                      a1*b2*f3 - a1*b3*f2 - a2*b1*f3 + a2*b3*f1 + a3*b1*f2 - a3*b2*f1, \
                                      -a1*b2*d3 + a1*b3*d2 + a2*b1*d3 - a2*b3*d1 - a3*b1*d2 + a3*b2*d1
    # Catch for F = 0.
    if F == 0:
        return
    # Instantiate our root arrays
    xs, ys, zs, Rs = [], [], [], []
    verts = []
    # Case 1:
    if ABC_rank == 3 and m_rank == 3 and f_rank == 3:
        # Calculate the radius polynomial coefficients
        a = ((F11 ** 2 + F21 ** 2 + F31 ** 2) / F ** 2) - 1
        b = (2 * (F10 * F11 + F20 * F21 + F30 * F31) / F ** 2) - 2 * R1
        c = ((F10 ** 2 + F20 ** 2 + F30 ** 2) / F ** 2) - R1 ** 2
        # If the discriminant is positive, find the real positive roots of the quadratic
        if -4*a*c + b**2 > 0:
            Rs = [R for R in np.roots([a, b, c]) if np.isreal(R) and R > 0]
        # Instantiate the verts array
        verts = []
        # Go through each radius and calculate the vertex
        for R in Rs:
            x, y, z = F10/F + R*F11/F, F20/F + R*F21/F, F30/F + R*F31/F
            # Move the vertex back to the actual location of the atoms
            verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))

    # Case 2:
    elif ABC_rank == 2 and m_rank == 3 and f_rank == 3:
        # Case 2 subcases:
        # Case 2.1
        if np.linalg.matrix_rank([A, B, d]) == 3:
            # Calculate the z value polynomial coefficients
            a = F**2 + F11**2 + F21**2 - F31**2
            b = 2*(F10*F11 + F20*F21 - F30*F31 - F*F31*R1)
            c = F10**2 + F20**2 - (F30 + F*R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                zs = [z for z in np.roots([a, b, c]) if np.isreal(z) and z > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for z in zs:
                x, y, R = F10 / F + z * F11 / F, F20 / F + z * F21 / F, F30 / F + z * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))

        # Case 2.2
        elif np.linalg.matrix_rank([A, d, C]) == 3:
            # Calculate the z value polynomial coefficients
            a = F ** 2 + F11 ** 2 - F21 ** 2 + F31 ** 2
            b = 2 * (F10 * F11 - F20 * F21 + F30 * F31 - F * F31 * R1)
            c = F10 ** 2 + F30 ** 2 - (F20 + F * R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                ys = [y for y in np.roots([a, b, c]) if np.isreal(y) and y > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for y in ys:
                x, R, z = F10 / F + y * F11 / F, F20 / F + y * F21 / F, F30 / F + y * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))

        # Case 2.3
        elif np.linalg.matrix_rank([d, B, C]):
            # Calculate the z value polynomial coefficients
            a = F ** 2 + F11 ** 2 + F21 ** 2 - F31 ** 2
            b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * R1)
            c = F10 ** 2 + F20 ** 2 - (F30 + F * R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                xs = [x for x in np.roots([a, b, c]) if np.isreal(x) and x > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for x in xs:
                R, y, z = F10 / F + x * F11 / F, F20 / F + x * F21 / F, F30 / F + x * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))
    # If no verts are found return None
    if not verts:
        return
    else:
        if len(verts) == 2 and verts[0].rad > verts[1].rad:
            return verts[1]
        return verts[0]



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



# Move sphere function. Takes in a location, an Atom object and a direction and updates the Atom's location
def move(loc, atom, to_home=False):
    # Change whether we are adding or subtracting the location to the sphere's location.
    d = 1
    if not to_home:
        d = -1
    # Update the atom's location
    atom.loc[0] = atom.loc[0] + d * loc[0]
    atom.loc[1] = atom.loc[1] + d * loc[1]
    atom.loc[2] = atom.loc[2] + d * loc[2]


# Sort by distance function. Sorts all atoms in the System by distance from COM of given atoms
def sortbyDist(atoms, net):
    # Find the point closest to each of the atoms
    loc = [0, 0, 0]
    for i in range(len(atoms)):
        loc = loc[0] + atoms[i].loc[0], loc[1] + atoms[i].loc[1], loc[2] + atoms[i].loc[2]
    loc = [loc[0]/len(atoms), loc[1]/len(atoms), loc[2]/len(atoms)]
    # Initialize the lists
    dist_list, atom_list = [], []
    # Go through all the atoms in the molecules
    for atom in net.atoms:
        # Don't include the atoms in our list of atom
        if atom in atoms:
            continue
        # Get the distance between the atoms and subtract their radii
        dist = calc_dist(loc, atom.loc) - atom.rad
        dist_list.append(dist)
        atom_list.append(atom)
    # Selection sort the atom list based off their distances from the point
    for i in range(len(dist_list)):
        low_in = i
        for j in range(i+1, len(dist_list)):
            if dist_list[low_in] > dist_list[j]:
                low_in = j
                dist_list[i], dist_list[low_in] = dist_list[low_in], dist_list[i]
                atom_list[i], atom_list[low_in] = atom_list[low_in], atom_list[i]
    # Return a list with the length specified
    return atom_list



 # Find subnetworks method. Divides the network up into smaller networks, finds their vertices and combines them
def find_subnets(self):
    # Instantiate the vertices
    self.verts = []
    # This gets us to average about 60 atoms per medium box
    n = max(int(np.cbrt(len(self.atoms) // 30)), 2)
    # The range to search for new vertices
    rnge = len(self.sub_boxes) // n
    # Get the atoms in the
    med_boxes = [[[[] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    # Set up the network list and the list of the atom indices from the main network for reference later
    test_nets, test_nets_real_ndxs = [], []
    # Calculate the overlap needed to prevent missing vertices
    overlap = int((self.max_vert + self.max_atom_rad) / min(self.sub_box_size)) + 1
    # Go through each of the boxes in the medium_boxes matrix to create networks
    for i in range(n):
        for j in range(n):
            for k in range(n):
                # Starting from the middle of the range search out the atoms within 0.75 + the increment
                med_boxes[i][j][k] = self.get_atoms([[int((i + .5) * rnge), int((j + .5) * rnge),
                                                      int((k + .5) * rnge)]], int(rnge // 2) + overlap)
                # Skip the net if there are no atoms in the quadrant
                if len(med_boxes) == 0:
                    continue
                # Get the increment to spread each net out by
                inc = overlap
                # Make sure there are enough atoms to do a valid search
                while len(med_boxes[i][j][k]) < 30:
                    # Starting from the middle of the range search out the atoms within 0.75 + the increment
                    med_boxes[i][j][k] = self.get_atoms([[int((i + .5) * rnge), int((j + .5) * rnge),
                                                          int((k + .5) * rnge)]], int(rnge // 2) + inc)
                    inc += 1
                # Set up the indices' tracker list
                new_atoms, old_ndxs = [], []
                # Create a copy of each atom in the network
                for atom in med_boxes[i][j][k]:
                    new_atoms.append(Atom(atom.loc, atom.rad))
                    old_ndxs.append(self.atoms.index(atom))
                # Create the network and add it to the list of subnetworks
                test_nets.append(Network(sys=self.sys, atoms=new_atoms, box_size=1.1, min_dist=self.min_dist,
                                         max_vert=self.max_vert, sol_verts=self.sol_verts))
                # Store the atoms' real indices
                test_nets_real_ndxs.append(old_ndxs)
    # Find the vertices for each of the networks
    for i in range(len(test_nets)):
        # Get the subnetwork
        net2 = test_nets[i]
        # Sort the atoms in the network
        net2.sort_atoms()
        # Set the name of the subnetwork and print the progress
        net2.name = "Subnetwork " + str(i + 1) + "/" + str(len(test_nets)) + " of " + self.name
        print("\rFinding " + net2.name + " - " + str(len(net2.atoms)) + " atoms", end="")
        # Find the vertices
        net2.find_verts(n)
    # Go through each of the vertices found in each of the subnetworks filtering out repeat/large verts + doublets
    for i in range(len(test_nets)):
        # Grab the subnetwork
        net2 = test_nets[i]
        # Sort through the vertices in each network
        for j in range(len(net2.verts)):
            # Grab the vertex
            vert = net2.verts[j]
            # Set the vertex's actual atoms and use them to set the vertex's index, then sort it
            vert.atoms = [self.atoms[test_nets_real_ndxs[i][ndx]] for ndx in vert.ndx]
            vert.ndx = [self.atoms.index(atom) for atom in vert.atoms]
            vert.ndx.sort()
            # Look for the vertex in the network
            v_ndx = search_verts(self.vert_ndxs, vert.ndx)
            # If found, verify it is not a doublet
            if 0 == len(self.vert_ndxs) or len(self.vert_ndxs) <= v_ndx:
                vert.myNet = self
                self.vert_ndxs.append(vert.ndx)
                self.verts.append(vert)
            # If the indices match we need to check if it is a doublet
            elif vert.ndx == self.vert_ndxs[v_ndx]:
                # If the location is different, and we haven't indicated the vertex as a doublet already add it
                if [round(vert.loc[k], 7) for k in range(3)] != \
                        [round(self.verts[v_ndx].loc[k], 7) for k in range(3)] \
                        and not self.verts[v_ndx].doublet:
                    # Add the doublet
                    self.verts[v_ndx].doublet = True
                    self.verts[v_ndx].loc2, self.verts[v_ndx].rad2 = vert.loc, vert.rad
            # If this is a novel vertex, we need to verify it
            else:
                if verify_site(vert, self):
                    # Add the vertex to the system
                    vert.myNet = self
                    self.vert_ndxs.insert(v_ndx, vert.ndx)
                    self.verts.insert(v_ndx, vert)



############################################ Filter System #############################################################


def filter_verts(net):
    """
    Filters out vertices that are repeated, outside the box or larger than the max vertex value

    :param net: network of unfiltered vertices

    """
    # Set up a list of vertex ndxs and vertices
    vert_ndxs = []
    verts = []
    # Check to see if no vertices have been made
    if net.verts is None:
        print("No vertices to filter")
        return

    # Go through the vertices
    for i in range(len(net.verts)):

        # Boolean for whether the vertex is inside the box or not
        loc_in_box = True
        vert = net.verts[i]
        # Check for None vertices
        if vert.loc is None:
            continue
        # Check if the vertex is inside the box
        for j in range(3):
            if vert.loc[j] < net.box[0][j] or vert.loc[j] > net.box[1][j]:
                loc_in_box = False

        # Search the list of vertices for the vertex
        if vert.ndx not in vert_ndxs and loc_in_box:
            vert_ndxs.append(net.verts[i].ndx)
            verts.append(net.verts[i])
            if net.verts[i].doublet is not None:
                verts.append(net.verts.append())

    # Set the networks vertices
    net.verts = verts
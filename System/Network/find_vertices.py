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


# Find site function. Takes in an edge and finds the only other vertex that does not overlap with other atoms
def find_site(net, edge_atoms, vn_1=None):
    # Instantiate the vertex
    if vn_1 is None:
        old_atoms = []
    else:
        old_atoms = vn_1.atoms
    myVert = None
    inc = 0
    min_vert = np.inf
    # Loop through the atoms to see if they create a vertex that doesn't overlap with other atoms or whole system
    while inc <= 7:
        # If no vert has been found yet, keep expanding
        if myVert is None or myVert.loc is None:
            # Keep adding boxes around the atoms to check
            atoms = net.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], inc)
            # Go through the atoms in the surrounding boxes and find the smallest vertex that can be created
            for atom1 in atoms:
                # This filters out any of the atoms in the edge or the remaining atom from the previous vertex
                if atom1 in old_atoms:
                    continue
                # Calculate the vertex with atom and pass if the vertex location is None
                vert = Vertex(atoms=edge_atoms + [atom1], net=net)
                # I need to fix Vertex to get to a point where I don't need this
                if vert is None or vert.loc is None:
                    continue
                # Sniff out the smallest vertex that can be made in the box and store it
                if vert.rad < min_vert:
                    # Replace the variables
                    min_vert = vert.rad
                    myVert = vert
        # If no vertex can be made restart the search with a larger set of boxes
        if myVert is None or myVert.loc is None:
            inc += 1
            continue
        # Find the box that the vertex would be in
        vi = int((myVert.loc[0] - net.box[0][0]) / net.sub_box_size[0])
        vj = int((myVert.loc[1] - net.box[0][1]) / net.sub_box_size[1])
        vk = int((myVert.loc[2] - net.box[0][2]) / net.sub_box_size[2])
        # Any atom that can overlap with this vertex is within the 'vertex's radius over the smallest box length'
        # plus the 'maximum radius of an atom over the smallest box length' # of boxes away from the vert.loc box
        atoms = net.get_atoms([[vi, vj, vk]],
                              int(myVert.rad / min(net.sub_box_size)) + int(5 / min(net.sub_box_size)) + 2)
        # Set up an overlap variable
        overlap = False
        min_rad = np.inf
        # Test the atoms in the new atom list to see if they overlap with the vertex
        for atom in atoms:
            # If the atom is one of the vert atoms move on
            if atom in myVert.atoms:
                continue
            # If the distance between the vertex and the atom is less than their radii create a new vertex and reset
            if myVert.loc and calc_dist(atom.loc, myVert.loc) - (atom.rad + myVert.rad) < 0:
                vert = Vertex(edge_atoms + [atom], net=net)
                overlap = True
                # Set the new vertex to test as the one that can be made with the atom that overlaps with it
                if myVert.rad < min_rad:
                    myVert = vert
                    min_rad = myVert.rad
        # Check to see if the vertex had no overlaps --> We found the vertex!
        if not overlap:
            return myVert
        # If we have made it this far increment the counter and keep expanding
        inc += 1


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
            edge, vert = e_stack.pop()
            # Find the next site in the network
            myVert = find_site(net, edge, vert)
            # If the vertex is none continue
            if myVert is None:
                continue
            # If the vertex exists in the network add the vertex to the edge and move on to the next edge
            found_vert = check_vert(set(myVert.atoms), net.verts)
            if not found_vert:
                vert_stack.append(myVert)
                net.verts.append(myVert)
                for atom in myVert.atoms:
                    atom.verts.append(myVert)

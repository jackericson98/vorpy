# if b0[0] <= loc[0] <= b1[0] and b0[1] <= loc[1] <= b1[1] and b0[2] <= loc[2] <= b1[2]:
#     self.loc, self.rad = loc, rad
# else:
#     return
# # Otherwise, choose the first
# else:
# loc, rad = verts[0][0], abs(verts[0][1])
# # Check to see if the vertex is in the box or not
# if b0[0] <= loc[0] <= b1[0] and b0[1] <= loc[1] <= b1[1] and b0[2] <= loc[2] <= b1[2]:
#     self.loc, self.rad = loc, rad
# else:
#     Worst
# case
# scenario
# try the Hu Method
# pass
# loc = self.fv2()
# if len(loc) > 0:
#     self.loc = loc
#     self.rad = np.linalg.norm(self.loc - self.atoms[0].loc) - self.atoms[0].rad

# Find site function. Takes in an edge and finds the only other vertex that does not overlap with other atoms
def find_site1(net, edge_atoms, vn_1=None):
    # Get the atoms that should not ba a part of the new vertex
    if vn_1 is None:
        vert_atoms = edge_atoms
    else:
        vert_atoms = vn_1.atoms
    # Instantiate the vertex, incrementer, minimum radius and overlap check
    min_rad = np.linalg.norm(np.array(net.box[0]) - np.array(net.box[1])) / 4
    myVert, inc, overlap = None, 0, False
    # Go through larger and larger search area looking for a vertex
    while myVert is None and inc < len(net.sub_boxes):
        # Grab atoms from the cells surrounding the edge atoms
        vert_test_atoms = net.get_atoms([edge_atoms[0].box, edge_atoms[1].box, edge_atoms[2].box], inc)
        # Go through each atom in the surrounding cells
        for atom in vert_test_atoms:
            # If the atom is a part of the vertex go to the next atom
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
            # Create the vertex to test against
            vert = Vertex(edge_atoms + [atom], net=net)
            if vert.rad is None or vert.rad > min_rad:
                continue
            # Otherwise, find the indices of the sub-box for the vertex
            vi = int((vert.loc[0] / net.sub_box_size[0]) - net.box[0][0])
            vj = int((vert.loc[1] / net.sub_box_size[1]) - net.box[0][1])
            vk = int((vert.loc[2] / net.sub_box_size[2]) - net.box[0][2])
            # Get the number of boxes that an overlapping atom could possibly be away from the vertex sub-box
            atom_range = int(vert.rad / min(net.sub_box_size)) + int(5 / min(net.sub_box_size)) + 1
            overlap_test_atoms = net.get_atoms([[vi, vj, vk]], atom_range)
            # Set the overlap to False so that we have an innocent assumption for our vertex
            for atom1 in overlap_test_atoms:
                # Skip any atoms in the vertex
                if atom1 in vert.atoms:
                    continue
                # Check to see if there is any overlap with the check atom
                if calc_dist(atom1.loc, vert.loc) <= atom1.rad + vert.rad:
                    overlap = True
                    break
            # Check to see if the vertex's radius exists and is less than the minimum radius
            if vert.rad < min_rad and not overlap:
                myVert = vert
                min_rad = vert.rad
        inc += 1
    return myVert
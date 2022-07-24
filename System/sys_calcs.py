import numpy as np


# Sort by distance function. Sorts all atoms in the System by distance from COM of given atoms
def sortbyDist(atoms, net, length=None):
    # If the length of the returned list is not specified return the whole list
    if length is None:
        length = len(net.atoms)
    # Find the point closest to each of the atoms
    loc = [0, 0, 0]
    for i in range(len(atoms)):
        f = i + 1
        loc = loc[0] + atoms[i].loc[0] / f, loc[1] + atoms[i].loc[1] / f, loc[2] + atoms[i].loc[2] / f
    # Initialize the lists
    dist_list = []
    atom_list = []
    # Go through all the atoms in the molecules
    for atom2 in net.atoms:
        # Don't include the atoms in our list of atom
        if atom2 in atoms:
            continue
        # Get the distance between the atoms and subtract their radii
        dist = calc_dist(loc, atom2.loc) - atom2.rad
        dist_list.append(dist)
        atom_list.append(atom2)
    # Selection sort the atom list based off their distances from the point
    for i in range(len(dist_list)):
        low_in = i
        for j in range(i+1, len(dist_list)):
            if dist_list[low_in] > dist_list[j]:
                low_in = j
                dist_list[i], dist_list[low_in] = dist_list[low_in], dist_list[i]
                atom_list[i], atom_list[low_in] = atom_list[low_in], atom_list[i]

    # Return a list with the length specified
    return atom_list[:length]


# Calculate distance function. Takes in 2 points and returns the distance between them
def calc_dist(l1, l2):
    d = np.sqrt((l1[0]-l2[0])**2+(l1[1]-l2[1])**2+(l1[2]-l2[2])**2)
    return d


# Calculate center of mass function. Takes in a set of points and returns the coordinates of the com
def calc_com(atoms):
    # Set the running sum for the x, y, z values to 0
    xtot, ytot, ztot = 0, 0, 0
    for atom in atoms:
        xtot = xtot + atom.loc[0]
        ytot = ytot + atom.loc[1]
        ztot = ztot + atom.loc[2]
    return xtot/len(atoms), ytot/len(atoms), ztot/len(atoms)


"""System checks"""


# Check surf function. Takes in a set of atoms and a list of surfs and returns the corresponding surf or None if no surf
def check_surf(s_atoms, surf_list):
    # Go through each surf in the surf list
    for surf in surf_list:
        # Check if the given atoms correspond to the atoms in the surf
        if s_atoms.issubset(surf.atoms):
            # Return the surf
            return surf
    return


# Check edge function. Takes in a set of atoms and a list of edges and returns the corresponding edge or None if no edge
def check_edge(e_atoms, edge_list):
    # Go through each edge in the edge list
    for edge in edge_list:
        # Check if the given atoms correspond to the atoms in the edge
        if e_atoms.issubset(edge.atoms):
            # Return the edge
            return edge
    return


# Check vert function. Takes in a set of atoms and a list of verts and returns the corresponding edge or None if no vert
def check_vert(v_atoms, vert_list):
    # Go through each edge in the edge list
    for vert in vert_list:
        # Check if the given atoms correspond to the atoms in the edge
        if v_atoms.issubset(vert.atoms):
            # Return the edge
            return vert
    return

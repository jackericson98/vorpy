import numpy as np
from System.system import System
from System.Network.network import Network


def ffind_near_atoms(net, a0, max_dist=20):
    # Get the closest atoms to a0
    max_atom_dist = 0
    dists = []
    a0_array = np.array(a0.loc)
    max_inc = int(max_dist/min(net.sub_box_size)) + 1
    # Get the atoms
    my_atoms = net.get_atoms([a0.box], reach=max_inc)
    # Get the atoms distance from a0
    for atom in my_atoms[len(dists):]:
        # Calculate the distance between the atoms
        my_dist = np.sqrt(sum(np.square(a0_array - np.array(atom.loc)))) - (a0.rad + atom.rad)
        # Replace the maximum distance if needed
        if my_dist > max_atom_dist:
            max_atom_dist = my_dist
            # Check that the atom's distance isn't larger than the specified maximum distance
            if my_dist < max_dist:
                dists.append(my_dist)
        else:
            dists.append(my_dist)
    # Sort the atoms
    sorted_atoms = [[x, _] for _, x in sorted(zip(dists, my_atoms), key=lambda pair: pair[0]) if x.num != a0.num]
    # Make sure the return doesn't trigger an index error
    return sorted_atoms


def make_cell(atom, surr_atoms):
    my_center = np.array(atom.loc)
    rns = []
    # Create the rn vectors
    for surr_atom in surr_atoms:
        r = my_center - np.array(surr_atom.loc)
        rns.append(r / np.linalg.norm(r))
    # Create a list of maximum dot products for each surr atom
    atom_dots = [[] for _ in range(len(surr_atoms))]
    atom_atoms = [[] for _ in range(len(surr_atoms))]
    for i in range(len(surr_atoms)):
        # Get the surrounding atoms and their rn's
        surr_atom, surr_atom_rn = surr_atoms[i], rns[i]
        # Go through the surfaces near this one
        for j in range(len(surr_atoms)):
            if i == j or surr_atoms[j] in atom_atoms[i]:
                continue
            # Calculate the dot product between the current surface and
            my_dot = np.dot(rns[i], rns[j])
            k = 0
            # Compare this dot product to the others and place if in the correct order
            while k < len(atom_dots[i]) and my_dot < atom_dots[i][k]:
                k += 1
            atom_dots[i].insert(k, my_dot)
            atom_atoms[i].insert(k, surr_atoms[j])
            # Find where to insert the other atom's info for this surface
            m = 0
            # Compare the dot product to the other dot products in the surrounding atom's list
            while m < len(atom_dots)

# Find sites function. Takes in an atom, a list of verified surfaces and a bank of additional atoms to choose from
def ffind_sites(net, atom, surfs, atoms_bank):
    max_vert_dist = 0
    # Keep cycling until the closest atom is further than maximum vertex




def ffind_verts(net, max_dist=20):
    # Create the network's ledger of surfaces and surface distances
    net_surfs = [[] for _ in range(len(net.atoms))]
    net_surfs_atoms = []
    net_surf_dists = [[] for _ in range(len(net.atoms))]
    # Go through the atoms in the network
    for i in range(len(net.atoms)):
        # Set up the current atom's variable
        my_atom = net.atoms[i]
        # Get the sorted closest max_surfs number of atoms, atoms
        surr_atoms = ffind_near_atoms(net=net, a0=my_atom, max_dist=max_dist)
        net_surfs_atoms.append([[_[0].num, _[1]] for _ in surr_atoms])
        # Add the atoms to the net surfs list
        for j in range(4):
            # Get the check atom variable
            check_atom, check_dist = surr_atoms[j]
            # Pass if we've searched this atom before
            if check_atom.num in net_surfs[i]:
                continue
            # Find where to insert the surface
            k = 0
            while k < len(net_surf_dists[i]) and check_dist > net_surf_dists[i][k]:
                k += 1
            # Insert the surface's other atom's index and distance
            net_surfs[i].insert(k, check_atom.num)
            net_surf_dists[i].insert(k, check_dist)
            # Find where to insert my_atom into the close atom's list
            m = 0
            while m < len(net_surf_dists[check_atom.num]) and check_dist > net_surf_dists[check_atom.num][m]:
                m += 1
            net_surfs[check_atom.num].insert(m, my_atom.num)
            net_surf_dists[check_atom.num].insert(m, check_dist)
    # Filter the atoms already found from the net_surfs_atoms lists
    net_atoms_bank = []
    for i in range(len(net.atoms)):
        net_atoms_bank.append([_ for _ in net_surfs_atoms if _[0] not in net_surfs[i]])

    ### At this point we have a list of surrounding surfaces for the atoms and a sorted list of potential surfaces
    for i in range(len(net.atoms)):
        # Find the sites for the current atom
        ffind_sites(net, net.atoms[i], net_surfs[i], net_atoms_bank[i])

if __name__ == '__main__':
    my_sys = System(file='C:/Users/jacke/PycharmProjects/vorpy/Data/test_data/Na5.pdb', root_dir='C:/Users/jacke/PycharmProjects/vorpy')
    my_sys.net = Network(sys=my_sys, atoms=my_sys.atoms)
    my_sys.net.sort_atoms()
    ffind_verts(my_sys.net)



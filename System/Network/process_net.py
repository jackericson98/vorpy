from System.calcs import *

# The point of this file is to pre-process all the information for the system


def find_sol_layers(net, mol=None):
    # We want a function that creates full solute layers
    # This variable holds the atoms for their respective layers
    layers_atoms = []
    # This variable holds the surfaces associated with the layers
    layers = []
    # Holds the current set of atoms being surrounded
    current_layer, next_layer = None, None
    sol_atoms = net.sys.sol.copy()
    if sol_atoms is None:
        return
    i = 0
    # Keep adding layers until the sol atoms are gone
    while len(sol_atoms) > 0:
        # Set the current layer to what was the next layer
        current_layer, next_layer = next_layer, []
        # If the current layer is none, this is the first layer and we are calculating the molecule's sol layer
        if current_layer is None:
            # If no molecule is provided go off of the system's molecule list
            if mol is None:
                current_layer = []
                for my_mol in net.sys.mols:
                    if net.sys.mol_names[net.sys.mols.index(my_mol)].lower() == 'sol':
                        continue
                    current_layer += my_mol
            else:
                current_layer = mol
        # Set up the shell storage for the surfaces
        shell = []
        # Check each surface to see if it is in the current layer
        for surf in net.surfs:
            if surf.atoms[0] in current_layer and surf.atoms[1] in sol_atoms:
                my_atom = sol_atoms.pop(sol_atoms.index(surf.atoms[1]))
                next_layer.append(my_atom)
                shell.append(surf)
        # Add the correct variables
        layers.append(shell)
        layers_atoms.append(current_layer)

        i += 1
    # Add the last layer of atoms
    layers.append(next_layer)
    return layers, layers_atoms




def calc_sas(net):
    # Go through every surface and calculate it's surface area
    for surf in net.surfs:
        surf.sa = calc_sa(surf)


def calc_vols(net):
    # Go through each atom in the system and calculate its volume
    for atom in net.atoms:
        calc_vol(atom)


def calc_curve(net):
    # Go through each surface and calculate the curvature
    for surf in net.surfs:
        pass

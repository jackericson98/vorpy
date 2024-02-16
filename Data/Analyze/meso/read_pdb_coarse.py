from System.sys_funcs.calcs.sorting import get_radius
from System.sys_objs.atom import make_atom
from System.sys_objs.residue import Residue
from System.sys_objs.chain import Chain, Sol
import os.path as path
import numpy as np
from pandas import DataFrame


def read_pdb(file_address):
    """
    Interprets pdb data into a system of atom objects
    :param sys: System to add the pdb information to
    :param file: .pdb file to be added to the system
    :return: list of tuples of locations and radii
    """

    # Get the file information and make sure to close the file when done
    with open(file_address, 'r') as f:
        my_file = f.readlines()
    # Set up the atom and the data lists
    atoms, data, atom_count = [], [], 0
    chains, residues = [], []
    chains, resids = {}, {}
    # Check if the file is a foam file
    if my_file[0].split()[1] == 'coarsify':
        coarse = True
    # Go through each line in the file and check if the first word is the word we are looking for
    for i in range(len(my_file)):
        # Check to make sure the line isn't empty
        if len(my_file[i]) == 0:
            continue
        # Pull the file line and first word
        line = my_file[i]
        word = line[:4].lower()
        # Check to see if the line is an atom line
        if line and word == 'atom':  # Check if the line starts with atom
            # Check for the "m" situation
            if line[76:78] == ' M':
                continue
            name = line[12:16]
            name.strip()
            res_seq = line[22:26]
            if line[22:26] == '    ':
                res_seq = 0
            # Create the atom
            atom = make_atom(location=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), system=sys,
                             element=line[76:78].strip(), res_seq=int(res_seq), name=name, seg_id=line[72:76],
                             index=atom_count)

            atom_count += 1
            # If no chain is specified, set the chain to 'None'
            res_str, chain_str = line[17:20].strip(), line[21]
            if chain_str == ' ':
                if res_str.lower() in {'sol', 'hoh', 'sod'}:
                    chain_str = 'SOL'
                elif res_str.lower() in {'cl', 'mg', 'na', 'k', 'ion'} and 'SOL' in chains:
                    chain_str = 'SOL'
                else:
                    chain_str = 'A'
            # Create the chain and residue dictionaries
            res_name, chn_name = line[17:20] + str(atom['res_seq']), chain_str
            # If the chain has been made before
            if chn_name in chains:
                # Get the chain from the dictionary and add the atom
                my_chn = chains[chn_name]
                my_chn.add_atom(atom['num'])
                atom['chn'] = my_chn
            # Create the chain
            else:
                # If the chain is the sol chain
                if res_str.lower() == 'sol':
                    my_chn = Sol(atoms=[atom['num']], residues=[], name=chn_name, sys=sys)
                    sol = my_chn
                # If the chain is not sol create a regular chain object
                else:
                    my_chn = Chain(atoms=[atom['num']], residues=[], name=chn_name, sys=sys)
                    chains.append(my_chn)
                # Set the chain in the dictionary and give the atom it's chain
                chains[chn_name] = my_chn
                atom['chn'] = my_chn

            # Assign the atoms and create the residues
            if res_name in resids:
                my_res = resids[res_name]
                my_res.atoms.append(atom['num'])
            else:
                my_res = Residue(sys=sys, atoms=[atom['num']], name=res_str, sequence=atom['res_seq'], chain=atom['chn'])
                atom['chn'].residues.append(my_res)
                resids[res_name] = my_res
                if res_str.lower() != 'sol':
                    sys.residues.append(my_res)
                else:
                    sys.sol.residues.append(my_res)
            # Assign the residue to the atom
            atom['res'] = my_res

            # Assign the radius
            if coarse:
                atom['rad'] = float(line[60:65])
                if atom['rad'] == 0:
                    atom['rad'] = 0.001
            else:
                atom['rad'] = get_radius(atom)
            # Add the atom to the atoms list
            atoms.append(atom)
        # If the line is not an atom line store the other data
        else:
            data.append(my_file[i].split())
    # Set the atoms and the data
    sys.atoms, sys.data = DataFrame(atoms), data
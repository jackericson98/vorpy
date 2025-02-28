import numpy as np
from pandas import DataFrame
from System.sys_objs.atom import make_atom
from System.sys_objs.chain import Sol, Chain
from System.sys_objs.residue import Residue
from System.chemistry_interpreter import residue_names, residue_atoms, my_masses
from System.sys_funcs.input.pdb import fix_sol


# Read cif function. Interprets the data in a cif file
def read_cif1(sys, file=None):
    # Check to see if the file is provided and use the bse file if not
    if file is None and sys.base_file[-3:] == 'cif':
        file = sys.base_file
    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Get the starting number for the line
    num = int(my_file[0].split()[0])
    # Go through each line of the file
    for i in range(len(my_file)):
        # Split the line
        line = my_file[i].split()
        # Add the atoms
        if line == int(num) and len(line) >= 7:
            sys.atoms.append(make_atom([line[9], line[10], line[11]], element=line[3], index=i))


def read_cif(sys, file):
    """
    Read crystallographic file format


    """
    # Check if there is a file input
    if file is None:
        file = sys.files['base_file']

    # Create the file dictionary
    file_dict = {'balls': [], 'Additional Information': []}

    # Open the file
    with open(file, 'r') as rf:
        # Set up the just in case occupancy doubled warning
        printed_occ_warn = False
        # Set up the atom counter
        atom_count, reset_checker = 0, 0
        chains, resids = {}, {}
        # Loop through the file
        for line in rf.readlines():
            # Split the line
            linfo = line.split()
            # Check if it is an atom line
            if linfo[0] == 'ATOM' or linfo[0] == 'HETATM':
                # Location
                loc = np.array([float(_) for _ in linfo[10:13]])
                # Get the occupancy assignment
                if linfo[4] != 'A' or linfo[4] != '.':
                    if not printed_occ_warn:
                        print("Warning! This molecule has multiple occupancy. Edit structure accordingly. "
                              "Program will default to occupancy \"A\"")
                        printed_occ_warn = True
                    continue

                # Create the ball
                ball = make_atom(location=loc, index=int(linfo[1]), element=linfo[2], name=linfo[3], occ_choice=linfo[4],
                                 res_name=linfo[5], chn_name=linfo[6], chn_id=int(linfo[7]), res_seq=int(linfo[8]),
                                 pdb_ins_code=linfo[9], occupancy=linfo[13], b_factor=linfo[14], charge=linfo[15],
                                 auth_seq_id=linfo[16], auth_comp_id=linfo[17], auth_asym_id=linfo[18],
                                 auth_atom_id=linfo[19], pdbx_PDB_model_num=linfo[20])
                res_str = linfo[5]
                res_seq = int(linfo[8])
                atom_count += 1
                if chain_str == ' ':
                    if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'}:
                        chain_str = 'SOL'
                    else:
                        chain_str = 'A'
                elif sys.type == 'foam' and res_str.lower() != 'bub' and chain_str != '0':
                    chain_str = 'SOL'

                # Create the chain and residue dictionaries
                res_name, chn_name = chain_str + '_' + line[17:20] + str(ball['res_seq']) + '_' + str(
                    reset_checker), chain_str
                # If the chain has been made before
                if chn_name in chains:
                    # Get the chain from the dictionary and add the atom
                    my_chn = chains[chn_name]
                    my_chn.add_atom(ball['num'])
                    ball['chn'] = my_chn
                # Create the chain
                else:
                    # If the chain is the sol chain
                    if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion',
                                           'cla'} or chn_name == 'SOL':
                        my_chn = Sol(atoms=[ball['num']], residues=[], name=chn_name, sys=sys)
                        sys.sol = my_chn
                    # If the chain is not sol create a regular chain object
                    else:
                        my_chn = Chain(atoms=[ball['num']], residues=[], name=chn_name, sys=sys)
                        sys.chains.append(my_chn)
                    # Set the chain in the dictionary and give the atom it's chain
                    chains[chn_name] = my_chn
                    ball['chn'] = my_chn

                # Assign the atoms and create the residues
                if res_name in resids:
                    my_res = resids[res_name]
                    my_res.atoms.append(ball['num'])
                else:
                    my_res = Residue(sys=sys, atoms=[ball['num']], name=res_str, sequence=ball['res_seq'],
                                     chain=ball['chn'])
                    resids[res_name] = my_res
                    if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion',
                                           'cla'} or chain_str == 'SOL':
                        sys.sol.residues.append(my_res)
                    else:
                        sys.residues.append(my_res)
                        ball['chn'].residues.append(my_res)
                # Assign the residue to the atom
                ball['res'] = my_res

                # Add the atom to the atoms list
                file_dict['balls'].append(ball)
                if res_seq == 9999:
                    reset_checker += 1
            # Otherwise add it to the date
            file_dict['Additional Information'].append(line)

        # Check that the sys.sol is not Noner
        if sys.sol is None:
            sys.sol = Sol(sys, [], [])
        # Set up the stuff
        for res in sys.residues:
            if res.name.lower() not in residue_names and res.chain.name != 'SOL':
                residue_names[res.name.lower()] = res.name.upper()
                residue_atoms[res.name.upper()] = {file_dict['balls'][_]['name'] for _ in res.atoms}

        # Set the atoms and the data
        sys.balls, sys.data = DataFrame(file_dict['ballls']), file_dict['Additional Information']
        # Adjust the SOL residues
        adjusted_residues = []
        for res in sys.sol.residues:
            if len(res.atoms) > 3:
                try:
                    adjusted_residues += fix_sol(sys, res)
                except TypeError:
                    print(res.atoms)
            else:
                adjusted_residues.append(res)

        sys.sol.residues = adjusted_residues


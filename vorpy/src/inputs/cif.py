import numpy as np
from pandas import DataFrame
from vorpy.src.objects import make_atom
from vorpy.src.objects import Sol, Chain
from vorpy.src.objects import Residue
from vorpy.src.chemistry import residue_names
from vorpy.src.chemistry import residue_atoms
from vorpy.src.inputs.fix_sol import fix_sol


def read_cif(sys, file):
    """
    Read and process a Crystallographic Information File (CIF) format into a system object.

    Parameters:
    -----------
    sys : System
        The system object to populate with CIF data
    file : str, optional
        Path to the CIF file. If None, uses sys.files['base_file']

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom objects from CIF data
        - Organizing atoms into residues and chains
        - Handling special cases like solvent molecules
        - Setting up system properties and metadata

    Notes:
    ------
    - Supports both ATOM and HETATM records
    - Handles multiple occupancy states (defaults to occupancy "A")
    - Special handling for solvent molecules (HOH, SOL, etc.)
    - Maintains residue and chain organization
    - Preserves atom metadata (B-factors, charges, etc.)
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

                # Create the ball object
                ball = make_atom(location=loc, index=int(linfo[1]), element=linfo[2], name=linfo[3], occ_choice=linfo[4],
                                 res_name=linfo[5], chn_name=linfo[6], chn_id=int(linfo[7]), res_seq=int(linfo[8]),
                                 pdb_ins_code=linfo[9], occupancy=linfo[13], b_factor=linfo[14], charge=linfo[15],
                                 auth_seq_id=linfo[16], auth_comp_id=linfo[17], auth_asym_id=linfo[18],
                                 auth_atom_id=linfo[19], pdbx_PDB_model_num=linfo[20])
                # Get the residue string
                res_str = linfo[5]
                # Get the residue sequence
                res_seq = int(linfo[8])
                # Increment the atom count
                atom_count += 1
                # If the chain string is empty
                if chain_str == ' ':
                    # If the residue string is a solvent
                    if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'}:
                        chain_str = 'SOL'
                    # Otherwise set the chain string to 'A'
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
                    # Get the residue from the dictionary and add the atom
                    my_res = resids[res_name]
                    my_res.atoms.append(ball['num'])
                else:
                    # Create the residue object
                    my_res = Residue(sys=sys, atoms=[ball['num']], name=res_str, sequence=ball['res_seq'],
                                     chain=ball['chn'])
                    # Add the residue to the dictionary
                    resids[res_name] = my_res
                    # If the residue string is a solvent
                    if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion',
                                           'cla'} or chain_str == 'SOL':
                        sys.sol.residues.append(my_res)
                    # Otherwise add the residue to the residues list
                    else:
                        sys.residues.append(my_res)
                        ball['chn'].residues.append(my_res)
                # Assign the residue to the atom
                ball['res'] = my_res

                # Add the atom to the atoms list
                file_dict['balls'].append(ball)
                # If the residue sequence is 9999
                if res_seq == 9999:
                    # Increment the reset checker
                    reset_checker += 1
            # Otherwise add it to the date
            file_dict['Additional Information'].append(line)

        # Check that the sys.sol is not None
        if sys.sol is None:
            # Create a new sol object
            sys.sol = Sol(sys, [], [])
        # Set up the residue names and atoms
        for res in sys.residues:
            # If the residue name is not in the residue names dictionary
            if res.name.lower() not in residue_names and res.chain.name != 'SOL':
                # Add the residue name to the dictionary
                residue_names[res.name.lower()] = res.name.upper()
                # Add the residue atoms to the dictionary
                residue_atoms[res.name.upper()] = {file_dict['balls'][_]['name'] for _ in res.atoms}

        # Set the atoms and the data
        sys.balls, sys.data = DataFrame(file_dict['ballls']), file_dict['Additional Information']
        # Adjust the SOL residues
        adjusted_residues = []
        for res in sys.sol.residues:
            # If the residue has more than 3 atoms
            if len(res.atoms) > 3:
                # Try to adjust the SOL residues
                try:
                    adjusted_residues += fix_sol(sys, res)
                except TypeError:
                    print(res.atoms)
            # Otherwise add the residue to the adjusted residues list
            else:
                adjusted_residues.append(res)
        # Set the sol residues
        sys.sol.residues = adjusted_residues

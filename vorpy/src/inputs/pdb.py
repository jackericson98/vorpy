import os
from vorpy.src.objects import make_atom
from vorpy.src.objects import Residue
from vorpy.src.objects import Chain, Sol
from vorpy.src.chemistry import residue_names
from vorpy.src.chemistry import residue_atoms
from vorpy.src.chemistry import my_masses
from vorpy.src.inputs.fix_sol import fix_sol
import os.path as path
import numpy as np
from pandas import DataFrame


def read_pdb(sys, file=None):
    """
    Read and process a PDB format file into a system object.

    This function parses PDB format files, which contain molecular structure data including:
    - Atom coordinates and properties
    - Chain and residue information
    - Secondary structure elements
    - Connectivity information

    Parameters:
    -----------
    sys : System
        The system object to populate with PDB data
    file : str, optional
        Path to the PDB file. If None, uses sys.files['base_file'] if it ends in .pdb

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom objects from PDB data
        - Organizing atoms into chains and residues
        - Storing atom properties and coordinates
        - Handling special cases like foam and coarse-grained systems
    """
    # Check to see if the file is provided and use the base file if not
    if file is None and sys.files['base_file'][-3:] == 'pdb':
        file = sys.files['base_file']
    # Check if the file exists and if it does use the vpy_dir
    if path.exists(file) and file[0] == '.' and sys.files['vpy_dir'] is not None:
        file_address = sys.files['vpy_dir'] + file[1:]
    # If the file exists and is not a relative path
    elif path.exists(file):
        file_address = file
    # If the file exists and the vpy_dir is not None
    elif sys.files['vpy_dir'] is not None and path.exists(sys.files['vpy_dir'] + file):
        file_address = sys.files['vpy_dir'] + file
    # If the file exists and the dir is not None
    elif sys.files['dir'] is not None and path.exists(sys.files['dir'] + file):
        file_address = sys.files['dir'] + file
    # If the file exists and the dir is not None and the file is a relative path
    elif sys.files['dir'] is not None and path.exists(sys.files['dir'] + file[1:]):
        file_address = sys.files['dir'] + file[1:]
    # If the file does not exist return
    else:
        return
    # Print a statement saying the file is being read
    print(f"\rReading File {os.path.basename(file)}", end="")
    # Get the file information and make sure to close the file when done
    with open(file_address, 'r') as f:
        my_file = f.readlines()
    # If the file is empty
    if len(my_file) == 0:
        # Open the file in the absolute path
        with open(os.path.abspath(file_address), 'r') as f:
            # Read the file
            my_file = f.readlines()

    # Add the system name and reset the atoms and data lists
    sys.name = path.basename(sys.files['base_file'])[:-4]
    # Set up the atom and the data lists
    atoms, data, atom_count, reset_checker = [], [], 0, 0
    previous_res_seq, previous_chain = None, None
    # Initialize the chains and residues lists
    sys.chains, sys.residues = [], []
    # Initialize the chains and residues dictionaries
    chains, resids = {}, {}
    # Check if the file is a foam file
    if my_file[0].split()[1] == 'foam_gen':
        # Set the system type to foam
        sys.type = 'foam'
        try:
            # Get the box width
            bw = float(my_file[0].split()[2])
            # Set the foam box
            sys.foam_box = [[0, 0, 0], [bw, bw, bw]]
            # Set the foam data
            sys.foam_data = my_file[0].split()[2:]
        except ValueError:
            # Get the box width
            bw = float(my_file[0].split()[5][:-1])
            # Set the foam box
            sys.foam_box = [[0, 0, 0], [bw, bw, bw]]
            # Set the foam data
            sys.foam_data = my_file[0].split()[5:]
    # Check if the file is a coarse file
    if my_file[0].split()[1] == 'coarsify':
        # Set the system type to coarse
        sys.type = 'coarse'
    # Go through each line in the file and check if the first word is the word we are looking for
    for i in range(len(my_file)):
        # Print a loading statement for longer files
        print(f"\rReading File {os.path.basename(file)} {round(100 * (i/len(my_file)), 2)} %", end="")
        # Check to make sure the line isn't empty
        if len(my_file[i]) == 0:
            continue
        # Pull the file line and first word
        line = my_file[i]
        word = line[:6].lower().strip()
        # Check to see if the line is an atom line
        if line and word in {'atom', 'hetatm'}:
            # Check for the "m" situation
            if line[76:78] == ' M':
                continue
            # Get the name
            name = line[12:16]
            # Get the residue sequence
            res_seq = line[22:26]
            if res_seq == '    ':
                res_seq = 0
            else:
                res_seq = int(res_seq)
            # If no chain is specified, set the chain to 'None'
            try:
                res_str, chain_str = line[17:20].strip(), line[21]
            except IndexError:
                continue
            # Assign the radius
            rad = None
            # If the system is a foam or coarse system
            if sys.type == 'foam' or sys.type == 'coarse':
                # Get the radius
                rad = float(line[60:66])
                # If the radius is 0
                if rad == 0:
                    # Set the radius to 0.001
                    rad = 0.001
            # Get the mass for the atom:
            if sys.type == 'mol' and line[76:78].strip().lower() in my_masses:
                # Get the mass
                mass = my_masses[line[76:78].strip().lower()]
            elif sys.type == 'foam':
                # Get the mass
                mass = (4 / 3) * np.pi * rad ** 3
            elif sys.type == 'coarse':
                # Get the mass
                mass = float(line[54:60])
            else:
                # Set the mass to 1
                mass = 1

            # Create the atom
            atom = make_atom(location=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                             system=sys, element=line[76:78].strip(), res_seq=res_seq, res_name=res_str,
                             chn_name=chain_str,
                             name=name.strip(), seg_id=line[72:76], index=atom_count, mass=mass, radius=rad)
            # Increment the atom count
            atom_count += 1
            # If the chain is empty
            if chain_str == ' ':
                # If the residue is a sol or hoh
                if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'}:
                    # Set the chain to SOL
                    chain_str = 'SOL'
                # Otherwise set the chain to A
                else:
                    chain_str = 'A'
            # If the system is a foam and the residue is not bub and the chain is not 0
            elif sys.type == 'foam' and res_str.lower() != 'bub' and chain_str != '0':
                # Set the chain to SOL
                chain_str = 'SOL'
            # Detect PDB residue-number rollover within the same chain.
            # PDB residue sequence numbers are limited to four columns, so very large
            # residue counts may eventually restart at a smaller sequence number.
            if previous_chain == chain_str and previous_res_seq is not None and res_seq < previous_res_seq:
                reset_checker += 1

            previous_chain, previous_res_seq = chain_str, res_seq

            # Create the chain and residue dictionaries
            res_name, chn_name = chain_str + '_' + line[17:20] + str(atom['res_seq']) + '_' + str(reset_checker), chain_str
            # If the chain has been made before
            if chn_name in chains:
                # Get the chain from the dictionary and add the atom
                my_chn = chains[chn_name]
                # Add the atom to the chain
                my_chn.add_atom(atom['num'])
                # Set the chain for the atom
                atom['chn'] = my_chn
            # Create the chain
            else:
                # If the chain is the sol chain
                if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'} or chn_name == 'SOL':
                    # Create the sol chain
                    my_chn = Sol(atoms=[atom['num']], residues=[], name=chn_name, sys=sys)
                    # Set the sol chain
                    sys.sol = my_chn
                # If the chain is not sol create a regular chain object
                else:
                    # Create the chain
                    my_chn = Chain(atoms=[atom['num']], residues=[], name=chn_name, sys=sys)
                    # Add the chain to the chains list
                    sys.chains.append(my_chn)
                # Set the chain in the dictionary and give the atom it's chain
                chains[chn_name] = my_chn
                # Set the chain for the atom
                atom['chn'] = my_chn

            # Assign the atoms and create the residues
            if res_name in resids:
                # Get the residue from the dictionary and add the atom
                my_res = resids[res_name]
                # Add the atom to the residue
                my_res.atoms.append(atom['num'])
            # If the residue is not in the dictionary
            else:
                # Create the residue
                my_res = Residue(sys=sys, atoms=[atom['num']], name=res_str, sequence=atom['res_seq'],
                                 chain=atom['chn'])
                # Add the residue to the dictionary
                resids[res_name] = my_res
                # If the residue is a sol or hoh or the chain is SOL
                if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'} or chain_str == 'SOL':
                    # Add the residue to the sol residues
                    sys.sol.residues.append(my_res)
                else:
                    # Add the residue to the residues list
                    sys.residues.append(my_res)
                    # Add the residue to the chain residues
                    atom['chn'].residues.append(my_res)
            # Assign the residue to the atom
            atom['res'] = my_res

            # Add the atom to the atoms list
            atoms.append(atom)
        # If the line is not an atom line store the other data
        else:
            data.append(my_file[i].split())
    # Check that the sys.sol is not Noner
    if sys.sol is None:
        sys.sol = Sol(sys, [], [])
    # Set up the stuff
    for res in sys.residues:
        # If the residue name is not in the residue names dictionary and the chain is not SOL
        if res.name.lower() not in residue_names and res.chain.name != 'SOL':
            # Add the residue name to the residue names dictionary
            residue_names[res.name.lower()] = res.name.upper()
            # Add the residue atoms to the residue atoms dictionary
            residue_atoms[res.name.upper()] = {atoms[_]['name'] for _ in res.atoms}


    # Set the atoms and the data
    sys.balls, sys.data = DataFrame(atoms), data
    # Adjust the SOL residues
    adjusted_residues = []
    # Go through the sol residues
    for res in sys.sol.residues:
        # If the residue has more than 3 atoms
        if len(res.atoms) > 3:
            try:
                # Fix the sol
                adjusted_residues += fix_sol(sys, res)
            except TypeError:
                # Print the residue atoms
                print(res.atoms)
        else:
            adjusted_residues.append(res)
    # Add the sol residues
    sys.sol.residues = adjusted_residues
    # Create a file read statement for the user
    print(f"\r{os.path.basename(file)} successfully added.", end="")
    return sys


def read_pdb_simple(file):
    """
    Read a PDB file and return a list of atoms.
    """
    # Read the file
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Get the atoms
    atoms = {}
    # Go through the file
    for line in my_file:
        # If the line is an atom line
        if line[0:6].strip() in ['ATOM', 'HETATM']:
            atom = read_pdb_line(line)
            atoms[atom['atom_serial_number']] = atom
            atoms[atom['atom_serial_number']]['radius'] = atom['temperature_factor']
    # Return the atoms
    return atoms


def read_pdb_line(pdb_line):
    """
    Parse a PDB format line into a dictionary of atom information.

    This function processes a single line from a PDB file and extracts all relevant
    atom information according to the PDB format specification. The function handles
    the fixed-width format of PDB files and converts string values to appropriate
    data types.

    Parameters:
    -----------
    pdb_line : str
        A single line from a PDB file containing atom information

    Returns:
    --------
    dict
        A dictionary containing the following keys:
        - record_type: Type of record (e.g., ATOM, HETATM)
        - atom_serial_number: Unique identifier for the atom
        - atom_name: Name of the atom
        - alternate_location: Alternate location indicator
        - residue_name: Name of the residue
        - chain_identifier: Chain identifier
        - residue_sequence_number: Sequence number of the residue
        - insertion_code: Insertion code for the residue
        - x_coordinate: X coordinate of the atom
        - y_coordinate: Y coordinate of the atom
        - z_coordinate: Z coordinate of the atom
        - occupancy: Occupancy value
        - temperature_factor: Temperature factor
        - segment_identifier: Segment identifier
        - element_symbol: Element symbol
        - charge: Charge on the atom

    Notes:
    ------
    - Follows the standard PDB format specification
    - Handles fixed-width fields with appropriate stripping
    - Converts numeric values to appropriate types
    - Preserves empty fields as empty strings
    """
    # Return the dictionary
    return {
        "record_type": pdb_line[0:6].strip(),
        "atom_serial_number": int(pdb_line[6:11].strip()),
        "atom_name": pdb_line[12:16].strip(),
        "alternate_location": pdb_line[16].strip(),
        "residue_name": pdb_line[17:20].strip(),
        "chain_identifier": pdb_line[21].strip(),
        "residue_sequence_number": int(pdb_line[22:26].strip()),
        "insertion_code": pdb_line[26].strip(),
        "x_coordinate": float(pdb_line[30:38].strip()),
        "y_coordinate": float(pdb_line[38:46].strip()),
        "z_coordinate": float(pdb_line[46:54].strip()),
        "occupancy": float(pdb_line[54:60].strip()),
        "temperature_factor": float(pdb_line[60:66].strip()),
        "segment_identifier": pdb_line[72:76].strip(),
        "element_symbol": pdb_line[76:78].strip(),
        "charge": pdb_line[78:80].strip()
    }

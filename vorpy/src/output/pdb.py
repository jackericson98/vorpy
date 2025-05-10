import os
import shutil
from shutil import SameFileError


def make_pdb_line(atom="ATOM", ser_num=0, name="", alt_loc=" ", res_name="", chain="A", res_seq=0, cfir="", x=0, y=0, z=0,
                  occ=1, tfact=0, seg_id="", elem="", charge=""):
    """
    Formats atom data into a properly formatted PDB file line.

    This function takes individual atom properties and formats them according to the PDB file format
    specification, ensuring proper spacing and alignment of all fields. The output string follows
    the standard PDB format with fixed-width columns for each field.

    Args:
        atom (str): Record type (default: "ATOM")
        ser_num (int): Atom serial number
        name (str): Atom name
        alt_loc (str): Alternate location indicator
        res_name (str): Residue name
        chain (str): Chain identifier
        res_seq (int): Residue sequence number
        cfir (str): Code for insertion of residues
        x (float): X coordinate
        y (float): Y coordinate
        z (float): Z coordinate
        occ (float): Occupancy
        tfact (float): Temperature factor
        seg_id (str): Segment identifier
        elem (str): Element symbol
        charge (str): Charge on the atom

    Returns:
        str: A properly formatted PDB file line string with all fields aligned according to PDB specifications
    """
    # Write the line for the file
    return "{:<6}{:>5} {:<4}{:1}{:>3} {:^1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4}{:>2}{}\n"\
        .format(atom, ser_num, name, alt_loc, res_name, chain[0], res_seq, cfir, x, y, z, occ, tfact, seg_id, elem, charge)


def write_pdb(atoms, file_name, sys, directory=None):
    """
    Writes a PDB (Protein Data Bank) file containing the specified atoms.

    This function creates a PDB file either by copying an existing base file (if available)
    or by manually constructing the PDB format from atom data. The output file will contain
    only the specified atoms while maintaining proper PDB formatting.

    Parameters:
        atoms (list): List of atoms to include in the PDB file. Can be either:
            - List of atom objects
            - List of integer indices corresponding to atoms in sys.balls
        file_name (str): Name of the output PDB file (without .pdb extension)
        sys (System): System object containing the full atom set and base PDB file reference
        directory (str, optional): Directory path where the PDB file should be written.
            If None, uses the current working directory.

    Returns:
        None: Writes a PDB file to the specified location

    Notes:
        - If sys.files['base_file'] exists and contains all atoms, the function will
          copy the relevant lines from the base file for efficiency
        - If no base file exists, the function will construct the PDB file manually
        - The output file will include a custom header with system name and group info
        - Empty atom lists will result in no file being created

    Examples:
        # Write all atoms from a system to a PDB file
        >>> atoms = list(range(len(sys.balls)))
        >>> write_pdb(atoms, "full_structure", sys)

        # Write specific atoms to a PDB file in a custom directory
        >>> selected_atoms = [0, 5, 10]  # Indices of atoms to include
        >>> write_pdb(selected_atoms, "subset", sys, directory="/path/to/output")
    """
    # Catch empty atoms cases
    if atoms is None or len(atoms) == 0:
        return
    # Make note of the starting directory
    start_dir = os.getcwd()
    # Change to the specified directory
    if directory is not None:
        os.chdir(directory)

    # Check to see if a system was provided
    if sys.files['base_file'] is not None:

        # If the output is all atoms just copy the pdb
        if len(atoms) == len(sys.balls):
            try:
                shutil.copy(sys.files['base_file'], os.getcwd() + '/' + file_name + '.pdb')
            except SameFileError:
                pass
            return

        # Open the file for writing
        with open(file_name + ".pdb", 'w') as pdb_file:

            # Open the base file and read the lines
            with open(sys.files['base_file'], 'r') as f:
                read_file = f.readlines()

            # Write a header for the pdb
            pdb_file.write("HEADER  vorpy output - " + sys.name + " group " + file_name + " atoms\n")
            # Figure out what lines the atoms start on
            offset = 0
            while read_file[offset][:6].lower().strip() not in {'atom', 'hetatm'}:
                offset += 1

            # Grab the lines from the initial pdb
            for j in range(len(sys.balls)):
                if sys.balls.iloc[j]['num'] in atoms:
                    pdb_file.write(read_file[j + offset])

    # Manually write the pdb file
    else:
        # Open the file for writing
        with open(file_name + ".pdb", 'w') as pdb_file:
            # Go through each atom in the system
            for i, a in enumerate(atoms):
                # Get the ball
                if type(a) is int:
                    a = sys.balls.iloc[a]
                # Get the location string
                x, y, z = a['loc']
                # Get the information from the atom in writable format
                tfact = 0
                if sys.type == 'foam' or sys.type == 'coarse':
                    tfact = a['rad']
                # Write the atom information
                pdb_file.write(make_pdb_line(ser_num=a['num'], name=a['name'], res_name=a['res'].name, chain=a['chn'].name,
                                             res_seq=a['res_seq'], x=x, y=y, z=z, tfact=tfact, elem=a['element']))
    # Change back to the starting directory
    os.chdir(start_dir)

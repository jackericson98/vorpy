import os
import shutil
import numpy as np
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
    if atoms is None or len(atoms) == 0:
        return

    start_dir = os.getcwd()

    try:
        if directory is not None:
            os.chdir(directory)

        # --------------------------------------------------------------
        # Existing base PDB
        # --------------------------------------------------------------

        if sys.files['base_file'] is not None:

            # Full system: simply copy original file.
            if len(atoms) == len(sys.balls):
                try:
                    shutil.copy(sys.files['base_file'], os.path.join(os.getcwd(), file_name + '.pdb'))
                except SameFileError:
                    pass
                except OSError:
                    pass
                return

            # ----------------------------------------------------------
            # Determine selected atom numbers
            # ----------------------------------------------------------

            first = atoms[0]

            if isinstance(first, (int, np.integer)):
                # Group atom lists normally contain indices into sys.balls.
                selected_nums = set(sys.balls.iloc[list(atoms)]['num'].tolist())
            else:
                # Atom/Series objects.
                selected_nums = {a['num'] for a in atoms}

            # ----------------------------------------------------------
            # Read base PDB
            # ----------------------------------------------------------

            with open(sys.files['base_file'], 'r') as f:
                read_file = f.readlines()

            offset = 0

            while offset < len(read_file) and read_file[offset][:6].lower().strip() not in {'atom', 'hetatm'}:
                offset += 1

            # ----------------------------------------------------------
            # Write selected atoms
            # ----------------------------------------------------------

            output_lines = [
                "HEADER  vorpy output - " + sys.name + " group " + file_name + " atoms\n"
            ]

            nums = sys.balls['num'].to_numpy()

            for j, num in enumerate(nums):
                if num in selected_nums:
                    output_lines.append(read_file[j + offset])

            with open(file_name + ".pdb", 'w', buffering=1024 * 1024) as pdb_file:
                pdb_file.writelines(output_lines)

        # --------------------------------------------------------------
        # No base PDB: construct manually
        # --------------------------------------------------------------

        else:
            output_lines = []

            for a in atoms:
                if isinstance(a, (int, np.integer)):
                    a = sys.balls.iloc[a]

                x, y, z = a['loc']
                tfact = a['rad'] if sys.type in {'foam', 'coarse'} else 0

                output_lines.append(
                    make_pdb_line(
                        ser_num=a['num'],
                        name=a['name'],
                        res_name=a['res'].name,
                        chain=a['chn'].name,
                        res_seq=a['res_seq'],
                        x=x,
                        y=y,
                        z=z,
                        tfact=tfact,
                        elem=a['element']
                    )
                )

            with open(file_name + ".pdb", 'w', buffering=1024 * 1024) as pdb_file:
                pdb_file.writelines(output_lines)

    finally:
        os.chdir(start_dir)

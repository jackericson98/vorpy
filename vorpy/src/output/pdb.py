import os
from numbers import Integral
from pathlib import Path
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
    # Field widths in str.format are minima, so reject serial overflow rather
    # than shifting every subsequent PDB column.
    if not isinstance(ser_num, Integral) or not 0 <= ser_num <= 99999:
        raise ValueError("PDB atom serial must be an integer between 0 and 99999.")
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
        - Without a source, molecular atoms require preserved original PDB serial metadata.
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
    from vorpy.src.inputs.frames import first_pdb_frame, is_pdb_virtual_site

    # Catch empty atoms cases
    if atoms is None or len(atoms) == 0:
        return
    destination = Path(directory or os.getcwd()) / (file_name + '.pdb')

    # Check to see if a system was provided
    if sys.files['base_file'] is not None:

        # If the output is all atoms just copy the pdb
        if len(atoms) == len(sys.balls) and getattr(sys, 'frame_count', 1) <= 1:
            try:
                shutil.copy(sys.files['base_file'], destination)
            except SameFileError:
                pass
            return

        # Open the file for writing
        with destination.open('w') as pdb_file:

            # Open the base file and read the lines.
            with open(sys.files['base_file'], 'r') as f:
                read_file = list(first_pdb_frame(f))

            # Write a header for the pdb.
            pdb_file.write(
                "HEADER  vorpy output - "
                + sys.name
                + " group "
                + file_name
                + " atoms\n"
            )

            # IMPORTANT:
            # System atom indices correspond to the ordered ATOM/HETATM
            # records, NOT to raw line positions in the PDB file.
            #
            # A PDB may contain TER, MODEL, ENDMDL, REMARK, ANISOU, etc.
            # between atom records.  The old ``read_file[j + offset]``
            # assumption therefore shifted exports after the first
            # interspersed non-atom record.  1BNA contains a TER before
            # chain B, which made system index 736 (N9) export the
            # preceding H1' record.
            atom_records = [
                line
                for line in read_file
                if line[:6].strip().upper() in {"ATOM", "HETATM"}
                and not is_pdb_virtual_site(line)
            ]

            selected = []
            for atom_index in atoms:
                try:
                    atom_index = int(atom_index)
                except (TypeError, ValueError):
                    continue

                if 0 <= atom_index < len(atom_records):
                    selected.append(atom_index)

            for atom_index in sorted(set(selected)):
                pdb_file.write(atom_records[atom_index])

    # Manually write the pdb file
    else:
        selected_atoms = [sys.balls.iloc[int(atom)] if isinstance(atom, Integral) else atom for atom in atoms]
        serials = []
        for atom in selected_atoms:
            serial = atom.get('pdb_serial')
            if not isinstance(serial, Integral):
                if getattr(sys, 'type', None) == 'mol':
                    raise ValueError(
                        'Original PDB serial identities are unavailable. Attach a verified source structure '
                        'before exporting molecular PDBs; internal atom indices are not original serials.')
                serial = atom['num']
            serials.append(serial)
        # Open the file for writing
        with destination.open('w') as pdb_file:
            # Go through each atom in the system
            for a, serial in zip(selected_atoms, serials):
                # Get the location string
                x, y, z = a['loc']
                # Get the information from the atom in writable format
                tfact = 0
                if sys.type == 'foam' or sys.type == 'coarse':
                    tfact = a['rad']
                # Write the atom information
                chain = a.get('chain_name')
                chain = (str(chain) or ' ') if chain is not None else a['chn'].name
                pdb_file.write(make_pdb_line(ser_num=serial, name=a['name'], res_name=a['res'].name, chain=chain,
                                             res_seq=a['res_seq'], x=x, y=y, z=z, tfact=tfact, elem=a['element']))

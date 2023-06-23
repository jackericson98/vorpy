import os
import shutil
from System.sys_funcs.output.surfs import write_surfs
from System.sys_funcs.output.edges import write_edges
from System.sys_funcs.output.verts import write_verts


def write_pdb(atoms, file_name, sys, directory=None):
    """
    Creates a pdb file type in the current working directory
    :param atoms: List of atom type objects for writing
    :param file_name: Name of the output file
    :param sys: System object used for writing the whole pbd file
    :param directory: Output directory for the file
    :return: Writes a pdb file for the set of atoms
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
    if sys.base_file is not None:

        # If the output is all atoms just copy the pdb
        if len(atoms) == len(sys.atoms):
            shutil.copy(sys.base_file, os.getcwd() + '/' + file_name + '.pdb')
            return

        # Open the file for writing
        with open(file_name + ".pdb", 'w') as pdb_file:

            # Open the base file and read the lines
            with open(sys.base_file, 'r') as f:
                read_file = f.readlines()

            # Write a header for the pdb
            pdb_file.write("HEADER  vorpy output - " + sys.name + " group " + file_name + " atoms\n")
            # Figure out what lines the atoms start on
            offset = 0
            while read_file[offset][:4].lower() != 'atom':
                offset += 1

            # Grab the lines from the initial pdb
            for j in range(len(sys.atoms)):
                if sys.atoms.iloc[j]['num'] in atoms:
                    pdb_file.write(read_file[j + offset])

    # Manually write the pdb file
    else:
        # Open the file for writing
        with open(file_name + ".pdb", 'w') as pdb_file:
            # Go through each atom in the system
            for i in atoms:
                if type(i) == int:
                    a = sys.net.atoms.iloc[i]
                else:
                    a = i
                    i = a['num']
                # Get the location string
                loc = ["{:.3f}".format(_) for _ in a['loc']]
                # Get the information from the atom in writable format
                ser_num = " " * (5 - len(str(i+1))) + str(i + 1)
                file_name = a['name'] + " " * (4 - len(a['name']))
                if 'residue' in a:
                    res = " " * (3 - len(a['residue'])) + a['residue']
                else:
                    res = "   "
                if 'chn' not in a or a['chn'].name.lower() == "zz" or a['chn'].name.lower() == 'mol' or a['chn'].name.lower() == 'sol':
                    chain = " "
                else:
                    chain = str(a.chn.name)
                if 'res_seq' in a:
                    res_seq = " " * (3 - len(str(a['res_seq']))) + str(a['res_seq'])
                else:
                    res_seq = "   "
                loc_strs = [" " * (7 - len(_)) + _ for _ in loc]
                occupancy = " " * 5
                t_fact = " " * 5
                if 'seg_id' in a:
                    seg_id = a['seg_id'] + " " * (3 - len(a['seg_id']))
                else:
                    seg_id = "   "
                if 'element' in a:
                    symbol = a['element']
                else:
                    symbol = 'h'
                charge = ''
                # Write the atom information
                pdb_file.write("ATOM  " + ser_num + " " + file_name + " " + res + " " + chain + res_seq + "    " +
                               " ".join(loc_strs) + occupancy + t_fact + "      " + seg_id + symbol + charge + "\n")
    # Change back to the starting directory
    os.chdir(start_dir)


def write_gro(atoms, file_name, sys=None, directory=None):
    """
    Writes a gro file for the atoms specified
    :param atoms: Atoms for writing
    :param file_name: Name of the output file
    :param sys: System to pull from
    :param directory: Output directory for the file
    :return: Outputs the file
    """
    # Change to the directory of specified
    if directory is not None and os.path.exists(directory):
        os.chdir(directory)

    # Create the title
    sys_name = sys.name

    # Copy the
    # Check to see if a system was provided
    if sys is not None and sys.base_file is not None:

        # If the output is all atoms just copy the pdb
        if len(atoms) == len(sys.atoms):
            shutil.copy(sys.base_file, os.getcwd() + file_name)
            return

    # Open the file
    with open(file_name + '.gro', 'w') as f:

        # Write the header
        f.write("{}\n{:5d}\n".format(sys_name, len(atoms)))
        # Write the atoms information
        for atom in atoms:
            f.write("{:5d}{:5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n"
                    .format(atom['res_seq'], atom['res'].name, atom['name'], atom['num'] + 1, *atom['loc']))
        # Write the box
        box = sys.net.box
        f.write("{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}\n".format(*box[0], *box[1]))


def write_atom_cells(net, atoms, directory=None, surfs=True, edges=False, verts=False):
    """
    Writes individual cell files for each of the atoms specified
    :param atoms: Atom objects for outputting
    :param directory: Output directory
    :param surfs: Bool for exporting the surfaces for the atoms
    :param edges: Bool for exporting the edges or not
    :param verts: Bool for exporting the vertices or not
    :return: None
    """
    # Change to the directory
    if directory is not None:
        os.chdir(directory)
    # Go through the atoms
    for i in atoms:
        atom = net.atoms.iloc[i]
        # Check if the surfaces should be exported
        if surfs:
            write_surfs(net, atom['asurfs'], directory=directory, file_name=str(atom['num']) + "_" + atom['name'])
        # Check for verts
        if verts:
            write_verts(net, atom['averts'], directory=directory, file_name=str(atom['num']) + "_" + atom['name'] + "_verts")
        # Check for edges
        if edges:
            write_edges(net, atom['aedges'], directory=directory, file_name=str(atom['num']) + "_" + atom['name'] + "_edges")

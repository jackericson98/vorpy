import os
import shutil
from datetime import datetime
from shutil import SameFileError

from System.sys_funcs.output.surfs import write_surfs
from System.sys_funcs.output.edges import write_edges
from System.sys_funcs.output.verts import write_off_verts


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
    if sys is not None and sys.files['base_file'] is not None:

        # If the output is all atoms just copy the pdb
        if len(atoms) == len(sys.atoms):
            shutil.copy(sys.files['base_file'], os.getcwd() + file_name)
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
        box = sys.net.box['verts']
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
        atom = net.balls.iloc[i]
        if not atom['complete']:
            continue
        # Check if the surfaces should be exported
        if surfs:
            write_surfs(net, atom['surfs'], directory=directory,
                        file_name='ball' + "_" + atom['name'].strip() + '_' + net.settings['net_type'],
                        color=(255, 0, 0) if net.settings['net_type'] == 'pow' else False)
        # Check for verts
        if verts:
            write_off_verts(net, atom['verts'], directory=directory,
                            file_name='ball_{}'.format(atom['name'].strip()) + "_" + net.settings['net_type'] + "_verts")
        # Check for edges
        if edges:
            write_edges(net, atom['edges'], directory=directory,
                        file_name='ball_{}'.format(atom['name'].strip()) + "_" + net.settings['net_type'] + "_edges")


def make_pdb_line(atom="ATOM", ser_num=0, name="", alt_loc=" ", res_name="", chain="A", res_seq=0, cfir="", x=0, y=0, z=0,
                  occ=1, tfact=0, seg_id="", elem="", charge=""):
    """
    Takes in values for a line in a pdb file and places them in the correct locations
    :return: String for each line
    """
    # Write the line for the file
    return "{:<6}{:>5} {:<4}{:1}{:>3} {:^1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4}{:>2}{}\n"\
        .format(atom, ser_num, name, alt_loc, res_name, chain[0], res_seq, cfir, x, y, z, occ, tfact, seg_id, elem, charge)


def write_atom_radii(my_sys, directory=None, file_name=None):
    # Check if a directory has been identified
    if directory is None:
        directory = my_sys.files['dir']
    # Check if the file_name has been specified
    if file_name is None:
        file_name = my_sys.name + '_atom_radii'
    # Open the file
    with open(directory + '/' + file_name + '.txt', 'w') as radii_file:
        # Write the header
        radii_file.write('{} solved at: {}\n\n'.format(my_sys.name, datetime.now()))
        # Write the elements header
        radii_file.write('Default Element Radii\n')
        # Loop through the elements
        for element in my_sys.element_radii:
            # Write the name of the element and the
            radii_file.write('{} = {} \u212B\n'.format(element, my_sys.element_radii[element]))
        # Write the special radii header
        radii_file.write('\nResidue Specific Radii\n')
        # Loop through the special radii
        for residue in my_sys.special_radii:
            for name in my_sys.special_radii[residue]:
                radii_file.write('{} {} = {} \u212B\n'.format(residue, name, my_sys.special_radii[residue][name]))

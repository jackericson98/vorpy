import os
from datetime import datetime
from vorpy.src.output.surfs import write_surfs
from vorpy.src.output.edges import write_edges
from vorpy.src.output.verts import write_off_verts


def write_atom_cells(net, atoms, directory=None, surfs=True, edges=False, verts=False, concave_colors=False):
    """
    Exports individual cell data files for specified atoms in a network.

    This function generates separate output files for each atom's cell components (surfaces, edges, and vertices)
    based on the network's Voronoi decomposition. Each atom's data can be exported as multiple file types
    depending on the specified parameters.

    Args:
        net: Network object containing the Voronoi decomposition data
        atoms: List of atom indices to process
        directory: Optional output directory path. If None, uses current directory
        surfs: If True, exports surface data for each atom (default: True)
        edges: If True, exports edge data for each atom (default: False)
        verts: If True, exports vertex data for each atom (default: False)
        concave_colors: If True, exports the concave colors for the surfaces. Default is False
    Returns:
        None: Creates individual files for each atom's cell components in the specified directory
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
                        color=(255, 0, 0) if net.settings['net_type'] == 'pow' else False,
                        concave_colors=concave_colors, ref_surfs=[i], universal_max=False)
        # Check for verts
        if verts:
            write_off_verts(net, atom['verts'], directory=directory,
                            file_name='ball_{}'.format(atom['name'].strip()) + "_" + net.settings['net_type'] + "_verts")
        # Check for edges
        if edges:
            write_edges(net, atom['edges'], directory=directory,
                        file_name='ball_{}'.format(atom['name'].strip()) + "_" + net.settings['net_type'] + "_edges")


def write_atom_radii(my_sys, directory=None, file_name=None):
    """
    Exports atom radii information to a text file.

    This function generates a text file containing detailed information about atom radii used in the system,
    including both default element radii and residue-specific radii. The output file provides a clear record
    of the radius values used for different elements and specific residues in the system.

    Args:
        my_sys: System object containing atom and radius information
        directory: Optional output directory path. If None, uses the system's default directory
        file_name: Optional name for the output file. If None, uses system name with '_atom_radii' suffix

    Returns:
        None: Creates a text file containing radius information in the specified directory

    Notes:
        - The output file includes a timestamp of when the radii were solved
        - Radii are written in Angstroms (Å)
        - The file is organized into two sections:
          1. Default Element Radii: Standard radii for each element
          2. Residue Specific Radii: Custom radii for specific residues
    """
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

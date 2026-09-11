import os
from datetime import datetime
from vorpy.src.output.surfs import write_surfs
from vorpy.src.output.edges import write_edges
from vorpy.src.output.verts import write_off_verts


def write_atom_cells(net, atoms, directory=None, surfs=True, edges=False, verts=False,
                     concave_colors=False, file_type='off', color_scheme=None,
                     color_map=None, color_limit=None):
    """Export each complete atom cell using cell-relative curvature colors.

    Individual atom cells have a well-defined outward orientation, so their
    surface/edge/vertex colors use that atom's own curvature perspective. A
    caller-supplied ``color_limit`` keeps all atom cells on the same scale as
    the corresponding group visualization.
    """
    if directory is not None:
        os.chdir(directory)

    scheme = net.settings.get('surf_scheme') if color_scheme is None else color_scheme
    cmap = net.settings.get('surf_col', 'coolwarm') if color_map is None else color_map

    for i in atoms:
        atom = net.balls.iloc[i]
        if not atom['complete']:
            continue

        cell_targets = [int(atom.get('num', i))]
        base_name = 'ball_' + atom['name'].strip() + '_' + net.settings['net_type']

        if surfs:
            write_surfs(
                net, atom['surfs'], directory=directory, file_name=base_name,
                color=(255, 0, 0) if net.settings['net_type'] == 'pow' else False,
                concave_colors=concave_colors, ref_surfs=cell_targets,
                universal_max=False, file_type=file_type,
                color_scheme=scheme, color_map=cmap, color_limit=color_limit,
                target_cells=cell_targets, color_mode='cell',
            )

        if verts:
            write_off_verts(
                net, atom['verts'], directory=directory,
                file_name=base_name + '_verts', file_type=file_type,
                color_scheme=scheme, color_map=cmap, color_limit=color_limit,
                target_cells=cell_targets, color_mode='cell',
            )

        if edges:
            write_edges(
                net, atom['edges'], directory=directory,
                file_name=base_name + '_edges', file_type=file_type,
                color_scheme=scheme, color_map=cmap, color_limit=color_limit,
                target_cells=cell_targets, color_mode='cell',
            )


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

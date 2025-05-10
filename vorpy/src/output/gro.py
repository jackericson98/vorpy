import os
import shutil


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


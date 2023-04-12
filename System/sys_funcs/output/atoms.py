import os
from System.sys_funcs.output.surfs import write_surfs
from System.sys_funcs.output.edges import write_edges
from System.sys_funcs.output.verts import write_verts


def write_pdb(atoms, name, sys=None, directory=None):
    """
    Creates a pdb file type in the current working directory
    :param directory:
    :param atoms: List of atom type objects for writing
    :param name: Name of the output file
    :param sys: System object used for writing the whole pbd file
    :return:
    """
    start_dir = None
    if directory is not None:
        start_dir = os.getcwd()
        os.chdir(directory)
    if atoms is None or len(atoms) == 0:
        return

    # Create the output file
    with open(name + ".pdb", 'w') as write_file:
        # Check to see if a system was provided
        if sys is not None and sys.base_file is not None:
            # Open the base file
            with open(sys.base_file, 'r') as f:
                read_file = f.readlines()
            # If the output is all atoms just copy the pdb
            if len(atoms) == len(sys.atoms):
                for line in read_file:
                    write_file.write(line)
                return
            # Otherwise, create a header and only export the relevant atoms
            else:
                # Write a header for the pdb
                write_file.write("HEADER  vorpy output - " + sys.name + " group " + name + " atoms\n")
                # Figure out what lines the atoms start on
                offset = 0
                while read_file[offset][:4].lower() != 'atom':
                    offset += 1
                # Grab the lines from the initial pdb
                for j in range(len(sys.atoms)):
                    if sys.atoms[j] in atoms:
                        write_file.write(read_file[j + offset])
        else:
            # Go through each atom in the system
            for i in range(len(atoms)):
                a = atoms[i]
                loc = ["{:.3f}".format(_) for _ in a.loc]
                # Get the information from the atom in writable format
                ser_num = " " * (5 - len(str(i+1))) + str(i + 1)
                name = a.name + " " * (4 - len(a.name))
                res = " " * (3 - len(a.residue)) + a.residue
                if a.chn.name.lower() == "zz" or a.chn.name.lower() == 'mol' or a.chn.name.lower() == 'sol':
                    chain = " "
                else:
                    chain = str(a.chn.name)
                res_seq = " " * (3 - len(str(a.res_seq))) + str(a.res_seq)
                loc_strs = [" " * (7 - len(_)) + _ for _ in loc]
                occupancy = " " * 5
                t_fact = " " * 5
                seg_id = a.seg_id + " " * (3 - len(a.seg_id))
                symbol = a.element
                charge = ''
                # Write the atom information
                write_file.write("ATOM  " + ser_num + " " + name + " " + res + " " + chain + res_seq + "    " +
                                 " ".join(loc_strs) + occupancy + t_fact + "      " + seg_id + symbol + charge + "\n")
    if start_dir is not None:
        os.chdir(start_dir)


def write_gro(atoms, name, sys, directory=None):
    # Change to the directory of specified
    if directory is not None and os.path.exists(directory):
        os.chdir(directory)
    # Create the title
    sys_name = sys.name
    # Open the file
    with open(name + '.gro', 'w') as f:
        # Write the header
        f.write("{}\n{:5d}\n".format(sys_name, len(atoms)))
        # Write the atoms information
        for atom in atoms:
            f.write("{:5d}{:5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(atom.res_seq, atom.res.name, atom.name,
                                                                         atom.num + 1, atom.loc[0], atom.loc[1],
                                                                         atom.loc[2]))
        # Write the box
        box = sys.net.box
        f.write("{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}\n".format(box[0][0], box[0][1], box[0][2], box[1][0],
                                                                            box[1][1], box[1][2]))


def write_atom_cells(atoms, directory, surfs=True, edges=False, verts=False):
    # Change to the directory
    os.chdir(directory)
    # Go through the atoms
    for i, atom in enumerate(atoms):
        # Check if the surfaces should be exported
        if surfs:
            # Write the surfaces
            write_surfs(atom.surfs, directory=directory, file_name=str(atom.num) + "_" + atom.name)
        # Check for verts
        if verts:
            write_verts(atom.verts, directory=directory, file_name=str(atom.num) + "_" + atom.name + "_verts")

        # Check for verts
        if edges:
            write_edges(atom.edges, directory=directory, file_name=str(atom.num) + "_" + atom.name + "_edges")

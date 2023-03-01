from System.sys_funcs.calcs import *
from System.sys_objs.atom import Atom, get_radius
from System.sys_funcs.draw import draw_edge, draw_vert
import os
import csv

####################################################### Exports ########################################################


def set_output_dir(sys, dir_name=None):
    """
    Sets the directory for the output data. If the directory exists add 1 to the end number
    :param sys: System to assign the output directory to
    :param dir_name: Name for the directory
    :return:
    """
    if not os.path.exists("/Data/User_data"):
        os.mkdir("./Data/User_data")
    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        if sys.vpy_dir is not None:
            dir_name = sys.vpy_dir + "/Data/User_data/" + sys.name
        else:
            dir_name = os.getcwd() + "/Data/User_data/" + sys.name
    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = '_' + str(i)
            # If no file with the system name exists change the string to empty
            if i == 0:
                i_str = ""
            # Try to create the directory
            os.mkdir(dir_name + i_str)
            break
        # If the file exists increment the counter and try creating the directory again
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    sys.dir = dir_name + i_str


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
                res = " " * (3 - len(a.mol_class)) + a.mol_class
                chain = str(a.chain) + " " * (1 - len(a.chain))
                if chain == "ZZ" or chain == 'MOL':
                    chain = "  "
                res_seq = " " * (3 - len(a.res_seq)) + a.res_seq
                loc_strs = [" " * (7 - len(_)) + _ for _ in loc]
                occupancy = " " * (5 - len(a.occupancy)) + a.occupancy
                t_fact = " " * (5 - len(a.t_fact)) + a.t_fact
                seg_id = a.seg_id + " " * (3 - len(a.seg_id))
                symbol = a.element
                charge = a.charge
                # Write the atom information
                write_file.write("ATOM  " + ser_num + " " + name + " " + res + " " + chain + res_seq + "    " + " ".join(loc_strs) +
                           occupancy + t_fact + "      " + seg_id + symbol + charge + "\n")
    if start_dir is not None:
        os.chdir(start_dir)


def write_verts(verts, file_name, atom_type=None, directory=None, pdb=False, spheres=False, color=None):
    """
    Creates a pdb file for vertex representation
    :param verts:
    :param file_name:
    :param atom_type:
    :param directory:
    :return:
    """
    # If no color is given, make the color random
    if color is None:
        color = [1, 0, 0]
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    if atom_type is None:
        atom_type = 'He'
    # If no surfaces are provided return
    if verts is None or len(verts) == 0:
        return
    if pdb:
        vert_atoms = [Atom(location=verts[i].loc, element=atom_type, mol_class="SOL", res_seq=str(i), name=atom_type, ) for i in range(len(verts))]
        # Write the pdb with the atom objects from the verts
        write_pdb(atoms=vert_atoms, name=file_name, directory=directory)
    else:

        num_verts, num_tris = 0, 0
        for vert in verts:
            draw_vert(vert, sphere=spheres)
            if spheres:
                # Go through and create each edge
                for i in range(len(verts)):
                    num_verts += len(verts[i].loc_points) + len(vert.sphere_points)
                    num_tris += len(verts[i].loc_tris) + len(vert.sphere_tris)
            else:
                num_verts = 8 * len(verts)
                num_tris = 8 * len(verts)
        # Create the file
        with open(file_name + ".off", 'w') as file:
            # Count the number of triangles and vertices there are
            # Write the numbers into the file
            file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
            # Go through the surfaces and add the points
            for i in range(len(verts)):
                sphere_points = []
                if spheres:
                    sphere_points = verts[i].sphere_points
                # Go through the points on the surface
                for point in verts[i].loc_points + sphere_points:
                    # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                    str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                    file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')

            num_verts, tri_count = 0, 0
            # Go through each surface and add the faces
            for i in range(len(verts)):
                vert = verts[i]
                sphere_tris = []
                if spheres:
                    sphere_tris = verts[i].sphere_tris
                tri_list = vert.loc_tris + sphere_tris
                # Go through the triangles in the surface
                for j in range(len(tri_list)):
                    # Get the triangle and colors
                    tri = tri_list[j]
                    # Add the triangle to the system file and the surface's file
                    str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                    file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                               str(color[1]) + " " + str(color[2]) + "\n")
                # Keep counting triangles for the system file
                if vert.sphere_points is None:
                    sphere_points_len = 0
                else:
                    sphere_points_len = len(vert.sphere_points)
                num_verts += len(vert.loc_points) + sphere_points_len


def write_edges(edges, file_name, color=None, directory=None):
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    # If no surfaces are provided return
    if edges is None or len(edges) == 0:
        return
    # If no color is given, make the color random
    if color is None:
        color = [1, 1, 1]
    # Check that the edge has been drawn
    for edge in edges:
        if edge.draw_points is None or edge.draw_tris is None:
            draw_edge(edge)
    num_verts, num_tris = 0, 0
    # Go through and create each edge
    for i in range(len(edges)):
        num_verts += len(edges[i].points) * 3
        num_tris += (len(edges[i].points) - 1) * 6
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for i in range(len(edges)):
            # Go through the points on the surface
            for point in edges[i].draw_points:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for i in range(len(edges)):
            edge = edges[i]
            # Go through the triangles in the surface
            for j in range(len(edge.draw_tris)):
                # Get the triangle and colors
                tri = edge.draw_tris[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(edge.draw_points)


def write_surfs(surfs, file_name, color_map='inferno', color_scheme='dist', color=False, directory=None):
    """
    Writes files given a list of surfaces into the current directory or the given one
    :param surfs: Surface object
    :param file_name: Name of the output file for the surfaces
    :param color: Color of the output surface
    :param directory:
    :return:
    """
    # Check to see if a directory is given
    if directory is not None:
        os.chdir(directory)
    # If no surfaces are provided return
    if surfs is None or len(surfs) == 0:
        return
    # If no color is given, make the color random
    if color is False:
        color = np.random.rand(3)
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        num_verts, num_tris = 0, 0
        for i in range(len(surfs)):
            if surfs[i].points is None:
                surfs[i].build()
            if color_map != surfs[i].color_map or color_scheme != surfs[i].color_scheme:
                surfs[i].color_tris(color_map=color_map, color_scheme=color_scheme)
            num_verts += len(surfs[i].points)
            num_tris += len(surfs[i].tris)
        # Write the numbers into the file
        file.write("OFF\n" + str(num_verts) + " " + str(num_tris) + " 0\n\n\n")
        # Go through the surfaces and add the points
        for i in range(len(surfs)):
            # Go through the points on the surface
            for point in surfs[i].points:
                # Add the point to the system file and the surface's file (rounded to 4 decimal points)
                str_point = [str(round(float(point[_]), 4)) for _ in range(3)]
                file.write(str_point[0] + " " + str_point[1] + " " + str_point[2] + '\n')
        num_verts, tri_count = 0, 0
        # Go through each surface and add the faces
        for i in range(len(surfs)):
            surf = surfs[i]
            # Go through the triangles in the surface
            for j in range(len(surfs[i].tris)):
                # Get the triangle and colors
                tri = surf.tris[j]
                colors = surf.tri_colors
                if colors is not None:
                    # If the surface is flat, average out the colors
                    color = colors[j]
                # Add the triangle to the system file and the surface's file
                str_tri = [str(tri[_] + num_verts) for _ in range(3)]
                file.write("3 " + str_tri[0] + " " + str_tri[1] + " " + str_tri[2] + " " + str(color[0]) + " " +
                           str(color[1]) + " " + str(color[2]) + "\n")
            # Keep counting triangles for the system file
            num_verts += len(surfs[i].points)


def export_mySys(sys):
    """
    Used to create and export the surfaces of a system as one file
    :param sys: System object
    :return:
    """
    # Write the surfaces
    write_surfs(sys.net.surfs, sys.name + "_system")


#################################################### Export Vertices ###################################################


def export_verts(net):
    """
    Exports a txt file with the vertex information for reloading later
    :param net: The network to interpret the vertex data from
    :return:
    """
    # Move to the correct output directory
    os.chdir(net.sys.dir)
    # Open the file for the vertices
    file = open(net.sys.name + "_verts.txt", 'w')
    # Create a header for the vertices file
    file.write(net.sys.name + " Vertices: \n")
    # Write the vertices
    for vert in net.verts:
        # Write the vertex
        file.write("VERT " + " ".join([str(_) for _ in vert.ndx]) + " " + " ".join([str(_) for _ in vert.loc]) + " " +
                   str(vert.rad) + "\n")
    # Write the end line for the file
    file.write("END")
    file.close()


#################################################### Export Network ####################################################


def export_net(net, output_surfs=True):
    # Create the file for export
    if net.sys.net_file is None:
        net.sys.net_file = net.sys.dir + "/" + net.sys.name + "_net.csv"
    # Create the file
    with open(net.sys.net_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write a separating line for the info and the surfaces points and tris
        writer.writerow(["Network", "Surface Resolution", "Maximum Vertex Resolution", "Box Size Multiplier",
                         "Calculate Surfaces?", "# of Vertices", "# of Edges", "# of Surfaces", "Surfaces Folder"])
        writer.writerow([net.sys.name] + [net.surf_res, net.max_vert, net.box_size, net.build_surfs,
                                       len(net.verts), len(net.edges), len(net.surfs), output_surfs])
        # Create a vertices header
        writer.writerow(["Vertex", "Loc - X", "Loc - Y", "Loc - Z", "Radius", "Atom 1", "Atom 2", "Atom 3", "Atom 4",
                         "Edge 1", "Edge 2", "Edge 3", "Edge 4", "Edge 5 (incorrect)", "Surface 1", "Surface 2", "Surface 3", "Surface 4",
                         "Surface 5", "Surface 6"])
        # Write the connections and location and radius for each vertex in the network
        for i in range(len(net.verts)):
            vert = net.verts[i]
            v_edges, v_surfs = [net.edges.index(_) for _ in vert.edges], [net.surfs.index(_) for _ in vert.surfs]
            writer.writerow([i] + [round(_, 3) for _ in vert.loc + [vert.rad]] + vert.ndx +
                            v_edges + [None] * (5 - len(v_edges)) + v_surfs + [None] * (6 - len(v_surfs)))

        # Create an edges header
        writer.writerow(["Edge", "Reference Surface", "Start Index", "End Index", "Atom 1", "Atom 2", "Atom 3",
                         "Vertex 1", "Vertex 2", "Surface 1", "Surface 2", "Surface 3"])
        # Write the connections and surface and points range information for each edge in the network
        for i in range(len(net.edges)):
            # Get the edge
            edge = net.edges[i]
            # Get the reference value for the edge
            edge_ref = [None, None, None]
            e_verts, e_surfs = [net.verts.index(_) for _ in edge.verts], [net.surfs.index(_) for _ in edge.surfs]
            # Write the edge information in the file
            writer.writerow([i] + edge_ref + edge.ndx + e_verts + [None] * (2 - len(e_verts)) + e_surfs +
                            [None] * (3 - len(e_surfs)))

        # Create a surfaces header
        writer.writerow(["Surface", "File", "Resolution", "Surface Area", "Curvature", "Atom 1", "Atom 2", "Function A", "Function B",
                         "Function C", "Function D", "Function E", "Function F", "Function G", "Function H",
                         "Function I", "Function J", "Function K", "Function d1", "Function d2", "Function d3"])
        # Write the connections and surface and points range information for each edge in the network
        for i in range(len(net.surfs)):
            # Get the surface
            surf = net.surfs[i]
            # Get the file address for the output points
            file_address = ""
            if surf.points is not None:
                file_address = "/surfs/" + "_".join([str(_) for _ in surf.ndx]) + ".off"
            if surf.res is None:
                surf.res = surf.net.surf_res
            if surf.sa is None:
                surf.sa = 0
            if surf.curv is None:
                surf.curv = 0
            # Write the surface information
            writer.writerow([i, file_address, surf.res, surf.sa, surf.curv, surf.ndx[0], surf.ndx[1]] + list(surf.func))
        # Check to see if the surfaces have been requested
        if output_surfs and net.build_surfs:
            # Create a surfaces folder and change to it
            if not os.path.exists(net.sys.dir + "/surfs"):
                os.mkdir(net.sys.dir + "/surfs")
            os.chdir(net.sys.dir + "/surfs")
            # Go through the surfaces 1 by one creating point files
            for surf in net.surfs:
                write_surfs([surf], "_".join([str(_) for _ in surf.ndx]))
    # Change back to the network file's directory
    os.chdir(net.sys.dir)


def export_net_info(net):
    # Open the file
    file = open(net.sys.name + "_net_info.txt", 'w')
    # Write the header
    file.write(net.sys.name + " Network")
    # Write the atom information
    for i in range(len(net.atoms)):
        file.write("{} - cell volume = {}, cell surface area {}\n".format(net.sys.atom_names[i], net.atoms[i].vol, net.atoms[i].sa))
    # write the surface information
    for i in range(len(net.surfs)):
        file.write("Surface {}-{} - Surface area = {}\n".format(net.surfs[i].ndx[0], net.surfs[i].ndx[1], net.surfs[i].sa))
    file.close()


############################################ Pymol Scripts #############################################################


def set_pymol_atoms(sys):
    """
    Creates a script to set the radii of the spheres in pymol
    :param sys:
    :return:
    """
    # Create the file
    file = open('set_atoms.pml', 'w')
    # Write the change radii script for the system's set atomic radii
    for i in range(len(sys.radii[0])):
        if sys.radii[1] is not None:
            file.write("alter (elem {}), vdw={}\n".format(sys.radii[0][i], sys.radii[1][i]))
    # Rebuild the system
    file.write("\nrebuild")
    file.close()

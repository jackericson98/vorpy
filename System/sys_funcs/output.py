from System.sys_funcs.calcs import *
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
    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        dir_name = sys.vorpy_directory + "/Data/User_data/" + sys.name
    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = str(i)
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


def write_pdb(atoms, name, sys=None):
    """
    Creates a pdb file type in the current working directory
    :param atoms: List of atom type objects for writing
    :param name: Name of the output file
    :param sys: System object used for writing the whole pbd file
    :return:
    """

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
                loc = [str(round(_, 3)) for _ in a.loc]
                # Get the information from the atom in writable format
                ser_num = " " * (5 - len(str(i+1))) + str(i + 1)
                name = a.name + " " * (4 - len(a.name))
                res = " " * (3 - len(a.res)) + a.res
                chain = str(a.mol) + " " * (1 - len(a.mol))
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


def write_surfs(surfs, file_name, color=None, directory=None):
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
    if color is None:
        color = np.random.rand(3)
    # Create the file
    with open(file_name + ".off", 'w') as file:
        # Count the number of triangles and vertices there are
        num_verts, num_tris = 0, 0
        for i in range(len(surfs)):
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
            for tri in surfs[i].tris:
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


def export_old_net(net, point_res=None):
    """
    Used to store pre-calculated surfaces in whatever directory the program is in
    :param net: The network object to export the data from
    :param point_res: Holds the number of decimal places to output. If None it is set to 3
    :return:
    """
    # Move to the output directory
    os.chdir(net.sys.dir)
    # Create the network file
    file = open(net.sys.name + "_network.txt", 'w')
    # Write the general information about the system
    file.write("NETW " + str(net.surf_res) + " " + str(net.max_vert) + " " + str(net.box_size) + " " + str(net.my_time)
               + " " + str(net.cpu_time) + " " + str(net.sol_verts) + " " + str(net.curved_faces) + " " +
               str(net.flat_faces) + " " + str(len(net.verts)) + " " + str(len(net.edges)) + " " +
               str(len(net.surfs)) + "\n")
    # Set the resolution of the network's output points
    if point_res is None:
        point_res = 3

    # Write Objects:

    # Write atoms
    # Go through the atoms in the network
    for atom in net.atoms:
        # Get the atom's box
        box = [str(_) for _ in atom.box]
        # Write atoms information: index, box, cell volume
        file.write("ATOM " + " " + str(net.atoms.index(atom)) + " " + ' '.join(box) + " " + str(atom.vol) + "\n")
        # Get the vertex, edge and surface index information
        vert_ndxs = [str(net.verts.index(vert)) for vert in atom.verts]
        edge_ndxs = [str(net.edges.index(edge)) for edge in atom.edges]
        surf_ndxs = [str(net.surfs.index(surf)) for surf in atom.surfs]
        # Write the object indices
        file.write("ACON " + " ".join(vert_ndxs) + "\n")
        file.write("ACON " + " ".join(edge_ndxs) + "\n")
        file.write("ACON " + " ".join(surf_ndxs) + "\n")
    # Write a separating line
    file.write("\n")

    # Write vertices
    for vert in net.verts:
        # Get the normal information
        loc, ndx = [str(round(_, point_res)) for _ in vert.loc], [str(_) for _ in vert.ndx]
        # Get the doublet information
        loc2, rad2 = [""], ""
        if vert.doublet:
            loc2, rad2 = [str(round(_, point_res)) for _ in vert.loc2], str(round(vert.rad2, point_res))
        # Write the vertex information
        file.write("VERT " + str(net.verts.index(vert)) + " " + " ".join(ndx) + " " + " ".join(loc) + " " +
                   str(vert.rad) + " " + str(vert.doublet) + " " + " ".join(loc2) + " " + rad2 + '\n')
        # Get the edge and surface index information
        edge_ndxs = [str(net.edges.index(edge)) for edge in vert.edges]
        surf_ndxs = [str(net.surfs.index(surf)) for surf in vert.surfs]
        # Write the connection information
        file.write("VCON " + " ".join(edge_ndxs) + "\n")
        file.write("VCON " + " ".join(surf_ndxs) + "\n")
    # Write a separating line
    file.write("\n")

    # Write edges
    for edge in net.edges:
        # If the edge location is None get a location
        if edge.loc is None:
            edge.loc, edge.rad = calc_circ(edge.atoms)
        # Get the atom's box
        ndx, loc = [str(_) for _ in edge.ndx], [str(round(_, point_res)) for _ in edge.loc]
        # Make sure the points are interpretable
        if edge.pv1 is None:
            edge.pv0, edge.pv1 = [np.inf, np.inf, np.inf], [np.inf, np.inf, np.inf]
        rad, pv0, pv1 = str(round(edge.rad, point_res)) + " ", [str(round(_, point_res)) for _ in edge.pv0], \
                        [str(round(_, point_res)) for _ in edge.pv1]
        # Write Edge information: index, location, radius, end points
        file.write("EDGE " + str(net.edges.index(edge)) + " " + " ".join(ndx) + " " + " ".join(loc) + " " + rad +
                   " ".join(pv0) + " " + " ".join(pv1) + " " + str(edge.doublet) + "\n")
        # Get the vertex, edge and surface index information
        vert_ndxs = [str(net.verts.index(vert)) for vert in edge.verts]
        surf_ndxs = [str(net.surfs.index(surf)) for surf in edge.surfs]
        # Write the object indices
        file.write("ECON " + " ".join(vert_ndxs) + "\n")
        file.write("ECON " + " ".join(surf_ndxs) + "\n")
        # Check to make sure the edge has points
        if edge.points is None:
            edge.build()
        # Go through the points along the edge
        for point in edge.points:
            # Add the points of the edge to the edge file
            file.write("EPNT " + str(round(point[0], point_res)) + " " + str(round(point[1], point_res)) + " " +
                       str(round(point[2], point_res)) + "\n")
    # Write a separating line
    file.write("\n")

    # Write surfaces
    for surf in net.surfs:
        # Write the main edge information
        ndx, rn = [str(_) for _ in surf.ndx], [str(_) for _ in surf.rn]
        file.write("SURF " + str(net.surfs.index(surf)) + " " + " ".join(ndx) + " " + " ".join(rn) + " " +
                   str(surf.sa) + " " + str(len(surf.perimeter)) + '\n')
        # Get the vertex, edge and surface index information
        vert_ndxs = [str(net.verts.index(vert)) for vert in surf.verts]
        edge_ndxs = [str(net.edges.index(edge)) for edge in surf.edges]
        # Write the object indices
        file.write("SCON " + " ".join(vert_ndxs) + "\n")
        file.write("SCON " + " ".join(edge_ndxs) + "\n")
        # Go through the points along the perimeter of the surface
        for i in range(len(surf.points)):
            # Add the points of the edge to the edge file
            file.write("SPNT " + " ".join([str(round(_, point_res)) for _ in surf.points[i]]) + "\n")
        # Go through the triangles in the surface's list of triangles
        for i in range(len(surf.tris)):
            # Add the triangles to the list of surface triangles
            file.write("STRI " + " ".join([str(_) + " " for _ in surf.tris[i]]) + "\n")
    # Write the end line
    file.write('END')
    file.close()


def export_csv_surfs(net):
    # Create a surfaces folder and change to it
    os.mkdir(net.sys.dir + "/csv_surfs")
    os.chdir(net.sys.dir + "/csv_surfs")
    # Go through the surfaces 1 by one creating point files
    for surf in net.surfs:
        # Create the surface file
        with open(os.getcwd() + "/" + "_".join([str(_) for _ in surf.ndx]) + ".off", 'w', newline='') as surf_file:
            surf_writer = csv.writer(surf_file)
            surf_writer.writerow(["Surface", "# of Points", "# of Triangles"])
            surf_writer.writerow([net.surfs.index(surf), len(surf.points), len(surf.tris)])
            # Write the header for the points
            surf_writer.writerow(["Point", "Loc - X", "Loc - Y", "Loc - Z"])
            # Go through the points
            for j in range(len(surf.points)):
                # Write the point information
                surf_writer.writerow([j, surf.points[j][0], surf.points[j][1], surf.points[j][2]])
            # Write the triangles header
            surf_writer.writerow(["Triangle", "Point 1", "Point 2", "Point 3"])
            # Go through the triangles
            for j in range(len(surf.tris)):
                # Write the triangle information
                surf_writer.writerow([j, surf.tris[j][0], surf.tris[j][1], surf.tris[j][2]])


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
        writer.writerow([net.sys.name] + [net.surf_res, net.max_vert, net.box_size, net.calc_surfs,
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
            if edge.ref is None:
                if len(edge.surfs) > 0:
                    try:
                        ndx_1 = edge.surfs[0].points.index(edge.points[0])
                        ndx_2 = ndx_1 + len(edge.points)
                        surf_ndx = net.surfs.index(edge.surfs[0])
                        # Set the reference value
                        edge.ref = [surf_ndx, ndx_1, ndx_2]
                    except IndexError:
                        edge.ref = [None, None, None]
                    except ValueError:
                        edge.ref = [None, None, None]
                    except AttributeError:
                        edge.ref = [None, None, None]
                else:
                    edge.ref = [None, None, None]
            e_verts, e_surfs = [net.verts.index(_) for _ in edge.verts], [net.surfs.index(_) for _ in edge.surfs]
            # Write the edge information in the file
            writer.writerow([i] + edge.ref + edge.ndx + e_verts + [None] * (2 - len(e_verts)) + e_surfs +
                            [None] * (3 - len(e_surfs)))

        # Create a surfaces header
        writer.writerow(["Surface", "File", "Surface Area", "Curvature", "Atom 1", "Atom 2", "Function A", "Function B",
                         "Function C", "Function D", "Function E", "Function F", "Function G", "Function H",
                         "Function I", "Function J", "Function K", "Function d1", "Function d2", "Function d3"])
        # Write the connections and surface and points range information for each edge in the network
        for i in range(len(net.surfs)):
            # Get the surface
            surf = net.surfs[i]
            # Get the file address for the output points
            file_address = ""
            if output_surfs and net.calc_surfs:
                file_address = net.sys.dir + "/surfs/" + "_".join([str(_) for _ in surf.ndx]) + ".off"
            # Write the surface information
            writer.writerow([i, file_address, surf.sa, surf.curv, surf.ndx[0], surf.ndx[1]] + list(surf.func[:11]) + list(surf.func[11]))
        # Check to see if the surfaces have been requested
        if output_surfs and net.calc_surfs:
            # Create a surfaces folder and change to it
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


def export_net1(net, verts_only=False):
    """
    An efficient storage of network information. 
    :param net:
    :param verts_only:
    :return:
    """



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

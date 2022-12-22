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
    file = open(name + ".pdb", 'w')
    # Check to see if a system was provided
    if sys is not None and sys.base_file is not None:
        try:
            # Open the base file
            with open(sys.base_file, 'r') as f:
                base_file = f.readlines()
            # Copy the lines
            for line in base_file:
                file.write(line)
            file.close()
            return
        except FileNotFoundError:
            return
    # Go through each atom in the system
    for i in range(len(atoms)):
        a = atoms[i]
        loc = [str(round(_, 3)) for _ in a.loc]
        # Get the information from the atom in writable format
        ser_num = " " * (5 - len(str(i+1))) + str(i + 1)
        name = a.name + " " * (4 - len(a.name))
        res = " " * (3 - len(a.res)) + a.res
        chain = str(a.mol) + " " * (1 - len(a.mol))
        if chain == "ZZ":
            chain = "  "
        res_seq = " " * (3 - len(a.res_seq)) + a.res_seq
        loc_strs = [" " * (7 - len(_)) + _ for _ in loc]
        occupancy = " " * (5 - len(a.occupancy)) + a.occupancy
        t_fact = " " * (5 - len(a.t_fact)) + a.t_fact
        seg_id = a.seg_id + " " * (3 - len(a.seg_id))
        symbol = a.element
        charge = a.charge
        # Write the atom information
        file.write("ATOM  " + ser_num + " " + name + " " + res + " " + chain + res_seq + "     " + " ".join(loc_strs) +
                   occupancy + t_fact + "        " + seg_id + symbol + charge + "\n")
    file.close()


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
    file = open(file_name + ".off", 'w')
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
    file.close()


def export_mySys(sys):
    """
    Used to create and export the surfaces of a system as one file
    :param sys: System object
    :return:
    """
    # Write the surfaces
    write_surfs(sys.net.surfs, sys.name + "_system")


#################################################### Main Requests #####################################################


def export_iface(groups, info_file=False, interface_atoms=False):
    """
    Exports the information from the given interface as a txt file
    :param groups: Group objects for the interface
    :param info_file: Whether to export a txt file with info on the interface or not
    :param interface_atoms: Whether to export a pdb file with the atoms around the interface or not
    :return:
    """
    # Get the groups
    g0, g1 = groups
    # Set the interface name
    interface_name = g0.name + "_" + g1.name + "_interface"
    # Move to the output directory
    os.chdir(g0.net.sys.dir)
    # Make sure the two groups are best friends
    if g0.bff != g1 or g1.bff != g0:
        # Set the groups as friends
        g0.bff, g1.bff = g1, g0
        # Get the information for the two atoms
        g0.get_info()
    # Create and move to the interface directory
    os.mkdir(os.getcwd() + "/" + interface_name)
    os.chdir(os.getcwd() + "/" + interface_name)
    # Write the surfaces for the interface
    write_surfs(g0.iface_surfs, interface_name)
    # Check to see of the user wants to export the interface's atoms
    if interface_atoms:
        # Get the two sets of interface atoms
        write_pdb(g0.iface_atoms, interface_name + "_" + g0.name + "_atoms")
        write_pdb(g1.iface_atoms, interface_name + "_" + g1.name + "_atoms")
    # Check to see if the user wants to export the interface's information
    if info_file:
        info = open(interface_name + "_info.txt", 'w')
        info.write("Interface between " + g0.name + " and " + g1.name + " : \n")
        info.write("Number of Surfaces: " + str(len(g0.iface_surfs)))
        info.write("Surface Area: " + str(g0.iface_sa))
        info.close()


def export_body(group, info_file=False, outer_atoms=False):
    """
    Exports the information from the given body as a txt file
    :param group: Group object
    :param info_file: Whether to export an information file about the group or not
    :param outer_atoms: Whether to export a pdb file of the outer atoms or not
    :return:
    """
    # Move to the output directory
    os.chdir(group.net.sys.dir)
    no_dir = True
    inc = 0
    while no_dir:
        # Create the directory for the storage of exported group data
        my_str = str(inc)
        if inc == 0:
            my_str = ""
        try:
            os.mkdir(group.name + my_str)
            no_dir = False
        except FileExistsError:
            inc += 1

    os.chdir(os.getcwd() + "/" + group.name)
    # If the group name is empty, name it
    if group.name is None:
        group.set_name()
    # Write the surfaces for the interface
    write_surfs(group.body_surfs, group.name)
    # Check to see of the user wants to export the interface's atoms
    if outer_atoms:
        write_pdb(group.outer_body_atoms, group.name + "_outside_atoms")
        write_pdb(group.surr_body_atoms, group.name + "_surrounding_atoms")
    # Check to see if the user wants to export the interface's information
    if info_file:
        info = open("cell_" + group.name + "_info.txt", 'w')
        info.write(group.name + " body: \n")
        info.write("Number of atoms: " + str(len(group.atoms)) + "\n")
        info.write("Volume: " + str(group.body_vol) + "\n")
        info.write("Surface Area: " + str(group.body_sa) + "\n")
        info.close()
    os.chdir("..")


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


def export_net(net, point_res=None):
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
        file.write("ATOM " + " " + str(net.atoms.index(atom)) + " " + ' '.join(box) + " " + str(atom.cell_vol) + "\n")
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


def export_net_fast(net, light=False, separate_surfs=False):
    # Create the file for export
    if net.sys.net_file is None:
        net.sys.net_file = net.sys.dir + "/" + net.sys.name + "_net.csv"
    # Create the file
    with open(net.sys.net_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write a separating line for the info and the surfaces points and tris
        writer.writerow(["Network", "Surface Resolution", "Maximum Vertex Resolution", "Box Size Multiplier",
                         "Calculate Solute Vertices?", "# of Atoms", "# of Vertices", "# of Edges", "# of Surfaces", "Light"])
        writer.writerow([net.sys.name] + [net.surf_res, net.max_vert, net.box_size, net.sol_verts, len(net.atoms),
                                       len(net.verts), len(net.edges), len(net.surfs), separate_surfs])
        # Create an atoms header
        writer.writerow(["Atom", "Loc - X", "Loc - Y", "Loc - Z", "Radius", "Vertices", "Edges", "Surfaces"])
        # Write the connections and location and radius information for each atom in the network
        for i in range(len(net.atoms)):
            atom = net.atoms[i]
            writer.writerow([net.atoms.index(atom)] + [round(_, 3) for _ in atom.loc + [atom.rad]] +
                            [[net.verts.index(_) for _ in atom.verts]] +
                            [[net.edges.index(_) for _ in atom.edges]] +
                            [[net.surfs.index(_) for _ in atom.surfs]])
        # Create a vertices header
        writer.writerow(["Vertex", "Loc - X", "Loc - Y", "Loc - Z", "Radius", "Atoms", "Edges", "Surfaces"])
        # Write the connections and location and radius for each vertex in the network
        for i in range(len(net.verts)):
            vert = net.verts[i]
            writer.writerow([i] + [round(_, 3) for _ in vert.loc + [vert.rad]] +
                            [[net.atoms.index(_) for _ in vert.atoms], [net.edges.index(_) for _ in vert.edges],
                             [net.surfs.index(_) for _ in vert.surfs]])
        # Create an edges header
        writer.writerow(["Edge", "Reference Surface", "Starting Surface Point Index", "Ending Surface Point Index",
                         "Atoms", "Vertices", "Surfaces"])
        # Write the connections and surface and points range information for each edge in the network
        for i in range(len(net.edges)):
            edge = net.edges[i]
            if edge.ref is None:
                try:
                    ndx_1 = edge.surfs[0].points.index[edge.points[0]]
                    ndx_2 = ndx_1 + len(edge.points)
                    surf_ndx = net.surfs.index(edge.surfs[0])
                except IndexError:
                    ndx_1, ndx_2, surf_ndx = None, None, None

                edge.ref = [surf_ndx, ndx_1, ndx_2]
            writer.writerow([i] + edge.ref + [[net.atoms.index(_) for _ in edge.atoms],
                            [net.verts.index(_) for _ in edge.verts], [net.surfs.index(_) for _ in edge.surfs]])
        # Create a surfaces header
        writer.writerow(["Surface", "Function Coefficient", "Atoms", "Vertices", "Edges"])
        # Write the connections and surface and points range information for each edge in the network
        for i in range(len(net.surfs)):
            surf = net.surfs[i]
            writer.writerow([i] + [surf.func, [net.atoms.index(_) for _ in surf.atoms],
                                   [net.verts.index(_) for _ in surf.verts], [net.edges.index(_) for _ in surf.edges]])
            if not light and not separate_surfs:
                writer.writerow(["Points: "] + [[round(_[j], 3) for j in range(3)] for _ in surf.points])
                writer.writerow(["Triangles:"] + [[_ for _ in tri] for tri in surf.tris])
    # If the user wants separate surfs help them
    if separate_surfs:
        # Create a directory for the separate network surface
        os.mkdir(surf.net.sys.dir + "/net_surfs")
        os.chdir(surf.net.sys.dir + "/net_surfs")
        # Create a file for each of the surfaces
        for surf in net.surfs:
            # Create the file
            with open("surf_" + str(surf.ndx[0]) + "_" + str(surf.ndx[1]) + ".csv") as sf:
                # Create the writer
                writer = csv.writer(sf)
                # Write each point
                for point in surf.points:
                    # Write the point
                    writer.writerow(["s"] + [round(_, 3) for _ in point])
                # Write each triangle
                for tri in surf.tris:
                    # Write the triangle
                    writer.writerow(['t'] + [round(_, 3) for _ in tri])
        # Change back to the system directory
        os.chdir(net.sys.dir)


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

from System.calcs import *
import os

####################################################### Exports ########################################################


# Set output directory function. prevents the system from making duplicate output directories
def set_output_dir(sys, dir_name=None):
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
    sys.output_directory = dir_name + i_str


# Create pdb method. Creates a pdb file type in the current working directory
def write_pdb(atoms, name, sys=None):
    # Create the output file
    file = open(name + ".pdb", 'w')
    # Check to see if a system was provided
    if sys is not None and sys.base_file is not None:
        # Open the base file
        with open(sys.base_file, 'r') as f:
            base_file = f.readlines()
        # Copy the lines
        for line in base_file:
            file.write(line)
        # Return the file
        return file
    # Go through each atom in the system
    for i in range(len(atoms)):
        a = atoms[i]
        loc = [str(round(_, 3)) for _ in a.loc]
        # Get the information from the atom in writable format
        ser_num = " " * (5 - len(str(i+1))) + str(i + 1)
        name = a.name + " " * (4 - len(a.name))
        res = " " * (3 - len(a.res)) + a.res
        chain = str(a.chain) + " " * (1 - len(a.chain))
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
        file.write("ATOM  " + ser_num + " " + name + " " + res + " " + chain + res_seq + "    " + " ".join(loc_strs) +
                   occupancy + t_fact + "        " + seg_id + symbol + charge + "\n")


# Write surfaces function. Writes files given a list of surfaces
def write_surfs(surfs, file_name, color=None):
    # If no color is given, make the color white
    if color is None:
        color = [1, 1, 1]
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


# Export system function. Used to create and export the surfaces of a system as one file
def export_mySys(sys):
    # Write the surfaces
    write_surfs(sys.net.surfs, sys.name + "_system")


#################################################### Main Requests #####################################################


# Export interface information function. Exports the information from the given interface as a txt file
def export_interface(groups, info_file=False, interface_atoms=False):
    # Get the groups
    g0, g1 = groups
    # Set the interface name
    interface_name = g0.name + "_" + g1.name + "_interface"
    # Move to the output directory
    os.chdir(g0.net.sys.output_directory)
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
    write_surfs(groups[0].interface_surfs, interface_name)
    # Check to see of the user wants to export the interface's atoms
    if interface_atoms:
        # Get the two sets of interface atoms
        write_pdb(g0.interface_atoms, interface_name + "_" + g0.name + "_atoms")
        write_pdb(g1.interface_atoms, interface_name + "_" + g1.name + "_atoms")
    # Check to see if the user wants to export the interface's information
    if info_file:
        info = open(interface_name + "_info.txt", 'w')
        info.write("Interface between " + g0.name + " and " + g1.name + " : \n")
        info.write("Number of Surfaces: " + str(len(g0.interface)))
        info.write("Surface Area: " + str(g0.interface_sa))


# Export interface information function. Exports the information from the given body as a txt file
def export_body(group, info_file=False, outer_atoms=False):
    # Move to the output directory
    os.chdir(group.net.sys.output_directory)
    # Write the surfaces for the interface
    write_surfs(group.body_surfs, group.name)
    # Check to see of the user wants to export the interface's atoms
    if outer_atoms:
        write_pdb(group.outer_body_atoms, group.name + "_outside_atoms")
        write_pdb(group.surr_body_atoms, group.name + "_surrounding_atoms")
    # Check to see if the user wants to export the interface's information
    if info_file:
        info = open(group + "_info.txt", 'w')
        info.write(group.name + " body: \n")
        info.write("Number of atoms: " + str(len(group.atoms)))
        info.write("Volume: " + str(group.body_vol))
        info.write("Surface Area: " + str(group.body_sa))


#################################################### Export Network ####################################################


# Export surfaces function. Used to store pre-calculated surfaces in whatever directory the program is in
def export_net(net):

    # Move to the output directory
    os.chdir(net.sys.output_directory)
    # Create the network file
    file = open(net.sys.name + "_network.txt", 'w')
    # Write the general information about the system
    file.write("NETW " + str(net.min_dist) + " " + str(net.beta_val) + " " + str(net.box_size) + " " + str(net.my_time)
               + " " + str(net.cpu_time) + " " + str(net.sol_verts) + " " + str(net.curved_faces) + " " +
               str(net.flat_faces) + "\n")

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
        loc, ndx = [str(_) for _ in vert.loc], [str(_) for _ in vert.ndx]
        # Get the doublet information
        loc2, rad2 = [""], ""
        if vert.doublet:
            loc2, rad2 = [str(_) for _ in vert.loc2], str(vert.rad2)
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
        ndx, loc = [str(_) for _ in edge.ndx], [str(_) for _ in edge.loc]
        # Make sure the points are interpretable
        if edge.pv1 is None:
            edge.pv0, edge.pv1 = [np.inf, np.inf, np.inf], [np.inf, np.inf, np.inf]
        rad, pv0, pv1 = str(edge.rad) + " ", [str(_) for _ in edge.pv0], [str(_) for _ in edge.pv1]
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
            file.write("EPNT " + str(point[0]) + " " + str(point[1]) + " " + str(point[2]) + "\n")
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
            file.write("SPNT " + " ".join([str(_) for _ in surf.points[i]]) + "\n")
        # Go through the triangles in the surface's list of triangles
        for i in range(len(surf.tris)):
            # Add the triangles to the list of surface triangles
            file.write("STRI " + " ".join([str(_) + " " for _ in surf.tris[i]]) + "\n")
    # Write the end line
    file.write('END')

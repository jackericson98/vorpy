import os
import csv
from System.sys_funcs.output.surfs import write_surfs


def export_net_logs(net):
    # Open the file
    with open(net.sys.name + "_logs.csv", 'w') as l:
        # Create the csv writer
        logs = csv.writer(l)
        # Write the atom header
        logs.writerow(["Atoms"])
        # Write the column labels
        logs.writerow(["index", "name", "volume", "surface area", "neighbors"])
        # Go through the atoms in the system
        for i, atom in enumerate(net.atoms):
            a_surfs = [_[0] if _[0] != atom.num else _[1] for _ in [_.ndx for _ in atom.surfs]]
            logs.writerow([i, atom.name, atom.vol, atom.sa, *a_surfs])
        # Write the surfaces header
        logs.writerow(["Surfaces"])
        # Write the surface column labels
        logs.writerow(["index", "atom0", "atom1", "surface area", "curvature", "vol a0", "vol a1"])
        # Go through the surfaces in the system and write their information
        for i, surf in enumerate(net.surfs):
            # Write the information for the surface
            logs.writerow([i, *surf.ndx, surf.sa, surf.curv, *surf.vols])
        # Write the edges header
        logs.writerow(["Edges"])
        # Write the edges headers
        logs.writerow(["index", "atom0", "atom1", "atom2", "length"])
        # Go through the edges in the network
        for i, edge in enumerate(net.edges):
            # Write the data for the edge
            logs.writerow([i, *edge.ndx, edge.length])
        # Write the vertices header
        logs.writerow(["Vertices"])
        # Write the vertices data labels
        logs.writerow(["index", "atom0", "atom1", "atom2", "atom3", "x", "y", "z", "r"])
        # Go through the vertices
        for i, vert in enumerate(net.verts):
            # Write the vertex information line
            logs.writerow([i, *vert.atoms, *vert.loc, vert.rad])
        # Write the connections header
        logs.writerow(["Vertex Connections"])
        # Write the connection data labels
        logs.writerow(["index", "edge0", "edge1", "edge2", "edge3", "surf0", "surf1", "surf2", "surf3", "surf4", "surf5"])
        # Go through the vertices
        for i, vert in enumerate(net.verts):
            # Write the data
            logs.writerow([i] + vert.edges + [-1] * (4 - len(vert.edges)) + vert.surfs)


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
                         "Edge 1", "Edge 2", "Edge 3", "Edge 4", "Edge 5 (incorrect)", "Surface 1", "Surface 2",
                         "Surface 3", "Surface 4", "Surface 5", "Surface 6"])
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
        edge_ref = [None, None, None]
        for i in range(len(net.edges)):
            # Get the edge
            edge = net.edges[i]
            # Get the reference value for the edge
            e_verts, e_surfs = [net.verts.index(_) for _ in edge.verts], [net.surfs.index(_) for _ in edge.surfs]
            # Write the edge information in the file
            writer.writerow([i] + edge_ref + edge.ndx + e_verts + [None] * (2 - len(e_verts)) + e_surfs +
                            [None] * (3 - len(e_surfs)))

        # Create a surfaces header
        writer.writerow(["Surface", "File", "Resolution", "Surface Area", "Curvature", "Atom 1", "Atom 2", "Function A",
                         "Function B", "Function C", "Function D", "Function E", "Function F", "Function G",
                         "Function H", "Function I", "Function J", "Function K", "Function d1", "Function d2",
                         "Function d3"])
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
from System.sys_funcs.calcs.calcs import round_func
import os
import csv


def write_logs(group, net_name=None, round_to=3):
    """
    Writes a file for the networks logs
    :param net: Network for logs
    :param round_to: Where to round to
    """
    net = group.net
    # Set the network name for the logs file to be exported
    net_name = '' if net_name is None else '_' + str(net_name)
    # Create the round function
    r = round_func(round_to)
    # Open the file
    with open(group.sys.name + net_name + "_logs.csv", 'w') as log_file:
        # Create the csv writer
        lg_fl = csv.writer(log_file, lineterminator='\n')
        # Write the build information header
        lg_fl.writerow(["build informaiton"])
        # Write the build information labels
        lg_fl.writerow(["name", "network type", "surface resolution", "box size", "max vert", "Total Time", "vert time",
                       "connect time", "surf time", "analysis time", "max vertex"])
        lg_fl.writerow([group.sys.name, net.settings['net_type'], net.settings['surf_res'], net.settings['box_size'],
                        net.settings['max_vert'], r(net.metrics['tot']), r(net.metrics['vert']), r(net.metrics['con']),
                        r(net.metrics['surf']), r(net.metrics['anal']), r(max(net.verts['rad']))])
        # Write the group information header
        lg_fl.writerow(["group information"])
        # Write the group information labels
        lg_fl.writerow(["Name", "Volume", "Surface Area", "Density", "Center of Mass", "VDW Volume",
                        "VDW Center of Mass"])
        # Write the group information
        if group.sa is None:
            group.get_info()
        lg_fl.writerow([group.name, r(group.vol), r(group.sa), r(group.density), [r(_) for _ in group.com],
                        r(group.vdw_vol), [r(_) for _ in group.vdw_com]])
        # Write the atom header
        lg_fl.writerow(["Atoms"])
        # Write the column labels
        lg_fl.writerow(["Index", "Name", "X", "Y", "Z", "Radius", "Volume", "Surface Area", "Complete Cell?",
                        "Maximum Curvature", "Average Surface Curvature", "Sphericity", "Isometric Quotient",
                        "Number of Neighbors", "Closest Neighbor", "Closest Neighbor Distance",
                        "Neighbor Distance Average", "Neighbor Distance RMSD", "Minimum Point Distance",
                        "Maximum Point Distance", "Number of Overlaps", "Contact Area", "Non-Overlap Volume",
                        "Overlap Volume", "Center of Mass", "Moment of Inertia Tensor", "Bounding Box", "neighbors"])
        # Go through the atoms in the system
        for i, atom in net.balls.iterrows():
            if atom['sa'] == 0:
                continue
            if atom['complete']:
                print(atom['moi'])
                nbrs = [satoms[0] if satoms[0] != atom['num'] else satoms[1] for satoms in [net.surfs['balls'][_] for _ in atom['surfs']]]
                lg_fl.writerow([i, atom['name'], atom['loc'][0], atom['loc'][1], atom['loc'][2], atom['rad'],
                                r(atom['vol']), r(atom['sa']), atom['complete'], r(atom['max_curv']),
                                r(atom['avg_surf_curv']), r(atom['sphericity']), r(atom['isometric_quotient']),
                                atom['neighbor_number'], atom['near_neighbor'], atom['near_neighbor_distance'],
                                r(atom['neighbor_distance_average']), r(atom['neighbor_distance_rmsd']),
                                r(atom['min_spike']), r(atom['max_spike']), atom['num_olaps'], r(atom['contact_area']),
                                r(atom['non_olap_vol']), r(atom['vdw_vol']), [r(_) for _ in atom['com']],
                                r(atom['moi']), [[r(_) for _ in atom['b_box'][0]], [r(_) for _ in atom['b_box'][1]]],
                                *nbrs])
        # Write the surfaces header
        lg_fl.writerow(["Surfaces"])
        # Write the surface column labels
        lg_fl.writerow(["Index", "Ball 1", "Ball 2", "Surface Area", "Curvature", "Ball 1 Volume Contribution",
                        "Ball 2 Volume Contribution"])
        # Go through the surfaces in the system and write their information
        for i, surf in net.surfs.iterrows():
            # Write the information for the surface
            lg_fl.writerow([i, *surf['balls'], r(surf['sa']), r(surf['curv']), r(surf['vols'][surf['balls'][0]]),
                            r(surf['vols'][surf['balls'][1]])])
        # Write the edges header
        lg_fl.writerow(["Edges"])
        # Write the edges headers
        lg_fl.writerow(["index", "atom0", "atom1", "atom2", "length"])
        # Go through the edges in the network
        for i, edge in net.edges.iterrows():
            # Write the data for the edge
            lg_fl.writerow([i, *edge['balls'], r(edge['length'])])
        # Write the vertices header
        lg_fl.writerow(["Vertices"])
        # Write the vertices data labels
        lg_fl.writerow(["index", "atom0", "atom1", "atom2", "atom3", "x", "y", "z", "r"])
        # Go through the vertices
        for i, vert in net.verts.iterrows():
            # Write the vertex information line
            lg_fl.writerow([i, *vert['balls'], *r(vert['loc']), r(vert['rad'])])


def write_net(group, file_name=None, round_to=3):
    """
    Exports a network checkpoint file to be loaded later
    :param net: Network object for export
    :param file_name: Name of the output file
    :param round_to: Number of decimal places to round the values of the network to
    :return: Outputs a .csv network file
    """
    net = group.net
    # Set up the round function
    r = round_func(round_to)
    # Create the file for export
    if file_name is None:
        file_name = group.sys.files['dir'] + "/" + group.sys.name + "_net.csv"
    group.sys.files['net_file'] = file_name
    # Create the file
    with open(file_name, 'w', newline='') as f:

        # Create the writer object
        nt_fl = csv.writer(f)
        # Write a separating line for the info and the surfaces points and tris
        nt_fl.writerow(["n", "nt", "sr", "mv", "bm", "vs", "es", "ss"])
        nt_fl.writerow([net.settings['net_type'], net.settings['surf_res'], net.settings['max_vert'],
                        net.settings['box_size'], len(net.verts), len(net.edges), len(net.surfs)])

        # Write the connections header
        nt_fl.writerow(["c", "e0a0", "e0a1", "e0a2", "e1a0", "e1a1", "e1a2", "e2a0", "e2a1",
                        "e2a2", "e3a0", "e3a1", "e3a2", "s0a0", "s0a1", "s1a0", "s1a1", "s2a0", "s2a1", "s3a0", "s3a1",
                        "s4a0", "s4a1", "s5a0", "s5a1"])
        # Write the connections
        for i, vert in net.verts.iterrows():
            # Reset the tracking variables
            edge_ndxs, surf_ndxs = [], []
            # Stupid dumb way
            for j in range(4):
                if j >= len(vert['edges']):
                    edge_ndxs += [-1, -1, -1]
                else:
                    edge_ndxs += vert['edges'][j]
            for j in range(6):
                if j >= len(vert['surfs']):
                    surf_ndxs += [-1, -1]
                else:
                    surf_ndxs += vert['surfs'][j]
            # Write the vertex connection data
            nt_fl.writerow([i, *edge_ndxs, *surf_ndxs])

        # Create a vertices header
        nt_fl.writerow(["v", "a0", "a1", "a2", "a3", "x", "y", "z", "r"])
        # Write the connections and location and radius for each vertex in the network
        for i, vert in net.verts.iterrows():
            nt_fl.writerow([i, *vert['balls'], *r(vert['loc']), r(vert['rad'])])

        # Create an edges header
        nt_fl.writerow(["e", "a0", "a1", "a2", "sa0", "sa1", "i_0", "i_n"])
        # Write the connections and surface and points range information for each edge in the network
        for i, edge in net.edges.iterrows():
            # Write the edge information in the file
            nt_fl.writerow([i, *edge['balls'], *edge['ref']['surf'], edge['ref']['i0'], edge['ref']['i1']])

        # Create a surfaces header
        nt_fl.writerow(["s", "a0", "a1", "pts/tris"])
        # Write the connections and surface and points range information for each edge in the network
        for i, surf in net.surfs.iterrows():
            # Combine the points
            surf_points = []
            for point in surf['points']:
                surf_points += list(point)
            # Combine the tris
            surf_tris = []
            for tri in surf['tris']:
                surf_tris += tri
            # Write the surface points
            nt_fl.writerow(["pts", *surf['balls'], *[r(_) for _ in surf_points]])
            # Write the surface triangles
            nt_fl.writerow(["tris", *surf['balls'], *surf_tris])

    # Change back to the network file's directory
    os.chdir(group.sys.files['dir'])


def write_verts(net):
    """
    Exports a txt file with the vertex information for reloading later
    :param net: The network to interpret the vertex data from
    """
    # Move to the 181L output directory
    os.chdir(net.settings['sys_dir'])

    # Open the file for the vertices
    with open(net.settings['sys_dir'] + "/verts.txt", 'w') as file:
        # Create a header for the vertices file
        file.write("Vertices - {} vertices, {} atoms, max vert = {}, Net type = {}\n"
                   .format(len(net.verts['balls']), len(net.group), max(net.verts['rad']),
                           net.settings['net_type']))
        # Write the vertices
        for i, vert in net.verts.iterrows():
            # Write the vertex
            file.write(" ".join([str(_) for _ in vert['balls']]) + " " + " ".join([str(_) for _ in vert['loc']]) +
                       " " + str(vert['rad']) + " " + str(vert['dub']) + "\n")
        # Write the end line for the file
        file.write("END")


def add_metrics(group):
    net = group.net
    with open(group.sys.files['root_dir'] + '/Data/user_data/metrics.csv', 'a') as metrics_file:
        # name, # atoms, # verts, # edges, # surfs, # grp atoms, grp vol, grp sa, doublets, type, surf_res, max_vert, grp dsty
        metrics_file.write('\n{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}'
                           .format(group.sys.name,
                                   net.metrics['tot'],
                                   len(net.balls['num']),
                                   len(net.verts['balls']),
                                   len(net.edges['balls']),
                                   len(net.surfs['balls']),
                                   len(group.ball_ndxs),
                                   group.vol,
                                   group.sa,
                                   sum(net.verts['dub']),
                                   net.settings['net_type'],
                                   net.settings['surf_res'],
                                   net.settings['max_vert'],
                                   net.metrics['splits'],
                                   group.density))

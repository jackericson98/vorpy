import os
import csv
from vorpy.src.calculations import round_func


def write_net(group, file_name=None, round_to=3):
    """
    Exports a network checkpoint file containing the complete state of a Voronoi network.
    
    This function creates a CSV file that captures all essential network data including:
    - Network settings (type, resolution, vertex limits)
    - Vertex information (locations, radii, connections)
    - Edge data (connections, surface associations)
    - Surface information (connections, areas)
    
    The file can be used to restore the network state later or for analysis purposes.
    
    Args:
        group: Group object containing the network to export
        file_name (str, optional): Name of the output file. If None, generates a name based on the system.
        round_to (int, optional): Number of decimal places to round numerical values to. Defaults to 3.
        
    Returns:
        None: Creates a .csv file containing the network data
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
    with open(net.settings['sys_dir'] + "/{}_verts.txt".format(net.settings['net_type']), 'w') as file:
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
    """
    Adds metrics to the metrics file
    :param group: The group to add the metrics to
    """
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

import csv
from datetime import datetime
from vorpy.src.calculations import round_func


def write_logs(group, net_name=None, round_to=3):
    """
    Exports a comprehensive log file containing detailed information about the network analysis.
    
    This function generates a CSV log file with multiple sections:
    
    Build Information:
    - Network name, location, and completion date
    - Network type and key parameters (surface resolution, box size, max vertices)
    - Performance metrics (total time, vertex processing time, connection time, etc.)
    - Maximum vertex radius found
    
    Group Information:
    - Basic properties (name, volume, surface area, mass, density)
    - Center of mass (both standard and VDW)
    - Moment of inertia tensors (standard and spatial)
    
    Atoms:
    - Basic atom properties (index, name, residue info, chain, mass)
    - Spatial information (coordinates, radius, volumes)
    - Curvature metrics (mean and Gaussian curvatures)
    - Topological properties (sphericity, isometric quotient)
    - Neighbor analysis (count, distances, overlaps)
    - Contact areas and volumes
    - Center of mass and moment of inertia
    
    Args:
        group: Group object containing the network and system information
        net_name (str, optional): Additional identifier for the log file name
        round_to (int, optional): Number of decimal places to round numerical values to. Defaults to 3.
    """
    net = group.net
    # Set the network name for the logs file to be exported
    net_name = '' if net_name is None else '_' + str(net_name)
    # Create the round function
    r = round_func(round_to)
    # Open the file
    with open(group.settings['net_type'] + "_logs.csv", 'w') as log_file:
        # Create the csv writer
        lg_fl = csv.writer(log_file, lineterminator='\n')
        # Write the build information header
        lg_fl.writerow(["build informaiton"])
        # Write the build information labels
        lg_fl.writerow(["Name", "Location", "Completion Date", "Network Type", "Surface Resolution", "Box Size", "Maximum Allowable Vertex", "Total Time", "Vertex Time",
                       "Connect Time", "Surface Building Time", "Analysis time", "Maximum Found Vertex", ])
        lg_fl.writerow([group.sys.name, group.sys.files['base_file'], datetime.now(), net.settings['net_type'], net.settings['surf_res'], net.settings['box_size'],
                        net.settings['max_vert'], r(net.metrics['tot']), r(net.metrics['vert']), r(net.metrics['con']),
                        r(net.metrics['surf']), r(net.metrics['anal']), r(max(net.verts['rad']))])
        # Write the group information header
        lg_fl.writerow(["group information"])
        # Write the group information labels
        lg_fl.writerow(["Name", "Volume", "Surface Area", "Mass", "Density", "Center of Mass", "VDW Volume",
                        "VDW Center of Mass", "Moment of Inertia", 'Spatial Moment of Inertia'])
        # Write the group information
        group.get_info()
        lg_fl.writerow([group.name, r(group.vol), r(group.sa), float(group.mass), r(group.density), [float(r(_)) for _ in group.com],
                        r(group.vdw_vol), [float(r(_)) for _ in group.vdw_com], [[float(r(__)) for __ in _] for _ in group.moi],
                        [[float(r(__)) for __ in _] for _ in group.spatial_moment]])
        # Write the atom header
        lg_fl.writerow(["Atoms"])
        # Write the column labels
        lg_fl.writerow(["Index", "Name", "Residue", "Residue Sequence", "Chain", "Mass", "X", "Y", "Z", "Radius",
                        "Volume", "Van Der Waals Volume", "Surface Area", "Complete Cell?",
                        "Maximum Mean Curvature", "Average Mean Surface Curvature", "Maximum Gaussian Curvature",
                        "Average Gaussian Surface Curvature", "Sphericity", "Isometric Quotient",
                        "Inner Ball?", "Number of Neighbors", "Closest Neighbor", "Closest Neighbor Distance",
                        "Layer Distance Average", "Layer Distance RMSD", "Minimum Point Distance",
                        "Maximum Point Distance", "Number of Overlaps", "Contact Area", "Non-Overlap Volume",
                        "Overlap Volume", "Center of Mass", "Moment of Inertia Tensor", "Bounding Box", "neighbors"])
        # Go through the atoms in the system
        sys_balls = group.sys.balls.iloc[[_ for _ in net.balls['num']]].to_dict(orient='records')
        for i, atom in net.balls.iterrows():
            sys_ball = sys_balls[i]
            if atom['sa'] == 0:
                continue
            if atom['complete']:
                nbrs = [satoms[0] if satoms[0] != atom['num'] else satoms[1] for satoms in [net.surfs['balls'][_] for _ in atom['surfs']]]
                lg_fl.writerow([i, sys_ball['name'], sys_ball['res_name'], sys_ball['res_seq'], sys_ball['chain_name'],
                                sys_ball['mass'], atom['loc'][0], atom['loc'][1], atom['loc'][2], atom['rad'],
                                r(atom['vol']), r(atom['vdw_vol']), r(atom['sa']), atom['complete'],
                                r(atom['max_mean_curv']), r(atom['avg_mean_surf_curv']), r(atom['max_gauss_curv']),
                                r(atom['avg_gauss_surf_curv']), r(atom['sphericity']), r(atom['isometric_quotient']),
                                atom['ball_inside'], atom['number_of_neighbors'], atom['nearest_neighbor'],
                                atom['nearest_neighbor_distance'],
                                [float(_) for _ in r(atom['neighbor_distance_average'])],
                                [float(_) for _ in r(atom['neighbor_distance_rmsd'])], r(atom['min_spike']),
                                r(atom['max_spike']), atom['number_of_olaps'], r(atom['contact_area']),
                                r(atom['olap_vol']), r(atom['vdw_vol']), [float(r(_)) for _ in atom['com']],
                                [[float(r(__)) for __ in _] for _ in atom['moi']],
                                [[float(r(_)) for _ in atom['bounding_box'][0]],
                                 [float(r(_)) for _ in atom['bounding_box'][1]]], nbrs])
        # Write the surfaces header
        lg_fl.writerow(["Surfaces"])
        # Write the surface column labels
        lg_fl.writerow(["Index", "Ball 1", "Ball 2", "Surface Area", "Mean Curvature", "Gaussian Curvature",
                        "Ball 1 Volume Contribution", "Ball 2 Volume Contribution", 'Contact Area', 'Overlap'])
        # Go through the surfaces in the system and write their information
        for i, surf in net.surfs.iterrows():
            # Write the information for the surface
            lg_fl.writerow([i, *surf['balls'], r(surf['sa']), r(surf['mean_curv']), r(surf['gauss_curv']),
                            r(surf['vols'][surf['balls'][0]]), r(surf['vols'][surf['balls'][1]]),
                            r(surf['contact_area']), r(surf['overlap'])])
        # Write the edges header
        lg_fl.writerow(["Edges"])
        # Write the edges headers
        lg_fl.writerow(["Index", "Ball 1", "Ball 2", "Ball 3", "Length"])
        # Go through the edges in the network
        for i, edge in net.edges.iterrows():
            # Write the data for the edge
            lg_fl.writerow([i, *edge['balls'], r(edge['length'])])
        # Write the vertices header
        lg_fl.writerow(["Vertices"])
        # Write the vertices data labels
        lg_fl.writerow(["Index", "Ball 1", "Ball 2", "Ball 3", "Ball 4", "x", "y", "z", "r"])
        # Go through the vertices
        for i, vert in net.verts.iterrows():
            # Write the vertex information line
            lg_fl.writerow([i, *vert['balls'], *r(vert['loc']), r(vert['rad'])])

from pandas import DataFrame
from System.sys_funcs.calcs.calcs import ndx_search
import csv
import os.path
import numpy as np


def read_net(net, file_name):
    """
    Reads the network file for the import network and interprets the lines
    :param net: Network object for integration into
    :param file_name: Name of the read file
    :return: The network object
    """
    # Check the file_name
    if os.path.exists(file_name):
        file = file_name
    else:
        return
    # Set up the data lists
    verts, edges, surfs, cons = np.array([]), np.array([]), np.array([]), np.array([])
    # Open the file
    with open(file, 'r') as net_file:
        # Get the file element array to read
        nt_fl = csv.reader(net_file, delimiter=",")
        # Set the read type to header
        reading = ""
        # Go through the lines in the read_file
        for i, line in enumerate(nt_fl):
            # Read the first line
            if line[0] in {'n', 'v', 'e', 's', 'c'}:
                reading = line[0]
                continue
            if reading == 'net':
                net.id, net.type, net.surf_res, net.max_vert, net.box_size = \
                    [int(line[0]), line[1], float(line[2]), float(line[3]), float(line[4])]
                # Read the verts
            elif reading == 'v':
                # Add the data
                np.append(verts, line)
            # Read the edges
            elif reading == 'e':
                # Add the edge data
                np.append(edges, line)
            # Read the surfaces
            elif reading == 's':
                # Surface points
                if line[0] == 'pts':
                    np.append(surfs, ({"atoms": {*line[1:3]}}))
                    surfs[-1]['points'] = line[3:]
                # Triangles
                elif line[0] == 'tris':
                    surfs[-1]['tris'] = line[3:]
            # Read the connections
            elif reading == 'c':
                # Add the connections
                np.append(cons, line)
    integrate_net(net, verts, edges, surfs, cons)
    return net


def integrate_verts(net, verts):
    """
    Integrates the vertex object inputs
    :param net: The network to be integrated into
    :param verts: Vertices to add
    """
    if net.vertsa is None:
        net.verts = []
    # Set up the vertex ndxs
    if net.vert_ndxs is None:
        net.vert_ndxs = []
    # Go through the vertices
    for i, vert in enumerate(verts):
        # Get the index for the vertex
        ndx = [int(_) for _ in vert[1:5]]
        vert_ndx = ndx_search(net.vert_ndxs, ndx)
        # Check if the vertex exists
        if vert_ndx >= len(net.vert_ndxs) or net.vert_ndxs[vert_ndx] != ndx:
            # Create the vertex
            my_vert = make_vert(net=net, location=np.array([float(_) for _ in vert[5:8]]), radius=float(vert[8]),
                                ndx=ndx, atoms=np.array([net.atoms[_] for _ in ndx]))
            net.verts.insert(vert_ndx, my_vert)
            net.vert_ndxs.insert(vert_ndx, ndx)
            for j in ndx:
                np.append(net.atoms[j]['verts'], my_vert)
    # Make a dataframe out of the vertices
    net.verts = DataFrame(net.verts)


def integrate_edges(net, edges):
    """
    Integrates the edges
    :param net: Network to be integrated into
    :param edges: Edges to integrate
    """
    # Set up the network lists
    if net.edges is None:
        net.edges = []
    if net.edge_ndxs is None:
        net.edge_ndxs = []
    # Go through the edges in the network
    for i, edge in enumerate(edges):
        # Get the index for the surface
        ndx = [int(_) for _ in edge[1:4]]
        edge_ndx = ndx_search(net.edge_ndxs, ndx)
        # Check if the Edge exists
        if edge_ndx >= len(net.edge_ndxs) or net.edge_ndxs[edge_ndx] != ndx:
            # Create the Edge
            surf_ndx1 = ndx_search(net.surf_ndxs, [int(edge[4]), int(edge[5])])
            surf = net.surfs[surf_ndx1]
            atoms = [net.atoms[_] for _ in ndx]
            if int(edge[6]) != -1:
                points = surf.points[int(edge[6]):int(edge[7])]
            else:
                points = None
            ref = {'surf': [int(_) for _ in edge[4:6]], 'i0': int(edge[6]), 'i1': int(edge[7])}
            my_edge = make_edge(net=net, atoms=atoms, ndx=ndx, points=points, ref=ref)
            net.edges.insert(edge_ndx, my_edge)
            net.edge_ndxs.insert(edge_ndx, ndx)
            for j in ndx:
                net.atoms[j].edges.append(my_edge)
    # Set the DataFrame
    net.edges = DataFrame(net.edges)


def integrate_surfs(net, surfs):
    """
    Integrates the surfaces
    :param net: Network object to be integrated into
    :param surfs: Surfaces for integration
    """
    # Set up the network lists
    if net.surfs is None:
        net.surfs = []
    if net.surf_ndxs is None:
        net.surf_ndxs = []
    # Go through the surfaces
    for i, surf in enumerate(surfs):
        # Get the index for the surface
        ndx = [int(_) for _ in surf['atoms']]
        surf_ndx = ndx_search(net.surf_ndxs, ndx)
        # Check if the Surface exists
        if surf_ndx >= len(net.surf_ndxs) or net.surf_ndxs[surf_ndx] != ndx:
            # Create the Surface
            points = [[float(_) for _ in point] for point in [surf['points'][a:a+3] for a in range(0, len(surf['points']), 3)]]
            tris = [[int(_) for _ in tri] for tri in [surf['tris'][a:a+3] for a in range(0, len(surf['tris']), 3)]]
            my_surf = make_surf(net=net, atoms=[net.atoms[_] for _ in ndx], ndx=ndx, points=points, tris=tris,
                                resolution=net.surf_res)
            net.surfs.insert(surf_ndx, my_surf)
            net.surf_ndxs.insert(surf_ndx, ndx)
            # Add the atoms to the surface
            for j in ndx:
                net.atoms[j].surfs.append(my_surf)
    # Set up the dataframe
    net.surfs = DataFrame(net.surfs)


def integrate_net(net, verts, edges, surfs, cons):
    """
    Integrates the new network file into the network
    :param net: Main network object to integrate the objects into
    :param verts: Vertex list
    :param edges: Edge list
    :param surfs: Surface list
    :param cons: Connection list
    """
    # Integrate verts, edges, surfs
    integrate_verts(net, verts)
    integrate_surfs(net, surfs)
    integrate_edges(net, edges)
    # Go through the vertices and interpret everything
    for i, vcon in enumerate(cons):
        # Get the atoms
        vert_ndx = ndx_search(net.vert_ndxs, np.array([int(_) for _ in verts[i][1:5]]))
        vert = net.verts[vert_ndx]
        # Get the edge indices
        my_edges = [vcon[a:a+3] for a in range(1, 11, 3)]
        v_edge_atoms = [[int(_) for _ in edge] for edge in my_edges if int(edge[0]) != -1]
        vert['edges'] = [net.edges[ndx_search(net.edge_ndxs, _)] for _ in v_edge_atoms]
        # Get the surface indices
        my_surfs = [vcon[a:a+2] for a in range(13, 24, 2)]
        v_surf_atoms = [[int(_) for _ in surf] for surf in my_surfs if int(surf[0]) != -1]
        vert.surfs = [net.surfs[ndx_search(net.surf_ndxs, _)] for _ in v_surf_atoms]
        # Add the edges and surfaces together
        for edge in vert.edges:
            if edge.verts is None:
                edge.verts = []
            edge.verts.append(vert)
            for surf in vert.surfs:
                if len([_ for _ in surf.ndx if _ in edge.ndx]) == 2:
                    if edge.surfs is None:
                        edge.surfs = []
                    edge.surfs.append(surf)
                    if surf.edges is None:
                        surf.edges = []
                    surf.edges.append(edge)
    net.metrics = {'tot': 0, 'vert': 0, 'con': 0, 'surf': 0, 'anal': 0}
    net.analyze()

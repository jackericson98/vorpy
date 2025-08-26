import numpy as np
from pandas import DataFrame
from vorpy.src.calculations import ndx_search


# def integrate_verts(net, verts):
#     """
#     Integrates vertex data into a network structure.

#     This function processes vertex objects and adds them to the network's vertex list.
#     It handles vertex creation, indexing, and linking with associated atoms.

#     Parameters:
#     -----------
#     net : Network
#         The network object to integrate vertices into
#     verts : list
#         List of vertex data to be integrated, where each vertex contains:
#         - Index information (4 values)
#         - Location coordinates (3 values)
#         - Radius value

#     Returns:
#     --------
#     None
#         Modifies the network object in place by:
#         - Adding new vertices to net.verts
#         - Updating vertex indices in net.vert_ndxs
#         - Linking vertices to associated atoms
#     """
#     # Initialize the vertex list
#     if net.verts is None:
#         net.verts = []
#     # Initialize the vertex index list
#     if net.vert_ndxs is None:
#         net.vert_ndxs = []
#     # Go through the vertices
#     for i, vert in enumerate(verts):
#         # Get the index for the vertex
#         ndx = [int(_) for _ in vert[1:5]]
#         vert_ndx = ndx_search(net.vert_ndxs, ndx)
#         # Check if the vertex exists
#         if vert_ndx >= len(net.vert_ndxs) or net.vert_ndxs[vert_ndx] != ndx:
#             # Create the vertex
#             my_vert = make_vert(net=net, location=np.array([float(_) for _ in vert[5:8]]), radius=float(vert[8]),
#                                 ndx=ndx, atoms=np.array([net.atoms[_] for _ in ndx]))
#             # Insert the vertex
#             net.verts.insert(vert_ndx, my_vert)
#             # Insert the vertex index
#             net.vert_ndxs.insert(vert_ndx, ndx)
#             # Add the vertex to the atoms
#             for j in ndx:
#                 # Check if the atom is already on the vertex
#                 net.atoms[j].verts.append(my_vert)
#     # Make a dataframe out of the vertices
#     net.verts = DataFrame(net.verts)


# def integrate_edges(net, edges):
#     """
#     Integrates edge data into a network structure.

#     This function processes edge objects and adds them to the network's edge list.
#     It handles edge creation, indexing, and linking with associated atoms and surfaces.

#     Parameters:
#     -----------
#     net : Network
#         The network object to integrate edges into
#     edges : list
#         List of edge data to be integrated, where each edge contains:
#         - Index information (3 values)
#         - Surface reference information
#         - Point indices for surface points

#     Returns:
#     --------
#     None
#         Modifies the network object in place by:
#         - Adding new edges to net.edges
#         - Updating edge indices in net.edge_ndxs
#         - Linking edges to associated atoms and surfaces
#     """
#     # Set up the network lists
#     if net.edges is None:
#         # Initialize the edge list
#         net.edges = []
#     if net.edge_ndxs is None:
#         # Initialize the edge index list
#         net.edge_ndxs = []
#     # Go through the edges in the network
#     for i, edge in enumerate(edges):
#         # Get the index for the surface
#         ndx = [int(_) for _ in edge[1:4]]
#         # Get the index for the edge
#         edge_ndx = ndx_search(net.edge_ndxs, ndx)
#         # Check if the Edge exists
#         if edge_ndx >= len(net.edge_ndxs) or net.edge_ndxs[edge_ndx] != ndx:
#             # Create the Edge
#             surf_ndx1 = ndx_search(net.surf_ndxs, [int(edge[4]), int(edge[5])])
#             # Get the surface
#             surf = net.surfs[surf_ndx1]
#             # Get the atoms
#             atoms = [net.atoms[_] for _ in ndx]
#             # Get the points
#             if int(edge[6]) != -1:
#                 points = surf.points[int(edge[6]):int(edge[7])]
#             else:
#                 points = None
#             # Get the reference
#             ref = {'surf': [int(_) for _ in edge[4:6]], 'i0': int(edge[6]), 'i1': int(edge[7])}
#             # Make the edge
#             my_edge = make_edge(net=net, atoms=atoms, ndx=ndx, points=points, ref=ref)
#             # Insert the edge
#             net.edges.insert(edge_ndx, my_edge)
#             # Insert the edge index
#             net.edge_ndxs.insert(edge_ndx, ndx)
#             # Add the edge to the atoms
#             for j in ndx:
#                 # Check if the atom is already on the edge
#                 net.atoms[j].edges.append(my_edge)
#     # Set the DataFrame
#     net.edges = DataFrame(net.edges)


# def integrate_surfs(net, surfs):
#     """
#     Integrates surface data into a network object.

#     Parameters:
#     -----------
#     net : Network
#         The network object to integrate surfaces into
#     surfs : list
#         List of surface data to be integrated, where each surface contains:
#         - Atom indices defining the surface
#         - Point coordinates for surface vertices
#         - Triangle indices for surface triangulation

#     Returns:
#     --------
#     None
#         Modifies the network object in place by:
#         - Adding new surfaces to net.surfs
#         - Updating surface indices in net.surf_ndxs
#         - Linking surfaces to associated atoms
#     """
#     # Set up the network lists
#     if net.surfs is None:
#         # Initialize the surface list
#         net.surfs = []
#     if net.surf_ndxs is None:
#         # Initialize the surface index list
#         net.surf_ndxs = []
#     # Go through the surfaces
#     for i, surf in enumerate(surfs):
#         # Get the index for the surface
#         ndx = [int(_) for _ in surf['atoms']]
#         surf_ndx = ndx_search(net.surf_ndxs, ndx)
#         # Check if the Surface exists
#         if surf_ndx >= len(net.surf_ndxs) or net.surf_ndxs[surf_ndx] != ndx:
#             # Create the Surface
#             points = [[float(_) for _ in point] for point in [surf['points'][a:a+3] for a in range(0, len(surf['points']), 3)]]
#             # Get the triangles
#             tris = [[int(_) for _ in tri] for tri in [surf['tris'][a:a+3] for a in range(0, len(surf['tris']), 3)]]
#             # Make the surface
#             my_surf = make_surf(net=net, atoms=[net.atoms[_] for _ in ndx], ndx=ndx, points=points, tris=tris,
#                                 resolution=net.settings['surf_res'])
#             # Insert the surface
#             net.surfs.insert(surf_ndx, my_surf)
#             # Insert the surface index
#             net.surf_ndxs.insert(surf_ndx, ndx)
#             # Add the atoms to the surface
#             for j in ndx:
#                 # Check if the atom is already on the surface
#                 net.atoms[j].surfs.append(my_surf)
#     # Set up the dataframe
#     net.surfs = DataFrame(net.surfs)


# def integrate_net(net, verts, edges, surfs, cons):
#     """
#     Integrates network data into a network object.

#     Parameters:
#     -----------
#     net : Network
#         The network object to integrate data into
#     verts : list
#         List of vertex data containing atom indices and coordinates
#     edges : list
#         List of edge data defining connections between vertices
#     surfs : list
#         List of surface data containing atom indices and triangulation
#     cons : list
#         List of connection data linking vertices to edges and surfaces

#     Returns:
#     --------
#     None
#         Modifies the network object in place by:
#         - Adding vertices, edges, and surfaces
#         - Establishing connections between network components
#         - Updating network metrics and analysis
#     """
#     # Integrate verts, edges, surfs
#     integrate_verts(net, verts)
#     integrate_surfs(net, surfs)
#     integrate_edges(net, edges)
#     # Go through the vertices and interpret everything
#     for i, vcon in enumerate(cons):
#         # Get the atoms
#         vert_ndx = ndx_search(net.vert_ndxs, np.array([int(_) for _ in verts[i][1:5]]))
#         vert = net.verts[vert_ndx]
#         # Get the edge indices
#         my_edges = [vcon[a:a+3] for a in range(1, 11, 3)]
#         # Get the edge atoms
#         v_edge_atoms = [[int(_) for _ in edge] for edge in my_edges if int(edge[0]) != -1]
#         # Get the edge objects
#         vert['edges'] = [net.edges[ndx_search(net.edge_ndxs, _)] for _ in v_edge_atoms]
#         # Get the surface indices
#         my_surfs = [vcon[a:a+2] for a in range(13, 24, 2)]
#         # Get the surface atoms
#         v_surf_atoms = [[int(_) for _ in surf] for surf in my_surfs if int(surf[0]) != -1]
#         # Get the surface objects
#         vert.surfs = [net.surfs[ndx_search(net.surf_ndxs, _)] for _ in v_surf_atoms]
#         # Add the edges and surfaces together
#         for edge in vert.edges:
#             # Add the vertex to the edge
#             if edge.verts is None:
#                 edge.verts = []
#             edge.verts.append(vert)
#             # Add the edge to the surfaces
#             for surf in vert.surfs:
#                 # Check if the edge is on the surface
#                 if len([_ for _ in surf.ndx if _ in edge.ndx]) == 2:
#                     # Add the edge to the surface
#                     if edge.surfs is None:
#                         edge.surfs = []
#                     edge.surfs.append(surf)
#                     # Add the surface to the edge
#                     if surf.edges is None:
#                         surf.edges = []
#                     surf.edges.append(edge)
#     # Update the network metrics
#     net.metrics = {'tot': 0, 'vert': 0, 'con': 0, 'surf': 0, 'anal': 0}
#     # Analyze the network
#     net.analyze()

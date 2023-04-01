import csv
import os.path

from System.Network.network import Network
from System.Network.net_objs.vertex import Vertex
from System.Network.net_objs.edge import Edge
from System.Network.net_objs.surface import Surface


def integrate_net(net, new_objs):
    # Set up the lists
    verts, edges, surfs = new_objs
    ndx = 0
    # Go through the vertices
    for vert in verts:
        # Keep checking the indices of the vertices in the network until we find one larger
        while net.vert_ndxs[ndx] < vert.ndx:
            ndx += 1
        # Check if the ndxs are the same
        if net.vert_ndxs[ndx] == vert.ndx:
            net.verts[ndx].add_vert_info(vert)
        # Insert the vertex
        net.verts.insert(ndx, vert)
        net.vert_ndxs.insert(ndx, vert.ndx)
    ndx = 0
    # Go through the edges
    for edge in edges:
        # Keep checking the indices of the vertices in the network until we find one larger
        while net.edge_ndxs[ndx] < edge.ndx:
            ndx += 1
        # Check if the ndxs are the same
        if net.edge_ndxs[ndx] == edge.ndx:
            net.edge[ndx].add_edge_info(edge)
        # Insert the vertex
        net.edges.insert(ndx, edge)
        net.edge_ndxs.insert(ndx, edge.ndx)
    ndx = 0
    # Go through the surfaces
    for surf in surfs:
        # Keep checking the indices of the surfaces in the network until we find one larger
        while net.surf_ndxs[ndx] < surf.ndx:
            ndx += 1
        # Check if the ndxs are the same
        if net.surf_ndxs[ndx] == surf.ndx:
            net.surf[ndx].add_surf_info(surf)
        # Insert the vertex
        net.surfs.insert(ndx, surf)
        net.surf_ndxs.insert(ndx, surf.ndx)



def read_net(sys, file=None, integrate=False):
    # Open the file
    if file is None:
        file = sys.net_file
        if sys.net_file is None:
            return
    # Get the directory for the surfaces
    net_dir = os.path.dirname(file)
    # Keep using the same directory, this will cut down on clutter
    sys.dir = net_dir
    # Open the file
    with open(file, 'r') as my_file:
        # Get the file element array to read
        read_file = list(csv.reader(my_file, delimiter=","))
        # Get the network information
        net_verts, net_edges, net_surfs = [int(_) for _ in read_file[1][5:8]]
        # Create the network if needed
        if sys.net is None:
            sys.net = Network(sys=sys, atoms=sys.atoms)
        # Create the blank objects
        my_verts = [Vertex(net=sys.net) for _ in range(net_verts)]
        my_edges = [Edge(net=sys.net) for _ in range(net_edges)]
        my_surfs = [Surface(net=sys.net) for _ in range(net_surfs)]
        # Add the settings
        surf_res, max_vert, box_size = [float(_) for _ in read_file[1][1:4]]
        build_surfs = bool(read_file[1][4])
        sys.net.calc_box()
        # Add the vertices
        for i in range(net_verts):
            print("\rloading vertices - {}%".format(round(100 * i / net_verts, 2)), end="")
            vert = my_verts[i]
            line = read_file[i + 3]
            vert.loc = [float(_) for _ in line[1:4]]
            vert.rad = float(line[4])
            vert.atoms = [sys.atoms[int(_)] for _ in line[5:9]]
            vert.ndx = [int(_) for _ in line[5:9]]
            vert.edges = [my_edges[int(_)] for _ in line[9:14] if _ != '']
            surf_ndxs = [int(_) for _ in line[14:] if _ != '']
            vert.surfs = [my_surfs[_] for _ in surf_ndxs]
            if i >= 1 and vert.ndx == my_verts[i - 1].ndx:
                vert.doublet = my_verts[i - 4]
                my_verts[i - 4].doublet = vert
            for atom in vert.atoms:
                atom.verts.append(vert)
            for surf in vert.surfs:
                if surf.verts is None:
                    surf.verts = []
                surf.verts.append(vert)
        # Add the edges
        for i in range(net_edges):
            print("\rloading edges - {}%".format(round(100 * i / net_edges, 2)), end="")
            edge = my_edges[i]
            line = read_file[i + 4 + net_verts]
            edge.atoms = [sys.atoms[int(_)] for _ in line[4:7]]
            edge.ndx = [int(_) for _ in line[4:7]]
            edge.verts = [my_verts[int(_)] for _ in line[7:9]]
            edge.surfs = [my_surfs[int(_)] for _ in line[9:] if _ != '']
            for atom in edge.atoms:
                atom.edges.append(edge)
            for surf in edge.surfs:
                if surf.edges is None:
                    surf.edges = []
                surf.edges.append(edge)
        # Add the surfaces
        # noinspection PyTypeChecker
        for i in range(net_surfs):
            print("\rloading surfaces - {}%".format(round(100 * i / net_surfs, 2)), end="")
            surf = my_surfs[i]
            line = read_file[i + 5 + net_verts + net_edges]
            surf.atoms = [sys.atoms[int(_)] for _ in line[5:7]]
            if surf.atoms[0].rad > surf.atoms[1].rad:
                surf.atoms[0], surf.atoms[1] = surf.atoms[1], surf.atoms[0]
            surf.ndx = [int(_) for _ in line[5:7]]
            if line[1] != '':
                surf.file = line[1]
            if line[2] != '':
                surf.res = float(line[2])
            if line[3] != '':
                surf.sa = float(line[3])
            if line[4].isdigit():
                surf.curv = float(line[4])
            if isinstance(line[16], tuple):
                surf.func = [float(_) for _ in line[7:16]] + [float(_) for _ in line[16:]]
            else:
                surf.func = [float(_) for _ in line[7:]]
            for atom in surf.atoms:
                atom.surfs.append(surf)
        if surf.atoms[0].rad > surf.atoms[1].rad:
            surf.atoms[0], surf.atoms[1] = surf.atoms[1], surf.atoms[0]
    # Set the network to connected
    sys.net.connect_net = False

    if integrate:
        # Integrate the network and the old network
        integrate_net(sys.net, [my_verts, my_edges, my_verts])
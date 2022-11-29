

def check_net_nums(net):
    net_good = True
    vert_atoms, vert_edges, vert_surfs = True, True, True
    for vert in net.verts:
        bad_vert = False
        if len(vert.atoms) > 4:
            bad_vert = True
            vert_atoms = False
        if len(vert.edges) > 4:
            bad_vert = True
            vert_edges = False
        if len(vert.surfs) > 6:
            bad_vert = True
            vert_surfs = False

        if bad_vert:
            net_good = False
            print("Bad Vertex: {}. Corrct number of:\n Atoms: {}\nEdges: {}\nSurfaces: {}"
                  .format(vert.ndx, vert_atoms, vert_edges, vert_surfs))


    edge_atoms, edge_verts, edge_surfs = True, True, True
    for edge in net.edges:
        bad_edge = False
        if len(edge.atoms) > 3:
            bad_edge = True
            edge_atoms = False
        if len(edge.verts) > 2:
            bad_edge = True
            edge_verts = False
        if len(edge.surfs) > 3:
            bad_edge = True
            edge_surfs = False

        if bad_edge:
            net_good = False
            print("Bad Edge: {}. Corrct number of:\n Atoms: {}\nVertices: {}\nSurfaces: {}"
                  .format(edge.ndx, edge_atoms, edge_verts, edge_surfs))

    surf_atoms, surf_verts, surf_edges = True, True, True
    for surf in net.surfs:
        bad_surf = False
        if len(surf.atoms) > 2:
            bad_surf = True
            surf_atoms = False
        if len(surf.verts) < 3:
            bad_surf = True
            surf_verts = False
        if len(surf.surfs) < 2:
            bad_surf = True
            surf_edges = False

        if bad_surf:
            net_good = False
            print("Bad Surface: {}. Corrct number of:\n Atoms: {}\nVertices: {}\nEdges: {}"
                  .format(surf.ndx, surf_atoms, surf_verts, surf_edges))

    if net_good:
        print("Network good")


def check_surf_cons(surf):
    """
    Checks the edges and vertices of the surface to make sure everything is in order
    :param surf:
    :return:
    """
    verts, edges = surf.verts.copy(), surf.edges.copy()
    vert = verts[0]
    next_vert = None
    while len(verts) > 0 or len(edges) > 0:
        next_edge = None
        for i in range(len(edges)):
            if len([0 for _ in edges[i].ndx if _ in vert.ndx]) == 3:
                next_edge = edges.pop(i)
                continue

        if next_edge is None:
            return False

        next_vert = None
        for j in range(len(verts)):
            if len([0 for _ in next_edge.ndx if _ in verts[j].ndx]) == 3:
                next_vert = verts.pop(j)
                continue

        if next_vert is None:
            return False

    if next_vert != vert:
        return False
    return True

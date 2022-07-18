from Meshes.fill_mesh import *

########################################################################################################################


# Build meshes function. Runs make_mesh on all surfaces in the network
def build_meshes(sys, min_dist=None):
    # Set the minimum distance
    if min_dist is None:
        min_dist = calc_dist(sys.net.edges[0].verts[0].loc, sys.net.edges[0].verts[1].loc) / 30
    # Make each surface
    for surf in sys.net.surfs:
        make_mesh(surf, min_dist)
    # Return the surfaces
    return sys.net.surfs

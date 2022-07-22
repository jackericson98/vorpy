from Meshes.fill_mesh import *

########################################################################################################################


# Build meshes function. Runs make_mesh on all surfaces in the network
def build_meshes(sys=None, meshes=None, min_dist=None):
    if sys:
        mesh_list = sys.net.surfs
    elif meshes:
        mesh_list = meshes
    else:
        return
    # Set the minimum distance
    if min_dist is None:
        min_dist = 0.5
    # Make each surface
    for surf in sys.net.surfs:
        make_mesh(surf, min_dist, vta=sys.Voronota)
    # Return the surfaces
    return sys.net.surfs

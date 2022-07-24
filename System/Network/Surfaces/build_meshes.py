from System.Network.Surfaces.fill_mesh import *

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
    num_surfs = len(sys.net.surfs)
    # Make each surface
    for i in range(num_surfs):
        # Calculate and print the running percentage for mesh calculations
        percentage = int(i/num_surfs * 100) + 1
        print("\rBuilding Surfaces: ", percentage, "%", end='')
        make_mesh(sys.net.surfs[i], min_dist, vta=sys.Voronota)
    # Return the surfaces
    return sys.net.surfs

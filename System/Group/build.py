import os
from System.sys_funcs.output.output import write_surfs
from System.sys_funcs.input.input import read_surf_file
from System.Network.net_funcs.build_surf import build_surf


def build_surfs(grp, resolution=None):
    """
    Checks the surfaces for points and allows for rebuilds surfaces with incorrect resolutions
    :param grp:
    :param resolution: If not None all surfs without this resolution will be rebuilt
    :return: All surfaces in the group will be constructed
    """
    # Get the resolution
    if resolution is None:
        resolution = grp.sys.net.surf_res
        grp.surf_res = resolution
    # Set up the build surfaces list
    bld_surfs = []
    # Go through the list of build surfaces checking for
    for surf in grp.surfs:
        # Check if the resolution is different from the set resolution or the surface has no points
        if surf.res != resolution:
            bld_surfs.append(surf)
        # Check if there is any sign of missing points or triangles
        elif surf.points is None or surf.tris is None or len(surf.points) <= 2 or len(surf.tris) == 0:
            # If it is possible to load the file
            if surf.file is not None and surf.file not in ["", " "]:
                test = read_surf_file(surf)
                if test is None:
                    bld_surfs.append(surf)
            # Worst case, add the surface to the list of surfaces to be built
            else:
                bld_surfs.append(surf)
    # Create the system's surface's file if needed
    if len(bld_surfs) > 0 and not os.path.exists(grp.sys.dir + "/surfs"):
        os.mkdir(grp.sys.dir + "/surfs")
        os.chdir(grp.sys.dir + '/surfs')
    # Build the surfaces
    for i in range(len(bld_surfs)):
        # Print the status of the surfaces being built
        print("\rbuilding " + grp.name + " surfaces " + " " * (len(str(len(grp.surfs) - 1)) - len(str(i + 1))) +
              str(i + 1) + "/" + str(len(grp.surfs)) + "                   ", end="")
        build_surf(bld_surfs[i], res=resolution)
        if bld_surfs[i].file is None:
            write_surfs([bld_surfs[i]], "_".join([str(_) for _ in bld_surfs[i].ndx]))
    # Change back
    os.chdir(grp.sys.dir)

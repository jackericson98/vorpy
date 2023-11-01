import os
from System.sys_funcs.input.input import read_surf_file
from System.Network.surfs.build_surf import build_surf


def build_surfs(grp, resolution=None):
    """
    Checks the surfaces for points and allows for rebuilds surfaces with incorrect resolutions
    :param grp: Group object to build surfs from
    :param resolution: If not None all surfs without this resolution will be rebuilt
    :return: All surfaces in the group will be constructed
    """
    # Get the resolution
    if resolution is not None:
        grp.surf_res = resolution
    else:
        grp.surf_res = grp.sys.net.surf_res
    # Go through the list of build surfaces checking for
    for i in grp.surfs:
        surf = grp.sys.net.surfs.iloc[i]
        # Print the status of the surfaces being built
        print("\rbuilding " + grp.name + " surfaces " + " " * (len(str(len(grp.surfs) - 1)) - len(str(int(i) + 1))) +
              str(i + 1) + "/" + str(len(grp.surfs)) + "                   ", end="")
        # Check if there is any sign of missing points or triangles
        if surf['points'] is None or surf['tris'] is None or len(surf['points']) <= 2 or len(surf['tris']) == 0:
            spoints, surf_tris, tri_curvs, surf_curv, sfunc, surf_com, flat = build_surf(
                alocs=[grp.sys.net.atoms['loc'][_] for _ in surf['satoms']],
                arads=[grp.sys.net.atoms['rad'][_] for _ in surf['satoms']],
                epnts=[grp.sys.net.edges['points'][_] for _ in surf['sedges']], res=grp.surf_res,
                net_type=grp.sys.net.type)
            # Set the value in sht dataframe
            try:
                grp.sys.net.surfs.loc[i:i, ['points', 'tris', 'tri_curvs', 'curv', 'func', 'com', 'flat']] = spoints, surf_tris, tri_curvs, surf_curv, sfunc, surf_com, flat
            except ValueError:
                pass
    # Change back
    if grp.sys.dir is not None:
        os.chdir(grp.sys.dir)

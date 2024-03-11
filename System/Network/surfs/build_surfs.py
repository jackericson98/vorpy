import time
from System.sys_funcs.calcs.calcs import get_time
from System.Network.surfs.build_surf import build_surf
from System.sys_funcs.calcs.surf import calc_surf_sa
from System.sys_funcs.calcs.calcs import calc_tetra_vol


def build_surfs(net):
    # Make each surface
    points, tris, tri_curvs, curvs, funcs, coms, flats, sas, vols = [], [], [], [], [], [], [], [], []
    for i, surf in net.surfs.iterrows():
        # Build the surfaces and print the progress
        my_time = time.perf_counter() - net.start_time
        h, m, s = get_time(my_time)
        print("\rRun Time = {:2}:{:2}:{:.2f} - Process: building surfaces {:.2f} %                                 "
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * round(i / len(net.surfs), 4))), end="")
        arads = [net.atoms['rad'][_] for _ in surf['satoms']]
        alocs = [net.atoms['loc'][_] for _ in surf['satoms']]
        anums = [net.atoms['num'][_] for _ in surf['satoms']]
        if arads[0] > arads[1]:
            arads, alocs, anums = [arads[1], arads[0]], [alocs[1], alocs[0]], [anums[1], anums[0]]
        my_surf = build_surf(alocs=alocs, arads=arads, epnts=[net.edges['points'][_] for _ in surf['sedges']],
                             res=net.surf_res, net_type=net.type)
        surf_points, surf_tris, surf_tri_curvs, surf_curv, surf_func, surf_com, surf_flat = my_surf
        sas.append(calc_surf_sa(edges=[net.edges['points'][_] for _ in surf['sedges']], com=surf_com, tris=surf_tris,
                                points=surf_points, flat=surf_flat))
        points.append(surf_points)
        tris.append(surf_tris)
        tri_curvs.append(surf_tri_curvs)
        curvs.append(surf_curv)
        funcs.append(surf_func)
        coms.append(surf_com)
        flats.append(surf_flat)
        # Get the surface Volumes
        sv0 = sum([calc_tetra_vol(alocs[0], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])
        sv1 = sum([calc_tetra_vol(alocs[1], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])
        vols.append({anums[0]: sv0, anums[1]: sv1})
    # Set the dataframe elements
    net.surfs['points'], net.surfs['tris'], net.surfs['tri_curvs'], net.surfs['curv'], net.surfs['func'], \
     net.surfs['com'], net.surfs['flat'], net.surfs['sa'], net.surfs['vols'] = \
        points, tris, tri_curvs, curvs, funcs, coms, flats, sas, vols
    # Get the curvature in the 95th percentile
    my_surf_curvs = net.surfs['curv'].to_list()
    my_surf_curvs.sort()
    try:
        net.max_curv = my_surf_curvs[min(int(0.99 * len(my_surf_curvs)), len(my_surf_curvs) - 1)]
    except IndexError:
        net.max_curv = 0
    print("\r                                                                                             ", end='')
    net.metrics['surf'] = time.perf_counter() - net.start_time - net.metrics['vert'] - net.metrics['con']
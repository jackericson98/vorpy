import time
from System.sys_funcs.calcs.calcs import get_time
from System.Network.surfs.build_surf import build_surf
from System.sys_funcs.calcs.surf import calc_surf_sa, calc_surf_tri_curvs


def build_surfs(net):
    # Make each surface
    points, tris, tri_curvs, curvs, funcs, coms, flats, sas = [], [], [], [], [], [], [], []
    for i, surf in net.surfs.iterrows():
        # Build the surfaces and print the progress
        my_time = time.perf_counter() - net.start_time
        h, m, s = get_time(my_time)
        print("\rRun Time = {:2}:{:2}:{:.2f} - Process: building surfaces {:.2f} %                                 "
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * round(i / len(net.surfs), 4))), end="")
        arads = [net.atoms['rad'][_] for _ in surf['satoms']]
        alocs = [net.atoms['loc'][_] for _ in surf['satoms']]
        if arads[0] > arads[1]:
            arads, alocs = [arads[1], arads[0]], [alocs[1], alocs[0]]
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
    # Set the dataframe elements
    net.surfs['points'], net.surfs['tris'], net.surfs['tri_curvs'], net.surfs['curv'], net.surfs['func'], \
     net.surfs['com'], net.surfs['flat'], net.surfs['sa'] = points, tris, tri_curvs, curvs, funcs, coms, flats, sas
    print("\r                                                                                             ", end='')
    net.metrics['surf'] = time.perf_counter() - net.start_time - net.metrics['vert'] - net.metrics['con']
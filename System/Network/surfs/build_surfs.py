import time
import pandas as pd
from io import StringIO
from csv import writer
from System.sys_funcs.calcs.calcs import get_time
from System.Network.surfs.build_surf import build_surf
from System.sys_funcs.calcs.surf import calc_surf_sa
from System.sys_funcs.calcs.calcs import calc_tetra_vol


# Putting a pause on this but should be much faster
def build_surfs1(net):
    output = StringIO()
    csv_writer = writer(output)
    # Make each surface
    for i, surf in net.surfs.iterrows():
        # Build the surfaces and print the progress
        my_time = time.perf_counter() - net.metrics['start']
        h, m, s = get_time(my_time)
        print("\rRun Time = {:2}:{:2}:{:.2f} - Process: building surfaces {:.2f} %                                 "
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * round(i / len(net.surfs), 4))), end="")
        arads = [net.balls['rad'][_] for _ in surf['balls']]
        alocs = [net.balls['loc'][_] for _ in surf['balls']]
        anums = [net.balls['num'][_] for _ in surf['balls']]
        if arads[0] > arads[1]:
            arads, alocs, anums = [arads[1], arads[0]], [alocs[1], alocs[0]], [anums[1], anums[0]]
        my_surf = build_surf(locs=alocs, rads=arads, epnts=[net.edges['points'][_] for _ in surf['edges']],
                             res=net.settings['surf_res'], net_type=net.settings['net_type'])
        surf_points, surf_tris, surf_tri_curvs, surf_curv, surf_func, surf_com, surf_flat = my_surf
        # Get the surface Volumes
        sv0 = sum([calc_tetra_vol(alocs[0], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])
        sv1 = sum([calc_tetra_vol(alocs[1], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])

        sa = calc_surf_sa(edges=[net.edges['points'][_] for _ in surf['edges']], com=surf_com, tris=surf_tris,
                          points=surf_points, flat=surf_flat)
        csv_writer.writerow([surf_points, surf_tris, surf_tri_curvs, surf_curv, surf_func, surf_com, surf_flat, sa,
                             {anums[0]: sv0, anums[1]: sv1}])

    output.seek(0)
    # net.surfs = pd.read_csv(output, names=['points', 'tris', 'tri_curvs', 'curv', 'func', 'com', 'flat', 'sa', 'vols'],
    #                         dtype={'points': pd.array(), 'tris'})
    # Set the dataframe elements
    # Get the curvature in the 95th percentile
    my_surf_curvs = net.surfs['curv'].to_list()
    my_surf_curvs.sort()
    try:
        net.max_curv = my_surf_curvs[min(int(0.99 * len(my_surf_curvs)), len(my_surf_curvs) - 1)]
    except IndexError:
        net.max_curv = 0
    print("\r                                                                                             ", end='')
    net.metrics['surf'] = time.perf_counter() - net.metrics['start'] - net.metrics['vert'] - net.metrics['con']


def build_surfs(net, store_points=True):
    # Instantiate the lists for storage
    points, tris, tri_curvs, curvs, funcs, coms, flats, sas, vols = [], [], [], [], [], [], [], [], []
    # Make each surface
    for i, surf in net.surfs.iterrows():
        # Build the surfaces and print the progress
        my_time = time.perf_counter() - net.metrics['start']
        h, m, s = get_time(my_time)
        print("\rRun Time = {:2}:{:2}:{:.2f} - Process: building surfaces {:.2f} %                                 "
              .format(int(h), int(m), round(s, 2), min(100.0, 100 * round(i / len(net.surfs), 4))), end="")
        rads = [net.balls['rad'][_] for _ in surf['balls']]
        locs = [net.balls['loc'][_] for _ in surf['balls']]
        nums = [net.balls['num'][_] for _ in surf['balls']]
        if rads[0] > rads[1]:
            rads, locs, nums = [rads[1], rads[0]], [locs[1], locs[0]], [nums[1], nums[0]]
        my_surf = build_surf(locs=locs, rads=rads, epnts=[net.edges['points'][_] for _ in surf['edges']],
                             res=net.settings['surf_res'], net_type=net.settings['net_type'])
        surf_points, surf_tris, surf_tri_curvs, surf_curv, surf_func, surf_com, surf_flat = my_surf
        # Get the surface Volumes
        sv0 = sum([calc_tetra_vol(locs[0], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])
        sv1 = sum([calc_tetra_vol(locs[1], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]]) for tri in
                   surf_tris])
        # Calculate the surface area of the surface
        sa = calc_surf_sa(edges=[net.edges['points'][_] for _ in surf['edges']], com=surf_com, tris=surf_tris,
                          points=surf_points, flat=surf_flat)
        # If we are doing a large export and will need the points later in the process for export and such
        if store_points:
            points.append(surf_points)
            tris.append(surf_tris)
            tri_curvs.append(surf_tri_curvs)
        else:
            points.append([])
            tris.append([])
            tri_curvs.append([])
        curvs.append(surf_curv)
        funcs.append(surf_func)
        coms.append(surf_com)
        flats.append(surf_flat)
        sas.append(sa)
        vols.append({nums[0]: sv0, nums[1]: sv1})
    (net.surfs['points'], net.surfs['tris'], net.surfs['tri_curvs'], net.surfs['curv'], net.surfs['func'],
     net.surfs['com'], net.surfs['flat'], net.surfs['sa'], net.surfs['vols']) = \
        points, tris, tri_curvs, curvs, funcs, coms, flats, sas, vols
    # Set the dataframe elements
    # Get the curvature in the 95th percentile
    my_surf_curvs = net.surfs['curv'].to_list()
    my_surf_curvs.sort()
    try:
        net.max_curv = my_surf_curvs[min(int(0.99 * len(my_surf_curvs)), len(my_surf_curvs) - 1)]
    except IndexError:
        net.max_curv = 0
    print("\r                                                                                             ", end='')
    net.metrics['surf'] = time.perf_counter() - net.metrics['start'] - net.metrics['vert'] - net.metrics['con']
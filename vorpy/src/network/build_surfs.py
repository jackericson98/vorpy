import time
from vorpy.src.calculations import calc_surf_sa
from vorpy.src.calculations import calc_tetra_vol
from vorpy.src.network.build_surf import build_surf
from vorpy.src.calculations.surface_energy import calc_surface_energy_geometry_from_curvatures


def _add_timing(timing, key, elapsed):
    timing[key] = timing.get(key, 0.0) + elapsed


def _print_surface_timing(total_elapsed, total_surfs, valid_surfs, invalid_surfs,
                          outer_timing, build_timing, total_points, total_tris):
    """Print a compact timing profile for the complete surface stage."""
    print()
    print('=' * 70)
    print('SURFACE BUILD TIMING')
    print('=' * 70)

    rows = [
        ('Data lookup', outer_timing.get('lookup', 0.0)),
        ('Build surface', outer_timing.get('build_surf', 0.0)),
        ('Volume atom 1', outer_timing.get('volume_0', 0.0)),
        ('Volume atom 2', outer_timing.get('volume_1', 0.0)),
        ('Surface area', outer_timing.get('surface_area', 0.0)),
        ('Energy geometry', outer_timing.get('energy_geometry', 0.0)),
        ('Result storage', outer_timing.get('storage', 0.0)),
        ('DataFrame assignment', outer_timing.get('dataframe_assignment', 0.0)),
        ('Curvature percentile', outer_timing.get('curvature_percentile', 0.0)),
    ]

    for name, elapsed in rows:
        pct = 100.0 * elapsed / total_elapsed if total_elapsed > 0 else 0.0
        print(f'{name:<28}{elapsed:>10.4f} s  {pct:>7.2f} %')

    accounted = sum(value for _, value in rows)
    other = max(0.0, total_elapsed - accounted)
    pct = 100.0 * other / total_elapsed if total_elapsed > 0 else 0.0
    print(f'{"Other / loop overhead":<28}{other:>10.4f} s  {pct:>7.2f} %')
    print('-' * 70)
    print(f'{"TOTAL":<28}{total_elapsed:>10.4f} s  {100.0:>7.2f} %')

    print()
    print('BUILD_SURF INTERNAL TIMING')
    print('-' * 70)
    internal_order = [
        ('Surface function', 'surf_func'),
        ('Flatness check', 'flat_check'),
        ('Perimeter', 'perimeter'),
        ('Surface COM', 'get_com'),
        ('Project perimeter', 'project_perimeter'),
        ('Project COM', 'project_com'),
        ('Triangulation', 'triangulate'),
        ('Project hyperboloid', 'project_hyperboloid'),
        ('Mean curvature', 'mean_curvature'),
        ('Gaussian curvature', 'gauss_curvature'),
        ('Flat unprojection', 'unproject_flat'),
        ('Flat curvature init', 'flat_curvature_init'),
    ]

    build_total = build_timing.get('total', 0.0)
    for name, key in internal_order:
        elapsed = build_timing.get(key, 0.0)
        if elapsed == 0.0:
            continue
        pct = 100.0 * elapsed / build_total if build_total > 0 else 0.0
        print(f'{name:<28}{elapsed:>10.4f} s  {pct:>7.2f} %')

    internal_accounted = sum(build_timing.get(key, 0.0) for _, key in internal_order)
    internal_other = max(0.0, build_total - internal_accounted)
    pct = 100.0 * internal_other / build_total if build_total > 0 else 0.0
    print(f'{"Internal overhead":<28}{internal_other:>10.4f} s  {pct:>7.2f} %')
    print('-' * 70)
    print(f'{"BUILD_SURF TOTAL":<28}{build_total:>10.4f} s  {100.0 if build_total else 0.0:>7.2f} %')

    print()
    print('SURFACE BUILD SIZE')
    print(f'Surfaces requested: {total_surfs:,}')
    print(f'Valid surfaces:     {valid_surfs:,}')
    print(f'Invalid surfaces:   {invalid_surfs:,}')
    print(f'Points:             {total_points:,}')
    print(f'Triangles:          {total_tris:,}')
    if valid_surfs:
        print(f'Points / surface:   {total_points / valid_surfs:,.1f}')
        print(f'Tris / surface:     {total_tris / valid_surfs:,.1f}')
    print('=' * 70)


def build_surfs(net, store_points=True):
    """Build all network surfaces while collecting detailed timing metrics."""
    stage_start = time.perf_counter()
    outer_timing = {}
    build_timing = {}

    points, tris, mean_tri_curvs, mean_curvs, avg_mean_curvs = [], [], [], [], []
    gauss_tri_curvs, gauss_curvs, avg_gauss_curvs = [], [], []
    int_mean_curvs, int_mean_curv_sqs, int_gauss_curvs, surf_energies = [], [], [], []
    funcs, coms, flats, sas, vols, surf_locs = [], [], [], [], [], []

    total_surfs = len(net.surfs)
    valid_surfs = 0
    invalid_surfs = 0
    total_points = 0
    total_tris = 0
    last_update = 0.0

    net.update_progress('Building surfaces | Initializing', 0.0)

    for count, (i, surf) in enumerate(net.surfs.iterrows(), start=1):
        lookup_start = time.perf_counter()
        rads = [net.balls['rad'][_] for _ in surf['balls']]
        locs = [net.balls['loc'][_] for _ in surf['balls']]
        nums = [net.balls['num'][_] for _ in surf['balls']]
        epnts = [net.edges['points'][_] for _ in surf['edges']]

        if rads[0] > rads[1]:
            rads, locs, nums = [rads[1], rads[0]], [locs[1], locs[0]], [nums[1], nums[0]]
        _add_timing(outer_timing, 'lookup', time.perf_counter() - lookup_start)

        build_start = time.perf_counter()
        my_surf = build_surf(
            locs=locs,
            rads=rads,
            epnts=epnts,
            res=net.settings['surf_res'],
            net_type=net.settings['net_type'],
            timing=build_timing,
        )
        _add_timing(outer_timing, 'build_surf', time.perf_counter() - build_start)

        if my_surf is None:
            invalid_surfs += 1
            net.surfs.drop(index=i, inplace=True)
            continue

        valid_surfs += 1
        (
            surf_points, surf_tris,
            mean_surf_tri_curvs, mean_surf_curv, avg_mean_surf_curv,
            gauss_surf_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv,
            surf_func, surf_com, surf_flat, surf_loc
        ) = my_surf

        total_points += len(surf_points)
        total_tris += len(surf_tris)

        timer = time.perf_counter()
        sv0 = sum(
            calc_tetra_vol(locs[0], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]])
            for tri in surf_tris
        )
        _add_timing(outer_timing, 'volume_0', time.perf_counter() - timer)

        timer = time.perf_counter()
        sv1 = sum(
            calc_tetra_vol(locs[1], surf_points[tri[0]], surf_points[tri[1]], surf_points[tri[2]])
            for tri in surf_tris
        )
        _add_timing(outer_timing, 'volume_1', time.perf_counter() - timer)

        timer = time.perf_counter()
        sa = calc_surf_sa(tris=surf_tris, points=surf_points)
        _add_timing(outer_timing, 'surface_area', time.perf_counter() - timer)

        timer = time.perf_counter()
        energy_geometry = calc_surface_energy_geometry_from_curvatures(
            points=surf_points,
            tris=surf_tris,
            mean_curvatures=mean_surf_tri_curvs,
            gaussian_curvatures=gauss_surf_tri_curvs,
            area=sa,
        )
        _add_timing(outer_timing, 'energy_geometry', time.perf_counter() - timer)

        current_time = time.perf_counter()
        if current_time - last_update >= 0.25 or count == total_surfs:
            percentage = 100.0 * count / max(total_surfs, 1)
            net.update_progress(
                f'Building surfaces: {count:,} / {total_surfs:,}',
                percentage,
            )
            last_update = current_time

        timer = time.perf_counter()
        if store_points:
            points.append(surf_points)
            tris.append(surf_tris)
            mean_tri_curvs.append(mean_surf_tri_curvs)
            gauss_tri_curvs.append(gauss_surf_tri_curvs)
        else:
            points.append([])
            tris.append([])
            mean_tri_curvs.append([])
            gauss_tri_curvs.append([])

        mean_curvs.append(mean_surf_curv)
        avg_mean_curvs.append(avg_mean_surf_curv)
        gauss_curvs.append(gauss_surf_curv)
        avg_gauss_curvs.append(avg_gauss_surf_curv)

        int_mean_curvs.append(energy_geometry['Integrated Mean Curvature'])
        int_mean_curv_sqs.append(energy_geometry['Integrated Mean Curvature Squared'])
        int_gauss_curvs.append(energy_geometry['Integrated Gaussian Curvature'])
        surf_energies.append(2.0 * energy_geometry['Integrated Mean Curvature Squared'])

        funcs.append(surf_func)
        coms.append(surf_com)
        flats.append(surf_flat)
        sas.append(sa)
        vols.append({nums[0]: sv0, nums[1]: sv1})
        surf_locs.append(surf_loc)
        _add_timing(outer_timing, 'storage', time.perf_counter() - timer)

    timer = time.perf_counter()
    (
        net.surfs['points'], net.surfs['tris'],
        net.surfs['mean_tri_curvs'], net.surfs['mean_curv'], net.surfs['avg_mean_curv'],
        net.surfs['gauss_tri_curvs'], net.surfs['gauss_curv'], net.surfs['avg_gauss_curv'],
        net.surfs['int_mean_curv'], net.surfs['int_mean_curv_sq'], net.surfs['int_gauss_curv'],
        net.surfs['surf_energy'], net.surfs['func'], net.surfs['com'], net.surfs['flat'],
        net.surfs['sa'], net.surfs['vols'], net.surfs['loc']
    ) = (
        points, tris,
        mean_tri_curvs, mean_curvs, avg_mean_curvs,
        gauss_tri_curvs, gauss_curvs, avg_gauss_curvs,
        int_mean_curvs, int_mean_curv_sqs, int_gauss_curvs,
        surf_energies, funcs, coms, flats, sas, vols, surf_locs
    )
    _add_timing(outer_timing, 'dataframe_assignment', time.perf_counter() - timer)

    timer = time.perf_counter()
    my_surf_curvs = net.surfs['mean_curv'].to_list()
    if net.settings['surf_scheme'] == 'gauss':
        my_surf_curvs = net.surfs['gauss_curv'].to_list()

    my_surf_curvs.sort()
    try:
        net.max_curv = my_surf_curvs[min(int(0.99 * len(my_surf_curvs)), len(my_surf_curvs) - 1)]
    except IndexError:
        net.max_curv = 0
    _add_timing(outer_timing, 'curvature_percentile', time.perf_counter() - timer)

    net.update_progress('Building surfaces', 100.0)
    net.metrics['surf'] = time.perf_counter() - net.metrics['start'] - net.metrics['vert'] - net.metrics['con']
    net.update_progress(
        f'Building surfaces: {len(net.surfs):,} / {len(net.surfs):,}',
        100.0,
    )

    total_elapsed = time.perf_counter() - stage_start
    _print_surface_timing(
        total_elapsed=total_elapsed,
        total_surfs=total_surfs,
        valid_surfs=valid_surfs,
        invalid_surfs=invalid_surfs,
        outer_timing=outer_timing,
        build_timing=build_timing,
        total_points=total_points,
        total_tris=total_tris,
    )
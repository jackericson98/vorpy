import time
import pandas as pd
from vorpy.src.calculations import box_search
from vorpy.src.calculations import get_balls
from vorpy.src.calculations import calc_dist
from vorpy.src.network.find_verts import find_verts
from vorpy.src.output import write_verts


def vertex_is_allowed(net, vertex_balls):
    """
    Determine whether a vertex belongs to the requested network scope.

    Normal networks retain a vertex when at least one defining ball belongs
    to ``net.group``. Interface networks retain a vertex only when its four
    defining balls include at least one ball from each interface side.

    Parameters
    ----------
    net : Network
        Network providing ``group`` and optional ``iface_grps`` constraints.
    vertex_balls : collection of int
        Defining ball indices for the candidate vertex.

    Returns
    -------
    bool
        True when the vertex belongs to the requested network scope.
    """
    vertex_balls = set(vertex_balls)

    if net.iface_grps is not None:
        return all(
            bool(vertex_balls.intersection(group_indices))
            for group_indices in net.iface_grps
        )

    if net.group is None:
        return True

    return bool(vertex_balls.intersection(net.group))


def _load_cached_vertex_state(net):
    """
    Reconstruct native vertex-search arrays from cached interface state.

    Cached records are filtered against the current network scope and converted
    back into the parallel arrays expected by ``find_verts``. Ball-to-vertex
    adjacency is rebuilt so cached vertices participate normally in duplicate
    detection and traversal.

    Returns
    -------
    tuple or None
        ``(vert_ndxs, vlocs, vrads, vloc2s, vrad2s, b_verts)`` when usable cached
        state exists, otherwise ``None``.
    """
    records = list(getattr(net, 'cached_interface_vertex_state', None) or [])
    if not records:
        return None

    records.sort(key=lambda record: record.get('source_index', 0))
    vert_ndxs = []
    vlocs = []
    vrads = []
    vloc2s = []
    vrad2s = []
    b_verts = [[] for _ in range(len(net.balls))]

    for record in records:
        row = record['data']
        balls = [int(ball) for ball in row.get('balls', [])]
        if not vertex_is_allowed(net, balls):
            continue

        vertex_index = len(vert_ndxs)
        vert_ndxs.append(balls)
        vlocs.append([float(value) for value in row.get('loc', [])])
        vrads.append(float(row.get('rad', 0.0)))

        loc2 = row.get('loc2', None)
        if loc2 is None or all(value is None for value in loc2):
            vloc2s.append([None, None, None])
            vrad2s.append(None)
        else:
            vloc2s.append([float(value) for value in loc2])
            rad2 = row.get('rad2', None)
            vrad2s.append(None if rad2 is None else float(rad2))

        for ball in balls:
            b_verts[ball].append(vertex_index)

    if not vert_ndxs:
        return None

    return vert_ndxs, vlocs, vrads, vloc2s, vrad2s, b_verts


def _store_native_vertex_state(net, vert_ndxs, vlocs, vrads, vloc2s, vrad2s):
    """
    Store the native vertex-search representation for interface reuse.

    The unexpanded primary/secondary vertex representation is preserved so future
    interface calculations can resume discovery without reconstructing doublets
    from exported dataframe rows.
    """
    net.interface_vertex_state = []
    for balls, loc, rad, loc2, rad2 in zip(
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s):
        net.interface_vertex_state.append({
            'balls': [int(ball) for ball in balls],
            'loc': [float(value) for value in loc],
            'rad': float(rad),
            'loc2': (
                [None, None, None]
                if loc2 is None or all(value is None for value in loc2)
                else [float(value) for value in loc2]
            ),
            'rad2': None if rad2 is None else float(rad2),
        })


def _get_interface_reseed_candidates(net, sphere_check_list):
    """
    Return unresolved balls that can geometrically participate in an interface seed.

    AW ``find_v0`` currently limits seed discovery to ``max_vert / 10``.
    Therefore two balls from opposite interface sides that define the same
    seed vertex cannot have a surface-to-surface separation greater than
    twice that seed radius.

    Candidates are taken from one interface side only because every valid
    interface component necessarily contains at least one ball from each side.
    """
    if net.iface_grps is None:
        return sphere_check_list

    side1, side2 = net.iface_grps
    unresolved = set(sphere_check_list)

    locs = net.balls['loc'].to_numpy()
    rads = net.balls['rad'].to_numpy()

    seed_max = net.settings['max_vert'] / 10
    max_surface_sep = 2 * seed_max

    side2 = set(side2)
    side2_max_rad = max(rads[ball] for ball in side2)

    candidates = []

    for ball in side1:
        if ball not in unresolved:
            continue

        # Any opposite-side ball satisfying the exact surface-separation
        # bound must lie inside this conservative center-distance search.
        search_dist = max_surface_sep + rads[ball] + side2_max_rad
        ball_box = box_search(locs[ball])
        nearby = get_balls([ball_box], dist=search_dist)

        min_surface_sep = None

        for other_ball in nearby:
            if other_ball not in side2:
                continue

            surface_sep = (
                calc_dist(locs[ball], locs[other_ball])
                - rads[ball]
                - rads[other_ball]
            )

            if surface_sep <= max_surface_sep:
                if min_surface_sep is None or surface_sep < min_surface_sep:
                    min_surface_sep = surface_sep

        if min_surface_sep is not None:
            candidates.append((min_surface_sep, ball))

    # Attack the most plausible interface seeds first.
    # Store farthest-to-nearest so pop() retrieves the closest candidate first.
    candidates.sort(key=lambda item: item[0], reverse=True)

    return [ball for _, ball in candidates]


def _print_find_verts_timing(net, timer, total):
    """Print a compact, meaningful vertex-search timing breakdown."""
    if not net.settings.get('verbose', False):
        return

    grouped = {
        'Setup + seed discovery': (
            timer.get('setup', 0.0)
            + timer.get('cache_load', 0.0)
            + timer.get('fv_setup', 0.0)
            + timer.get('seed_discovery', 0.0)
            + timer.get('seed_insert', 0.0)
        ),
        'Site-container search': timer.get('site_container', 0.0),
        'Vertex bookkeeping': (
            timer.get('edge_prep', 0.0)
            + timer.get('candidate_validation', 0.0)
            + timer.get('state_append', 0.0)
            + timer.get('adjacency_updates', 0.0)
        ),
        'Progress reporting': timer.get('progress', 0.0),
        'Post-search processing': (
            timer.get('encapsulation', 0.0)
            + timer.get('interface_reseed_setup', 0.0)
            + timer.get('cache_store', 0.0)
            + timer.get('doublets', 0.0)
            + timer.get('dataframe', 0.0)
        ),
        'Vertex export': timer.get('export', 0.0),
    }

    print("\n" + "=" * 70)
    print("FIND VERTICES TIMING")
    print("=" * 70)

    for label, elapsed in grouped.items():
        pct = 100.0 * elapsed / total if total > 0 else 0.0
        print(f"{label:<28} {elapsed:10.4f} s  {pct:6.2f} %")

    measured = sum(grouped.values())
    other = max(total - measured, 0.0)
    pct = 100.0 * other / total if total > 0 else 0.0
    print(f"{'Other / unmeasured':<28} {other:10.4f} s  {pct:6.2f} %")
    print("-" * 70)
    print(f"{'TOTAL':<28} {total:10.4f} s  100.00 %")

    print("\nFIND VERTICES SIZE")
    print(f"Balls requested:    {len(net.group or []):,}")
    print(f"Vertices found:     {0 if net.verts is None else len(net.verts):,}")
    if net.verts is not None and 'dub' in net.verts:
        print(f"Doublet rows:       {sum(1 for value in net.verts['dub'] if value in (1, 2)):,}")
    print(f"Site searches:      {int(timer.get('edge_search_calls', 0)):,}")
    print(f"Accepted vertices:  {int(timer.get('accepted_vertices', 0)):,}")
    print(f"No-site results:    {int(timer.get('rejected_none', 0)):,}")
    print("=" * 70)


def find_net_verts(net):
    """
    Find and store all vertices belonging to a network.

    Detailed timings are always collected in ``net.vert_timing``.
    They are printed only when ``net.settings['verbose']`` is True.
    """
    vert_start = time.perf_counter()
    timer = {
        'setup': 0.0,
        'cache_load': 0.0,
        'encapsulation': 0.0,
        'interface_reseed_setup': 0.0,
        'cache_store': 0.0,
        'doublets': 0.0,
        'dataframe': 0.0,
        'export': 0.0,
    }

    # --------------------------------------------------------------
    # Setup
    # --------------------------------------------------------------
    t = time.perf_counter()
    if net.group is None:
        net.group = net.balls['num'].tolist()

    sphere_check_list = net.group.copy()
    net.update_progress("Finding vertices | Initializing", 0.0)
    timer['setup'] += time.perf_counter() - t

    # --------------------------------------------------------------
    # Cached vertex state
    # --------------------------------------------------------------
    t = time.perf_counter()
    cached_state = _load_cached_vertex_state(net)
    if cached_state is None:
        vert_ndxs = vlocs = vrads = vloc2s = vrad2s = averts = None
    else:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, averts = cached_state
    timer['cache_load'] += time.perf_counter() - t

    # --------------------------------------------------------------
    # Initial vertex search
    # --------------------------------------------------------------
    my_guuy = find_verts(
        net=net,
        locs=net.balls['loc'].to_numpy(),
        rads=net.balls['rad'].to_numpy(),
        max_vert=net.settings['max_vert'],
        net_type=net.settings['net_type'],
        check_ndxs=sphere_check_list,
        my_group=net.group,
        iface_grps=net.iface_grps,
        vert_ndxs=vert_ndxs,
        vlocs=vlocs,
        vrads=vrads,
        vloc2s=vloc2s,
        vrad2s=vrad2s,
        b_verts=averts,
        start_time=net.metrics['start'],
        vert_box=net.settings['foam_box'],
        box=net.box['verts'],
        timing=timer
    )

    if my_guuy is not None:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, sphere_check_list, averts = my_guuy
    elif cached_state is None:
        total = time.perf_counter() - vert_start
        net.vert_timing = timer.copy()
        net.vert_timing['total'] = total
        net.metrics['vert'] = time.perf_counter() - net.metrics['start']
        _print_find_verts_timing(net, timer, total)
        return

    # --------------------------------------------------------------
    # Encapsulation checks
    # --------------------------------------------------------------
    t = time.perf_counter()
    if len(sphere_check_list) > 0 and net.iface_grps is None:
        skip_nums = []
        max_ball_rad = max(net.balls['rad'])

        for sphere in sphere_check_list:
            sphere_rad = net.balls['rad'][sphere]
            sphere_loc = net.balls['loc'][sphere]
            sphere_box = box_search(sphere_loc)
            close_spheres = get_balls(
                sphere_box,
                dist=max_ball_rad - sphere_rad
            )

            if close_spheres is not None:
                for sphere2 in close_spheres:
                    if calc_dist(
                        sphere_loc,
                        net.balls['loc'][sphere2]
                    ) < abs(net.balls['rad'][sphere2] - sphere_rad):
                        print(
                            "\nUh oh! Ball # {} is fully encapsulated by ball # {}! Skipping {}"
                            .format(sphere, sphere2, sphere)
                        )
                        skip_nums.append(sphere)
                        break

        for sphere in skip_nums:
            sphere_check_list.pop(sphere_check_list.index(sphere))

    timer['encapsulation'] += time.perf_counter() - t

    # --------------------------------------------------------------
    # Interface reseed preparation
    # --------------------------------------------------------------
    t = time.perf_counter()
    if net.iface_grps is not None:
        sphere_check_list[:] = _get_interface_reseed_candidates(
            net,
            sphere_check_list
        )
    timer['interface_reseed_setup'] += time.perf_counter() - t

    # --------------------------------------------------------------
    # Additional seed searches
    # --------------------------------------------------------------
    reseed_count = 0

    while sphere_check_list:
        seed_ball = sphere_check_list.pop()
        reseed_count += 1

        my_guuy = find_verts(
            b0=seed_ball,
            locs=net.balls['loc'].to_numpy(),
            rads=net.balls['rad'].to_numpy(),
            max_vert=net.settings['max_vert'],
            net_type=net.settings['net_type'],
            check_ndxs=sphere_check_list,
            my_group=net.group,
            iface_grps=net.iface_grps,
            vert_ndxs=vert_ndxs,
            vlocs=vlocs,
            vrads=vrads,
            vloc2s=vloc2s,
            vrad2s=vrad2s,
            start_time=net.metrics['start'],
            vert_box=net.settings['foam_box'],
            b_verts=averts,
            box=net.box['verts'],
            seed_timeout=0.05,
            net=net,
            timing=timer
        )

        if my_guuy is not None:
            (
                vert_ndxs,
                vlocs,
                vrads,
                vloc2s,
                vrad2s,
                sphere_check_list,
                averts
            ) = my_guuy

        if (
            'ball_type' in net.settings
            and net.settings['ball_type'] == 'foam'
            and len(sphere_check_list) <= 0.25 * len(net.balls['loc'])
        ):
            print(f'Missing Ball Indices:\n{sphere_check_list}\n')
            break

    # --------------------------------------------------------------
    # Store native interface state
    # --------------------------------------------------------------
    t = time.perf_counter()
    if net.iface_grps is not None:
        _store_native_vertex_state(
            net,
            vert_ndxs,
            vlocs,
            vrads,
            vloc2s,
            vrad2s
        )
    timer['cache_store'] += time.perf_counter() - t

    # --------------------------------------------------------------
    # Doublet expansion
    # --------------------------------------------------------------
    t = time.perf_counter()
    doublets = [0] * len(vert_ndxs)
    i = 0

    while i < len(vlocs):
        if vrad2s[i] is not None:
            vert_ndxs.insert(i + 1, vert_ndxs[i])
            vlocs.insert(i + 1, vloc2s[i])
            vrads.insert(i + 1, vrad2s[i])
            doublets[i] = 2
            doublets.insert(i + 1, 1)
            vrad2s.insert(i + 1, None)
            vloc2s.insert(i + 1, [None, None, None])
        i += 1

    timer['doublets'] += time.perf_counter() - t

    # --------------------------------------------------------------
    # Vertex DataFrame
    # --------------------------------------------------------------
    t = time.perf_counter()
    net.verts = pd.DataFrame({
        "balls": vert_ndxs,
        'loc': vlocs,
        'rad': vrads,
        'dub': doublets
    })
    timer['dataframe'] += time.perf_counter() - t

    net.update_progress(
        f"Finding vertices: {len(net.verts):,} / {len(net.verts):,}",
        100.0
    )

    # Historical metric remains system-relative for compatibility.
    net.metrics['vert'] = time.perf_counter() - net.metrics['start']

    # --------------------------------------------------------------
    # Vertex output
    # --------------------------------------------------------------
    t = time.perf_counter()
    write_verts(net)
    timer['export'] += time.perf_counter() - t

    total = time.perf_counter() - vert_start

    net.vert_timing = timer.copy()
    net.vert_timing['total'] = total
    net.vert_timing['reseed_count'] = reseed_count

    _print_find_verts_timing(net, timer, total)
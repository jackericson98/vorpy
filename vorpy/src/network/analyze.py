import time
import numpy as np
from numpy import pi, sqrt
from time import perf_counter as now

from vorpy.src.calculations import calc_sphericity
from vorpy.src.calculations import calc_isoperimetric_quotient
from vorpy.src.calculations import calc_contacts_cached
from vorpy.src.calculations import calc_cell_point_properties_cached
from vorpy.src.calculations import calc_cell_mass_properties_cached


def _rmsd(values, average):
    """Return the RMS deviation of values from average."""
    if not values:
        return 0.0
    return sqrt(sum((value - average) ** 2 for value in values) / len(values))


def _print_analysis_timing(timer, total):
    """Print a detailed analysis timing breakdown."""
    print("\n" + "=" * 70)
    print("ANALYSIS TIMING")
    print("=" * 70)

    labels = [
        ('setup', 'Setup / cache construction'),
        ('surface_gather', 'Surface gathering'),
        ('completeness', 'Cell completeness'),
        ('basic', 'Surface area + volume'),
        ('curvs', 'Curvature'),
        ('geometric', 'Geometric metrics'),
        ('neighbors_1', 'First neighbors'),
        ('neighbors_2', 'Second neighbors'),
        ('spikes', 'Spikes + bounding box'),
        ('contacts', 'Contacts'),
        ('com', 'COM + moment of inertia'),
        ('moi', 'Moment of inertia only'),
        ('b_box', 'Bounding box only'),
        ('surface_assign', 'Surface assignment'),
        ('ball_assign', 'Ball assignment'),
    ]

    for key, label in labels:
        elapsed = timer.get(key, 0.0)
        pct = 100.0 * elapsed / total if total > 0 else 0.0
        print(f"{label:<28} {elapsed:10.4f} s  {pct:6.2f} %")

    measured = sum(timer.get(key, 0.0) for key, _ in labels)
    other = max(total - measured, 0.0)
    pct = 100.0 * other / total if total > 0 else 0.0
    print(f"{'Other / loop overhead':<28} {other:10.4f} s  {pct:6.2f} %")
    print("-" * 70)
    print(f"{'TOTAL':<28} {total:10.4f} s  100.00 %")


def analyze(
    net,
    complicated=True,
    spikes=None,
    contacts=None,
    second_neighbors=None,
    com=None,
    moi=None,
    bounding_box=None,
):
    """
    Analyze the cells in ``net.group``.

    This implementation keeps the historical output columns and default
    behavior while removing most pandas work from the hot per-cell loop.

    Parameters
    ----------
    net
        VorPy network object.
    complicated : bool, default=True
        Backward-compatible master switch for the expensive metrics.
    spikes, contacts, second_neighbors, com, moi, bounding_box : bool or None
        Optional per-feature switches. ``None`` inherits ``complicated``.
    """
    analysis_start = now()
    net.update_progress("Analyzing network | Initializing", 0.0)

    # Preserve the old complicated=True/False behavior unless a feature is
    # explicitly overridden by the caller.
    spikes = complicated if spikes is None else spikes
    contacts = complicated if contacts is None else contacts
    second_neighbors = complicated if second_neighbors is None else second_neighbors
    com = complicated if com is None else com
    moi = complicated if moi is None else moi
    bounding_box = complicated if bounding_box is None else bounding_box

    timer = {
        'setup': 0.0,
        'surface_gather': 0.0,
        'completeness': 0.0,
        'basic': 0.0,
        'curvs': 0.0,
        'geometric': 0.0,
        'neighbors_1': 0.0,
        'neighbors_2': 0.0,
        'spikes': 0.0,
        'contacts': 0.0,
        'com': 0.0,
        'moi': 0.0,
        'b_box': 0.0,
        'surface_assign': 0.0,
        'ball_assign': 0.0,
    }

    setup_start = now()

    # ------------------------------------------------------------------
    # Cache DataFrame columns once. Object-valued columns remain object
    # arrays, but this still avoids Series creation, .iloc row access, and
    # repeated DataFrame slicing throughout the hot loop.
    # ------------------------------------------------------------------
    n_balls = len(net.balls)
    n_surfs = len(net.surfs)

    ball_nums = net.balls['num'].to_numpy()
    ball_surfs_all = net.balls['surfs'].to_numpy()
    ball_verts_all = net.balls['verts'].to_numpy()
    ball_edges_all = net.balls['edges'].to_numpy()
    ball_locs = net.balls['loc'].to_numpy()
    ball_locs_matrix = np.asarray([np.asarray(loc, dtype=float) for loc in ball_locs])
    ball_rads = net.balls['rad'].to_numpy(dtype=float)
    ball_index = net.balls.index.to_numpy()

    surf_sa = net.surfs['sa'].to_numpy()
    surf_vols = net.surfs['vols'].to_numpy()
    surf_balls = net.surfs['balls'].to_numpy()
    surf_mean_curv = net.surfs['mean_curv'].to_numpy()
    surf_gauss_curv = net.surfs['gauss_curv'].to_numpy()
    surf_int_mean = net.surfs['int_mean_curv'].to_numpy()
    surf_int_mean_sq = net.surfs['int_mean_curv_sq'].to_numpy()
    surf_int_gauss = net.surfs['int_gauss_curv'].to_numpy()

    # Normalize surface geometry once. Every analyzed cell can now reuse these
    # arrays directly without DataFrame slicing, dictionary construction, or
    # repeated np.asarray conversion inside geometry helpers.
    raw_surf_points = net.surfs['points'].to_numpy()
    raw_surf_tris = net.surfs['tris'].to_numpy()
    surf_points = np.empty(n_surfs, dtype=object)
    surf_tris = np.empty(n_surfs, dtype=object)
    for surf_id in range(n_surfs):
        points = np.asarray(raw_surf_points[surf_id], dtype=np.float64)
        if points.size == 0:
            points = np.empty((0, 3), dtype=np.float64)
        elif points.ndim != 2:
            points = points.reshape((-1, 3))
        surf_points[surf_id] = np.ascontiguousarray(points)

        tris = np.asarray(raw_surf_tris[surf_id], dtype=np.int64)
        if tris.size == 0:
            tris = np.empty((0, 3), dtype=np.int64)
        elif tris.ndim != 2:
            tris = tris.reshape((-1, 3))
        surf_tris[surf_id] = np.ascontiguousarray(tris)

    edge_balls = net.edges['balls'].to_numpy()
    vert_edges = net.verts['edges'].to_numpy()
    vert_balls = net.verts['balls'].to_numpy()

    # Ball numbers have historically matched DataFrame positions, but keeping
    # an explicit lookup costs little and makes this routine safer.
    num_to_pos = {int(num): pos for pos, num in enumerate(ball_nums)}

    # Group-only iteration: irrelevant balls stay at their preallocated zero
    # defaults and are never visited by the expensive loop.
    group_nums = [int(num) for num in net.group if int(num) in num_to_pos]
    group_set = set(group_nums)
    n_group = len(group_nums)

    # Precompute immutable edge-ball sets for the completeness fallback.
    edge_ball_sets = [frozenset(balls) for balls in edge_balls]

    # ------------------------------------------------------------------
    # Build the atom adjacency graph once from solved surfaces.
    # Each entry contains (neighbor_ball_num, surface_id).
    # ------------------------------------------------------------------
    ball_neighbors = [[] for _ in range(n_balls)]
    neighbor_sets = [set() for _ in range(n_balls)]

    for surf_id, balls in enumerate(surf_balls):
        if balls is None or len(balls) < 2:
            continue

        # Normal Voronoi surfaces are pairwise. The nested loop also behaves
        # sensibly if a future representation contains >2 balls.
        for i, ball_num in enumerate(balls):
            ball_num = int(ball_num)
            ball_pos = num_to_pos.get(ball_num)
            if ball_pos is None:
                continue

            for j, neighbor_num in enumerate(balls):
                if i == j:
                    continue
                neighbor_num = int(neighbor_num)
                if neighbor_num not in num_to_pos:
                    continue
                ball_neighbors[ball_pos].append((neighbor_num, surf_id))
                neighbor_sets[ball_pos].add(neighbor_num)

    # ------------------------------------------------------------------
    # Preallocate every output for the full ball table. This removes the large
    # append_0 blocks and allows direct assignment by DataFrame position.
    # ------------------------------------------------------------------
    b_vols = [0.0] * n_balls
    b_sas = [0.0] * n_balls
    b_cell = [0] * n_balls

    b_max_mean_curvs = [0.0] * n_balls
    b_avg_mean_surf_curvs = [0.0] * n_balls
    b_max_gauss_curvs = [0.0] * n_balls
    b_avg_gauss_surf_curvs = [0.0] * n_balls

    b_int_mean_curvs = [0.0] * n_balls
    b_int_mean_curv_sqs = [0.0] * n_balls
    b_int_gauss_curvs = [0.0] * n_balls

    b_sphrctys = [0.0] * n_balls
    b_isopmqs = [0.0] * n_balls

    num_nbors = [0] * n_balls
    near_nbors = [0] * n_balls
    near_nbor_dists = [0.0] * n_balls
    nbor_lyr_rmsds = [0] * n_balls
    nbor_dst_avgs = [0] * n_balls
    b_inner = [0] * n_balls

    b_min_spikes = [0.0] * n_balls
    b_max_spikes = [0.0] * n_balls

    contact_areas = [0.0] * n_balls
    non_olap_vols = [0.0] * n_balls
    olap_vols = [0.0] * n_balls
    num_olaps = [0] * n_balls

    coms = [0.0] * n_balls
    mois = [0.0] * n_balls
    b_boxs = [0.0] * n_balls

    # Array-based surface tracking replaces nested dictionaries and repeated
    # .loc writes at the end of analysis.
    surface_contact_area = np.zeros(n_surfs, dtype=float)
    surface_overlap = np.zeros(n_surfs, dtype=float)

    timer['setup'] += now() - setup_start

    # Progress counter counts only balls that are actually analyzed.
    count = 0
    last_update = now()

    # ==================================================================
    # HOT LOOP: ONLY BALLS IN net.group
    # ==================================================================
    for ball_num in group_nums:
        ball_pos = num_to_pos[ball_num]
        surf_ids = ball_surfs_all[ball_pos]

        # Keep historical zero defaults for cells with no surfaces.
        if surf_ids is None or len(surf_ids) == 0:
            count += 1
            continue

        # --------------------------------------------------------------
        # Gather all basic surface/cell quantities in one surface pass.
        # --------------------------------------------------------------
        t = now()

        sa = 0.0
        volume = 0.0
        max_mean_curv = -float('inf')
        max_gauss_curv = -float('inf')
        int_mean_curv = 0.0
        int_mean_curv_sq = 0.0
        int_gauss_curv = 0.0

        for surf_id in surf_ids:
            surf_id = int(surf_id)
            sa_i = surf_sa[surf_id]
            sa += sa_i

            vols_i = surf_vols[surf_id]
            try:
                volume += vols_i[ball_num]
            except (IndexError, KeyError, TypeError):
                # Fallback for non-positional volume containers.
                volume += vols_i[ball_pos]

            mean_i = surf_mean_curv[surf_id]
            gauss_i = surf_gauss_curv[surf_id]
            if mean_i > max_mean_curv:
                max_mean_curv = mean_i
            if gauss_i > max_gauss_curv:
                max_gauss_curv = gauss_i

            int_mean_curv += surf_int_mean[surf_id]
            int_mean_curv_sq += surf_int_mean_sq[surf_id]
            int_gauss_curv += surf_int_gauss[surf_id]

        timer['surface_gather'] += now() - t

        if sa == 0:
            count += 1
            continue

        # --------------------------------------------------------------
        # Completeness check using cached arrays/frozensets.
        # --------------------------------------------------------------
        t = now()
        ball_verts = ball_verts_all[ball_pos]
        ball_edges = ball_edges_all[ball_pos]
        complete = True

        ball_index_value = ball_index[ball_pos]
        for vert in ball_verts:
            vert = int(vert)
            owning_edge_count = 0
            for edge in vert_edges[vert]:
                # Preserve the historical check, which used the DataFrame
                # row index (k) rather than ball['num'] here.
                if ball_index_value in edge_ball_sets[int(edge)]:
                    owning_edge_count += 1

            if owning_edge_count != 3:
                vert_ball_set = set(vert_balls[vert])
                new_count = 0
                for edge in ball_edges:
                    if edge_ball_sets[int(edge)].issubset(vert_ball_set):
                        new_count += 1
                if new_count < 3:
                    complete = False
                    break

        if len(ball_verts) < 3 or len(ball_edges) < 4 or len(surf_ids) < 3:
            complete = False

        b_cell[ball_pos] = complete
        timer['completeness'] += now() - t

        # --------------------------------------------------------------
        # Basic cell values were gathered above; assignment is separated
        # in timing so the profiler does not hide completeness cost.
        # --------------------------------------------------------------
        t = now()
        b_sas[ball_pos] = sa
        b_vols[ball_pos] = volume
        timer['basic'] += now() - t

        # --------------------------------------------------------------
        # Curvature
        # --------------------------------------------------------------
        t = now()
        b_max_mean_curvs[ball_pos] = max_mean_curv
        b_max_gauss_curvs[ball_pos] = max_gauss_curv
        b_int_mean_curvs[ball_pos] = int_mean_curv
        b_int_mean_curv_sqs[ball_pos] = int_mean_curv_sq
        b_int_gauss_curvs[ball_pos] = int_gauss_curv
        b_avg_mean_surf_curvs[ball_pos] = int_mean_curv / sa
        b_avg_gauss_surf_curvs[ball_pos] = int_gauss_curv / sa
        timer['curvs'] += now() - t

        # --------------------------------------------------------------
        # Shape metrics
        # --------------------------------------------------------------
        t = now()
        b_sphrctys[ball_pos] = calc_sphericity(volume=volume, surface_area=sa)
        b_isopmqs[ball_pos] = calc_isoperimetric_quotient(volume=volume, surface_area=sa)
        timer['geometric'] += now() - t

        # --------------------------------------------------------------
        # First neighbor layer from the prebuilt adjacency graph.
        # --------------------------------------------------------------
        t = now()
        ball_loc = ball_locs[ball_pos]
        ball_rad = ball_rads[ball_pos]
        adjacency = ball_neighbors[ball_pos]

        neighbors_nums = []
        neighbor_dists = []

        for neighbor_num, surf_id in adjacency:
            neighbor_pos = num_to_pos[neighbor_num]
            delta = ball_locs_matrix[neighbor_pos] - ball_locs_matrix[ball_pos]
            center_dist = sqrt(float(np.dot(delta, delta)))
            neighbor_dist = center_dist - ball_rad - ball_rads[neighbor_pos]

            neighbors_nums.append(neighbor_num)
            neighbor_dists.append(neighbor_dist)

            overlap_dist = max(-neighbor_dist, 0.0)
            if overlap_dist > surface_overlap[surf_id]:
                surface_overlap[surf_id] = overlap_dist

        # A solved cell should have neighbors if it has surfaces; keep safe
        # defaults if malformed topology reaches this point.
        if neighbor_dists:
            b_inner[ball_pos] = group_set.issuperset(neighbors_nums)
            num_nbors[ball_pos] = len(neighbors_nums)

            min_index = min(range(len(neighbor_dists)), key=neighbor_dists.__getitem__)
            min_dist = neighbor_dists[min_index]
            near_nbor_dists[ball_pos] = min_dist
            near_nbors[ball_pos] = neighbors_nums[min_index]

            nbor_dist_avg = sum(neighbor_dists) / len(neighbor_dists)
            nbor_dst_avgs[ball_pos] = [nbor_dist_avg]
            nbor_lyr_rmsds[ball_pos] = [_rmsd(neighbor_dists, nbor_dist_avg)]
        else:
            b_inner[ball_pos] = True
            nbor_dst_avgs[ball_pos] = [0.0]
            nbor_lyr_rmsds[ball_pos] = [0.0]

        timer['neighbors_1'] += now() - t

        # --------------------------------------------------------------
        # Point properties: spikes + bounding box in one compiled pass.
        # --------------------------------------------------------------
        if spikes or bounding_box:
            t = now()
            min_spike, max_spike, box = calc_cell_point_properties_cached(ball_loc, surf_ids, surf_points)
            elapsed = now() - t

            if spikes:
                b_min_spikes[ball_pos] = min_spike
                b_max_spikes[ball_pos] = max_spike

            if bounding_box:
                b_boxs[ball_pos] = box

            # When both are requested (the normal complicated=True path), the
            # single traversal is charged to this combined timing category.
            timer['spikes'] += elapsed

        # --------------------------------------------------------------
        # Contacts / overlap volume
        # --------------------------------------------------------------
        if contacts:
            t = now()
            contact_area, vdw_vol = calc_contacts_cached(ball_loc, ball_rad, surf_ids, surf_points, surf_tris)

            num_olaps[ball_pos] = sum(1 for dist in neighbor_dists if dist < 0)
            contact_areas[ball_pos] = sum(contact_area.values())

            for surf_id in surf_ids:
                surf_id = int(surf_id)
                surface_contact_area[surf_id] = contact_area[surf_id]

            non_olap_vols[ball_pos] = vdw_vol
            olap_vols[ball_pos] = (4.0 / 3.0) * pi * ball_rad ** 3 - vdw_vol
            timer['contacts'] += now() - t

        # --------------------------------------------------------------
        # Second neighbor layer from cached adjacency sets.
        # This replaces get_next_layer() and all pandas/dict conversion in
        # the previous nested neighbor search.
        # --------------------------------------------------------------
        if second_neighbors:
            t = now()
            first_set = neighbor_sets[ball_pos]
            second_set = set()

            # Match the historical behavior: expand the first layer plus the
            # central ball, then remove every first-layer/central atom.
            previous = set(first_set)
            previous.add(ball_num)

            for first_num in previous:
                first_pos = num_to_pos.get(first_num)
                if first_pos is not None:
                    second_set.update(neighbor_sets[first_pos])

            second_set.difference_update(previous)

            if second_set:
                layer2_dists = []
                ball_loc_array = ball_locs_matrix[ball_pos]
                for ball_2 in second_set:
                    pos_2 = num_to_pos[ball_2]
                    delta = ball_locs_matrix[pos_2] - ball_loc_array
                    layer2_dists.append(sqrt(float(np.dot(delta, delta))))

                lyr2_dist_avg = sum(layer2_dists) / len(layer2_dists)
                nbor_dst_avgs[ball_pos].append(lyr2_dist_avg)
                nbor_lyr_rmsds[ball_pos].append(_rmsd(layer2_dists, lyr2_dist_avg))
            else:
                nbor_dst_avgs[ball_pos].append(0.0)
                nbor_lyr_rmsds[ball_pos].append(0.0)

            timer['neighbors_2'] += now() - t

        # If second neighbors were explicitly disabled, preserve the shape of
        # the historical complicated=False result (one first-layer value).

        # --------------------------------------------------------------
        # Center of mass / moment of inertia. Calculate both in one
        # tetrahedral traversal when both are requested.
        # --------------------------------------------------------------
        if com and moi:
            t = now()
            com_val, moi_val = calc_cell_mass_properties_cached(ball_loc, surf_ids, surf_points, surf_tris, volume)
            elapsed = now() - t
            coms[ball_pos] = com_val
            mois[ball_pos] = moi_val
            timer['com'] += elapsed

        elif com:
            t = now()
            com_val, _ = calc_cell_mass_properties_cached(ball_loc, surf_ids, surf_points, surf_tris, volume)
            coms[ball_pos] = com_val
            timer['com'] += now() - t

        elif moi:
            t = now()
            _, moi_val = calc_cell_mass_properties_cached(ball_loc, surf_ids, surf_points, surf_tris, volume)
            mois[ball_pos] = moi_val
            timer['moi'] += now() - t

        count += 1
        current_time = now()
        if current_time - last_update >= 0.25 or count == n_group:
            percentage = 100.0 * count / max(n_group, 1)
            net.update_progress("Analyzing network", percentage)
            last_update = current_time

    # ------------------------------------------------------------------
    # Bulk DataFrame assignments.
    # ------------------------------------------------------------------
    t = now()
    net.balls = net.balls.assign(
        vol=b_vols,
        sa=b_sas,
        max_mean_curv=b_max_mean_curvs,
        complete=b_cell,
        max_gauss_curv=b_max_gauss_curvs,
        avg_mean_surf_curv=b_avg_mean_surf_curvs,
        avg_gauss_surf_curv=b_avg_gauss_surf_curvs,
        int_mean_curv=b_int_mean_curvs,
        int_mean_curv_sq=b_int_mean_curv_sqs,
        int_gauss_curv=b_int_gauss_curvs,
        sphericity=b_sphrctys,
        isometric_quotient=b_isopmqs,
        ball_inside=b_inner,
        number_of_neighbors=num_nbors,
        nearest_neighbor=near_nbors,
        nearest_neighbor_distance=near_nbor_dists,
        neighbor_distance_average=nbor_dst_avgs,
        neighbor_distance_rmsd=nbor_lyr_rmsds,
        number_of_olaps=num_olaps,
        min_spike=b_min_spikes,
        max_spike=b_max_spikes,
        contact_area=contact_areas,
        olap_vol=olap_vols,
        vdw_vol=non_olap_vols,
        com=coms,
        moi=mois,
        bounding_box=b_boxs,
    )
    timer['ball_assign'] += now() - t

    t = now()
    net.surfs = net.surfs.assign(
        contact_area=surface_contact_area,
        overlap=surface_overlap,
    )
    timer['surface_assign'] += now() - t

    analysis_total = now() - analysis_start

    # Keep the existing top-level analysis metric calculation intact for
    # compatibility with the rest of VorPy's system-wide timing code.
    net.metrics['anal'] = (
        now()
        - net.metrics['start']
        - net.metrics['surf']
        - net.metrics['con']
        - net.metrics['vert']
    )

    # Detailed profiler is deliberately kept separate from the historical
    # scalar metrics so existing metric consumers are not forced to change.
    net.analysis_timing = timer.copy()
    net.analysis_timing['total'] = analysis_total

    _print_analysis_timing(timer, analysis_total)
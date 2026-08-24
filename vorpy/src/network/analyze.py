import time
from numpy import pi, sqrt
from vorpy.src.calculations import calc_sphericity
from vorpy.src.calculations import calc_isoperimetric_quotient
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import calc_spikes
from vorpy.src.calculations import calc_contacts
from vorpy.src.calculations import calc_cell_box
from vorpy.src.calculations import calc_cell_com
from vorpy.src.calculations import calc_cell_moi
from time import perf_counter as now


def append_0(*lists):
    """
    Appends a 0 to each list provided as an argument.
    """
    for my_list in lists:
        my_list.append(0)
    return lists


def get_next_layer(net, prev_layer):
    """
    Identifies and returns the next layer of balls based on the previous layer within a network.
    """
    layer2 = []
    prev_ndxs = [_['num'] for _ in prev_layer]
    for ball in prev_layer:
        ball_surfs = net.surfs.iloc[ball['surfs']].to_dict(orient='records')
        for surf in ball_surfs:
            other_ball = [_ for _ in surf['balls'] if _ != ball['num']][0]
            if other_ball not in layer2 and other_ball not in prev_ndxs:
                layer2.append(other_ball)
    return layer2


def analyze(net, complicated=True):
    """
    Performs a comprehensive analysis of a network, calculating various physical, geometrical,
    and topological properties of the cells (balls) within the network.
    """
    net.update_progress("Analyzing network", 0.0)
    # Precompute for speed
    group_set = set(net.group)
    n_group = len(group_set)

    # Set up the balls' volumes, surface areas, and completion variables
    b_vols, b_sas, b_cell = [], [], []

    # Set up the curvature variables
    b_max_mean_curvs, b_avg_mean_surf_curvs = [], []
    b_max_gauss_curvs, b_avg_gauss_surf_curvs = [], []

    # Set up the geometric variables
    b_sphrctys, b_isopmqs = [], []

    # Set up the neighbors variables
    num_nbors, near_nbors, near_nbor_dists = [], [], []
    nbor_lyr_rmsds, nbor_dst_avgs, b_inner = [], [], []

    # Set up the spike variables
    b_min_spikes, b_max_spikes = [], []

    # Set up the contacts lists
    contact_areas, non_olap_vols, olap_vols, num_olaps = [], [], [], []

    # Physical values
    coms, mois = [], []

    # Bounding boxes
    b_boxs = []

    # Timer
    timer = {
        'basic': 0, 'curvs': 0, 'geometric': 0,
        'nbors': 0, 'spikes': 0, 'contacts': 0,
        'com': 0, 'moi': 0, 'b_box': 0,
    }
    time_start = time.perf_counter()

    # Surfaces tracker
    surfaces_tracker = {}

    # Progress counter (only for balls actually analyzed)
    count = 0

    for k, ball in net.balls.iterrows():
        ball_num = ball['num']

        # FAST PATH: balls not in the group get zeros and are skipped early.
        # This is the key optimization when net.balls is large but net.group is small.
        if ball_num not in group_set:
            (b_vols, b_sas, b_cell,
             b_max_mean_curvs, b_avg_mean_surf_curvs,
             b_max_gauss_curvs, b_avg_gauss_surf_curvs,
             b_sphrctys, b_isopmqs, b_inner,
             num_nbors, near_nbors, near_nbor_dists,
             nbor_lyr_rmsds, num_olaps, nbor_dst_avgs,
             b_min_spikes, b_max_spikes,
             contact_areas, olap_vols, non_olap_vols,
             coms, mois, b_boxs) = append_0(
                b_vols, b_sas, b_cell,
                b_max_mean_curvs, b_avg_mean_surf_curvs,
                b_max_gauss_curvs, b_avg_gauss_surf_curvs,
                b_sphrctys, b_isopmqs, b_inner,
                num_nbors, near_nbors, near_nbor_dists,
                nbor_lyr_rmsds, num_olaps, nbor_dst_avgs,
                b_min_spikes, b_max_spikes,
                contact_areas, olap_vols, non_olap_vols,
                coms, mois, b_boxs,
            )
            continue

        count += 1

        if count == 1 or count % 100 == 0 or count == n_group:
            percentage = 100.0 * count / max(n_group, 1)
            net.update_progress("Analyzing network", percentage)

        # Get the ball's surfaces once
        surf_ids = ball['surfs']
        if not surf_ids:
            # No surfaces at all -> incomplete / zeroed
            (b_vols, b_sas, b_cell,
             b_max_mean_curvs, b_avg_mean_surf_curvs,
             b_max_gauss_curvs, b_avg_gauss_surf_curvs,
             b_sphrctys, b_isopmqs, b_inner,
             num_nbors, near_nbors, near_nbor_dists,
             nbor_lyr_rmsds, num_olaps, nbor_dst_avgs,
             b_min_spikes, b_max_spikes,
             contact_areas, olap_vols, non_olap_vols,
             coms, mois, b_boxs) = append_0(
                b_vols, b_sas, b_cell,
                b_max_mean_curvs, b_avg_mean_surf_curvs,
                b_max_gauss_curvs, b_avg_gauss_surf_curvs,
                b_sphrctys, b_isopmqs, b_inner,
                num_nbors, near_nbors, near_nbor_dists,
                nbor_lyr_rmsds, num_olaps, nbor_dst_avgs,
                b_min_spikes, b_max_spikes,
                contact_areas, olap_vols, non_olap_vols,
                coms, mois, b_boxs,
            )
            continue

        ball_surfs = net.surfs.iloc[surf_ids].to_dict(orient='records')

        # Quick test for pathological case: zero total surface area
        b_sa_list = [_['sa'] for _ in ball_surfs]
        sa_total = sum(b_sa_list)
        if sa_total == 0:
            (b_vols, b_sas, b_cell,
             b_max_mean_curvs, b_avg_mean_surf_curvs,
             b_max_gauss_curvs, b_avg_gauss_surf_curvs,
             b_sphrctys, b_isopmqs, b_inner,
             num_nbors, near_nbors, near_nbor_dists,
             nbor_lyr_rmsds, num_olaps, nbor_dst_avgs,
             b_min_spikes, b_max_spikes,
             contact_areas, olap_vols, non_olap_vols,
             coms, mois, b_boxs) = append_0(
                b_vols, b_sas, b_cell,
                b_max_mean_curvs, b_avg_mean_surf_curvs,
                b_max_gauss_curvs, b_avg_gauss_surf_curvs,
                b_sphrctys, b_isopmqs, b_inner,
                num_nbors, near_nbors, near_nbor_dists,
                nbor_lyr_rmsds, num_olaps, nbor_dst_avgs,
                b_min_spikes, b_max_spikes,
                contact_areas, olap_vols, non_olap_vols,
                coms, mois, b_boxs,
            )
            continue

        time1 = time.perf_counter()

        # Check for complete cells
        complete = True
        for vert in ball['verts']:
            # Count edges at this vertex that belong to this ball
            edge_balls_for_vert = [net.edges['balls'][e] for e in net.verts['edges'][vert]]
            if len([eb for eb in edge_balls_for_vert if k in eb]) != 3:
                if ball_num in group_set:
                    new_count = 0
                    for edge in ball['edges']:
                        if set(net.edges['balls'][edge]).issubset(set(net.verts['balls'][vert])):
                            new_count += 1
                    if new_count < 3:
                        complete = False
                else:
                    complete = False

        if len(ball['verts']) < 3 or len(ball['edges']) < 4 or len(ball_surfs) < 3:
            complete = False

        b_cell.append(complete)

        # Basic SA & volume
        sa = sa_total
        b_sas.append(sa)

        volume = sum(_['vols'][ball_num] for _ in ball_surfs)
        b_vols.append(volume)

        time2 = time.perf_counter()
        timer['basic'] += time2 - time1

        # Curvature metrics

        b_max_mean_curvs.append(max([_['mean_curv'] for _ in ball_surfs]))
        b_avg_mean_surf_curvs.append(
            sum(s['sa'] * s['avg_mean_curv'] for s in ball_surfs) / sa
        )

        b_max_gauss_curvs.append(max([_['gauss_curv'] for _ in ball_surfs]))
        b_avg_gauss_surf_curvs.append(
            sum(s['sa'] * s['avg_gauss_curv'] for s in ball_surfs) / sa
        )

        time3 = time.perf_counter()
        timer['curvs'] += time3 - time2

        # Geometric metrics
        b_sphrctys.append(calc_sphericity(volume=volume, surface_area=sa))
        b_isopmqs.append(calc_isoperimetric_quotient(volume=volume, surface_area=sa))

        time4 = time.perf_counter()
        timer['geometric'] += time4 - time3

        # Neighbor metrics
        neighbors, neighbors_nums, neighbor_dists = [], [], []
        ball_loc = ball['loc']
        ball_rad = ball['rad']

        for i, surf in enumerate(ball_surfs):
            neighbor_num = [_ for _ in surf['balls'] if _ != ball_num][0]
            neighbors_nums.append(neighbor_num)
            neighbor = net.balls.iloc[neighbor_num]
            neighbors.append(neighbor)

            neighbor_dist = calc_dist(ball_loc, neighbor['loc']) - ball_rad - neighbor['rad']
            neighbor_dists.append(neighbor_dist)

            surf_id = surf_ids[i]
            if surf_id not in surfaces_tracker:
                surfaces_tracker[surf_id] = {'olap_dist': max(-neighbor_dist, 0)}
            else:
                # Keep the max overlap distance per surface
                surfaces_tracker[surf_id]['olap_dist'] = max(
                    surfaces_tracker[surf_id]['olap_dist'],
                    max(-neighbor_dist, 0),
                )

        # Inner / outer
        b_inner.append(group_set.issuperset(neighbors_nums))

        # Neighbor stats
        num_nbors.append(len(neighbors))
        min_dist = min(neighbor_dists)
        near_nbor_dists.append(min_dist)
        near_nbors.append(neighbors[neighbor_dists.index(min_dist)]['num'])

        nbor_dist_avg = sum(neighbor_dists) / len(neighbor_dists)
        nbor_dst_avgs.append([nbor_dist_avg])
        nbor_lyr_rmsds.append([
            sqrt(sum((d - nbor_dist_avg) ** 2 for d in neighbor_dists) / len(neighbor_dists))
        ])

        time5 = time.perf_counter()
        timer['nbors'] += time5 - time4

        # Complicated / optional metrics
        if complicated:
            # Spikes
            min_spike, max_spike = calc_spikes(ball_loc, ball_surfs)
            b_min_spikes.append(min_spike)
            b_max_spikes.append(max_spike)

            time6 = time.perf_counter()
            timer['spikes'] += time6 - time5

            # Contacts
            contact_area, vdw_vol = calc_contacts(ball_loc, ball_rad, ball_surfs, surf_ids)
            num_olaps.append(len([d for d in neighbor_dists if d < 0]))
            contact_areas.append(sum(contact_area.values()))

            for surf_id in surf_ids:
                if surf_id not in surfaces_tracker:
                    surfaces_tracker[surf_id] = {'contact_area': contact_area[surf_id]}
                else:
                    surfaces_tracker[surf_id]['contact_area'] = contact_area[surf_id]

            non_olap_vols.append(vdw_vol)
            olap_vols.append((4.0 / 3.0) * pi * ball_rad ** 3 - vdw_vol)

            time6a = time.perf_counter()
            timer['contacts'] += time6a - time6

            # Second neighbor layer
            layer2 = get_next_layer(net, neighbors + [ball])
            layer2_dists = []
            for ball_2 in layer2:
                neighbor2 = net.balls.iloc[ball_2]
                layer2_dists.append(calc_dist(ball_loc, neighbor2['loc']))

            if layer2_dists:
                lyr2_dist_avg = sum(layer2_dists) / len(layer2_dists)
                nbor_dst_avgs[-1].append(lyr2_dist_avg)
                nbor_lyr_rmsds[-1].append(
                    sqrt(sum((d - lyr2_dist_avg) ** 2 for d in layer2_dists) / len(layer2_dists))
                )
            else:
                nbor_dst_avgs[-1].append(0.0)
                nbor_lyr_rmsds[-1].append(0.0)

            time7 = time.perf_counter()
            timer['nbors'] += time7 - time6a

            # Center of mass
            coms.append(calc_cell_com(ball_loc, ball_surfs, volume))
            time7a = time.perf_counter()
            timer['com'] += time7a - time7

            # Moment of inertia
            mois.append(calc_cell_moi(ball_loc, ball_surfs, volume))
            time8 = time.perf_counter()
            timer['moi'] += time8 - time7a

            # Bounding box
            b_boxs.append(calc_cell_box(ball_surfs))
            time9 = time.perf_counter()
            timer['b_box'] += time9 - time8

        else:
            # Cheap defaults when complicated=False
            b_min_spikes.append(0.0)
            b_max_spikes.append(0.0)
            contact_areas.append(0.0)
            non_olap_vols.append(0.0)
            olap_vols.append(0.0)
            num_olaps.append(0)
            coms.append(0.0)
            mois.append(0.0)
            b_boxs.append(0.0)

    # Assign the balls values
    net.balls = net.balls.assign(
        vol=b_vols, sa=b_sas, max_mean_curv=b_max_mean_curvs, complete=b_cell,
        max_gauss_curv=b_max_gauss_curvs, avg_mean_surf_curv=b_avg_mean_surf_curvs,
        avg_gauss_surf_curv=b_avg_gauss_surf_curvs, sphericity=b_sphrctys,
        isometric_quotient=b_isopmqs, ball_inside=b_inner, number_of_neighbors=num_nbors,
        nearest_neighbor=near_nbors, nearest_neighbor_distance=near_nbor_dists,
        neighbor_distance_average=nbor_dst_avgs, neighbor_distance_rmsd=nbor_lyr_rmsds,
        number_of_olaps=num_olaps, min_spike=b_min_spikes, max_spike=b_max_spikes,
        contact_area=contact_areas, olap_vol=olap_vols, vdw_vol=non_olap_vols,
        com=coms, moi=mois, bounding_box=b_boxs,
    )

    # Initialize surface columns
    net.surfs = net.surfs.assign(
        contact_area=[0.0 for _ in range(len(net.surfs))],
        overlap=[0.0 for _ in range(len(net.surfs))],
    )
    for surf_id, info in surfaces_tracker.items():
        if 'contact_area' in info:
            net.surfs.loc[surf_id, 'contact_area'] = info['contact_area']
        if 'olap_dist' in info:
            net.surfs.loc[surf_id, 'overlap'] = info['olap_dist']

    net.metrics['anal'] = (
        now()
        - net.metrics['start']
        - net.metrics['surf']
        - net.metrics['con']
        - net.metrics['vert']
    )

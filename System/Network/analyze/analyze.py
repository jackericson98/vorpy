import time

from numpy import pi, sqrt
from System.sys_funcs.calcs.calcs import (calc_sphericity, get_time, calc_isoperimetric_quotient, calc_dist, calc_spikes,
                                          calc_contacts, calc_cell_box, calc_cell_com, calc_cell_moi)
from time import perf_counter as now


def append_0(*lists):
    for my_list in lists:
        my_list.append(0)
    return lists


def get_next_layer(net, prev_layer):
    # Set up the layer2 list
    layer2 = []
    # Gather the ball numbers for the previous layer
    prev_ndxs = [_['num'] for _ in prev_layer]
    # Loop through the current layer
    for ball in prev_layer:
        ball_surfs = net.surfs.iloc[ball['surfs']].to_dict(orient='records')
        for surf in ball_surfs:
            other_ball = [_ for _ in surf['balls'] if _ != ball['num']][0]
            if other_ball not in layer2 and other_ball not in prev_ndxs:
                layer2.append(other_ball)
    return layer2


def analyze(net, complicated=True):
    """
     Analyzes the output surfaces, cells and solute vertices for the network for later reference
     """
    # Create the group set for the group indexes
    group_set = set(net.group)
    # Set up the balls' volumes, surface areas, and completion variables
    b_vols, b_sas, b_cell = [], [], []
    # Set up the curvature variables
    b_max_curvs, b_avg_surf_curvs = [], []
    # Set up the geometric variables
    b_sphrctys, b_isopmqs = [], []
    # Set up the neighbors variables
    num_nbors, near_nbors, near_nbor_dists, nbor_lyr_rmsds, nbor_dst_avgs, b_inner = [], [], [], [], [], []
    # Set up the spike variables
    b_min_spikes, b_max_spikes = [], []
    # Set up the contacts lists
    contact_areas, non_olap_vols, olap_vols, num_olaps = [], [], [], []
    # Physical values
    coms, mois = [], []
    # Bounding boxes
    b_boxs = []
    # Set up the timer
    timer = {'basic': 0, 'curvs': 0, 'geometric': 0, 'nbors': 0, 'spikes': 0, 'contacts': 0, 'com': 0, 'moi': 0, 'b_box': 0}
    time_start = time.perf_counter()
    # Set up the surfaces tracker
    surfaces_tracker = {}
    # Go through each ball in the system and find the volume
    count = 0
    for k, ball in net.balls.iterrows():

        # Get the percentage for printing
        percentage = int(count / len(net.group) * 100)
        # Print the actions
        my_time = now() - net.metrics['start']
        h, m, s = get_time(my_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                 "
              .format(int(h), int(m), round(s, 2), percentage), end="")

        # Get the ball surfs
        ball_surfs = net.surfs.iloc[ball['surfs']].to_dict(orient='records')

        # Initial test for completeness
        if len(ball['surfs']) == 0 or sum([_['sa'] for _ in ball_surfs]) == 0 or ball['num'] not in group_set:
            (b_vols, b_sas, b_cell, b_max_curvs, b_avg_surf_curvs, b_sphrctys, b_isopmqs, b_inner, num_nbors,
             near_nbors, near_nbor_dists, nbor_lyr_rmsds, num_olaps, nbor_dst_avgs, b_min_spikes, b_max_spikes,
             contact_areas, olap_vols, non_olap_vols, coms, mois, b_boxs) = (
                append_0(b_vols, b_sas, b_cell, b_max_curvs, b_avg_surf_curvs, b_sphrctys, b_isopmqs, b_inner,
                         num_nbors, near_nbors, near_nbor_dists, nbor_lyr_rmsds, num_olaps, nbor_dst_avgs, b_min_spikes,
                         b_max_spikes, contact_areas, olap_vols, non_olap_vols, coms, mois, b_boxs))
            continue
        count += 1
        time1 = time.perf_counter()
        # Check for complete cells in the balls
        complete = True
        # Go through each of the vertices in the ball
        for vert in ball['verts']:
            # Check the number of edges from the vertex that hold
            if len([_ for _ in [net.edges['balls'][__] for __ in net.verts['edges'][vert]] if k in _]) != 3:
                complete = False
        # Additional catch for any ball that doesn't have the 181L number of network elements associated with it
        if len(ball['verts']) < 3 or len(ball['edges']) < 4 or len(ball_surfs) < 3:
            complete = False
        # Add the complete designation for the cell
        b_cell.append(complete)

        # Calculate the surface area of the ball by summing the surface areas of all it's surfaces
        sa = sum([_['sa'] for _ in ball_surfs])
        b_sas.append(sa)

        # Calculate the volume of the ball by the previously stored volume data
        volume = sum([_['vols'][ball['num']] for _ in ball_surfs])
        b_vols.append(volume)

        time2 = time.perf_counter()
        timer['basic'] += time2 - time1

        # Go through the ball's surfaces
        b_max_curvs.append(max([_['curv'] for _ in ball_surfs]))
        b_avg_surf_curvs.append(sum(_['sa'] * _['curv'] for _ in ball_surfs) / sa)

        time3 = time.perf_counter()
        timer['curvs'] += time3 - time2

        # Calculate the sphericity
        b_sphrctys.append(calc_sphericity(volume=volume, surface_area=sa))

        # Calculate the isoperimetric quotient
        b_isopmqs.append(calc_isoperimetric_quotient(volume=volume, surface_area=sa))

        time4 = time.perf_counter()
        timer['geometric'] += time4 - time3

        # Gather the neighbors
        neighbors, neighbors_nums, neighbor_dists = [], [], []
        for i, surf in enumerate(ball_surfs):
            neighbor_num = [_ for _ in surf['balls'] if _ != ball['num']][0]

            neighbors_nums.append(neighbor_num)
            neighbor = net.balls.iloc[neighbor_num]
            neighbor_dist = calc_dist(ball['loc'], neighbor['loc']) - ball['rad'] - neighbor['rad']
            neighbor_dists.append(neighbor_dist)
            neighbors.append(neighbor)
            if ball['surfs'][i] not in surfaces_tracker:
                surfaces_tracker[ball['surfs'][i]] = {'olap_dist': max(-neighbor_dist, 0)}
        # Check to see if the ball is inside or outside
        if group_set.issuperset(neighbors_nums):
            b_inner.append(True)
        else:
            b_inner.append(False)

        # Add the number of neighbors and the nearest neighbor
        num_nbors.append(len(neighbors))
        near_nbor_dists.append(min(neighbor_dists))
        near_nbors.append(neighbors[neighbor_dists.index(near_nbor_dists[-1])]['num'])
        nbor_dist_avg = sum(neighbor_dists) / len(neighbor_dists)
        nbor_dst_avgs.append([nbor_dist_avg])
        nbor_lyr_rmsds.append([sqrt(sum([(_ - nbor_dist_avg) ** 2 for _ in neighbor_dists]) / len(neighbor_dists))])

        time5 = time.perf_counter()
        timer['nbors'] += time5 - time4

        # The more complicated/time_consuming calculations happen here
        if complicated:
            # Add the spike variables
            min_spike, max_spike = calc_spikes(ball['loc'], ball_surfs)
            # Add them to the respective lists
            b_min_spikes.append(min_spike)
            b_max_spikes.append(max_spike)

            time6 = time.perf_counter()
            timer['spikes'] += time6 - time5

            # Get the contact information

            contact_area, vdw_vol = calc_contacts(ball['loc'], ball['rad'], ball_surfs, ball['surfs'])
            num_olaps.append(len([_ for _ in neighbor_dists if _ < 0]))
            contact_areas.append(sum([contact_area[_] for _ in contact_area]))
            for i, surf in enumerate(ball_surfs):
                surface_key = ball['surfs'][i]
                if surface_key not in surfaces_tracker:
                    surfaces_tracker[surface_key] = {'contact_area': contact_area[surface_key]}
                else:
                    surfaces_tracker[surface_key]['contact_area'] = contact_area[surface_key]
            non_olap_vols.append(vdw_vol)
            olap_vols.append((4/3) * pi * ball['rad'] ** 3 - vdw_vol)

            time6a = time.perf_counter()
            timer['contacts'] += time6a - time6

            # Get the additional layers
            layer2 = get_next_layer(net, neighbors + [ball])
            layer2_balls, layer2_dists = [], []
            # get the distances from the ball for the next layer
            for ball_2 in layer2:
                neighbor2 = net.balls.iloc[ball_2]
                layer2_dists.append(calc_dist(ball['loc'], neighbor2['loc']))
                layer2_balls.append(neighbor2)
            if len(layer2_dists) > 0:
                lyr2_dist_avg = sum(layer2_dists) / len(layer2_dists)
                nbor_dst_avgs[-1].append(lyr2_dist_avg)
                nbor_lyr_rmsds[-1].append(sqrt(sum([(_ - lyr2_dist_avg) ** 2 for _ in layer2_dists]) / len(layer2_dists)))
            else:
                nbor_dst_avgs[-1].append(0)
                nbor_lyr_rmsds[-1].append(0)

            time7 = time.perf_counter()
            timer['nbors'] += time7 - time6a

            # get the center of mass
            coms.append(calc_cell_com(ball['loc'], ball_surfs, volume))

            time7a = time.perf_counter()
            timer['com'] += time7a - time7
            # Get the Moment of inertia
            mois.append(calc_cell_moi(ball['loc'], ball_surfs, volume))

            time8 = time.perf_counter()
            timer['moi'] += time8 - time7a

            # Get the bounding box
            b_boxs.append(calc_cell_box(ball_surfs))

            time9 = time.perf_counter()
            timer['b_box'] += time9 - time8
        else:
            b_min_spikes.append(0)
            b_max_spikes.append(0)
            contact_areas.append(0)
            non_olap_vols.append(0)
            coms.append(0)
            mois.append(0)
            b_boxs.append(0)
    # Assign the balls values
    net.balls = net.balls.assign(vol=b_vols, sa=b_sas, max_curv=b_max_curvs, complete=b_cell,
                                 avg_surf_curv=b_avg_surf_curvs, sphericity=b_sphrctys, isometric_quotient=b_isopmqs,
                                 ball_inside=b_inner, number_of_neighbors=num_nbors, nearest_neighbor=near_nbors,
                                 nearest_neighbor_distance=near_nbor_dists, neighbor_distance_average=nbor_dst_avgs,
                                 neighbor_distance_rmsd=nbor_lyr_rmsds, number_of_olaps=num_olaps,
                                 min_spike=b_min_spikes, max_spike=b_max_spikes, contact_area=contact_areas,
                                 non_olap_vol=non_olap_vols, vdw_vol=olap_vols, com=coms, moi=mois, bounding_box=b_boxs)
    # First make the surfaces columns for contact area and overlap volume
    net.surfs = net.surfs.assign(contact_area=[0.0 for _ in range(len(net.surfs))],
                                 overlap=[0.0 for _ in range(len(net.surfs))])
    for surf in surfaces_tracker:
        net.surfs.loc[surf, 'contact_area'] = surfaces_tracker[surf]['contact_area']
        net.surfs.loc[surf, 'overlap'] = surfaces_tracker[surf]['olap_dist']

    for _ in timer:
        print(_, timer[_])
    print('total = {} s\n'.format(time.perf_counter() - time_start))

    net.metrics['anal'] = now() - net.metrics['start'] - net.metrics['surf'] - net.metrics['con'] - net.metrics['vert']

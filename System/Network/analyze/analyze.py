import time

from numpy import pi, sqrt
from System.sys_funcs.calcs.calcs import (calc_sphericity, get_time, calc_isoperimetric_quotient, calc_dist, calc_spikes,
                                          calc_contacts, calc_cell_box, calc_cell_com, calc_cell_moi)
from time import perf_counter as now


def append_0(*lists):
    for my_list in lists:
        my_list.append(0)
    return lists


def analyze(net, complicated=True):
    """
     Analyzes the output surfaces, cells and solute vertices for the network for later reference
     """
    # Set up the balls' volumes, surface areas, and completion variables
    b_vols, b_sas, b_cell = [], [], []
    # Set up the curvature variables
    b_max_curvs, b_avg_surf_curvs = [], []
    # Set up the geometric variables
    b_sphrctys, b_isopmqs = [], []
    # Set up the neighbors variables
    num_nbors, near_nbors, near_nbor_dists, nbor_dst_rmsds, nbor_dst_avgs = [], [], [], [], []
    # Set up the spike variables
    b_min_spikes, b_max_spikes = [], []
    # Set up the contacts lists
    contact_areas, non_olap_vols, olap_vols, num_olaps = [], [], [], []
    # Physical values
    coms, mois = [], []
    # Bounding boxes
    b_boxs = []
    # Set up the timer
    timer = {'basic': 0, 'curvs': 0, 'geometric': 0, 'nbors': 0, 'spikes': 0, 'contacts': 0, 'physical': 0, 'b_box': 0}
    # Go through each ball in the system and find the volume
    for k, ball in net.balls.iterrows():

        # Get the percentage for printing
        percentage = int(k / len(net.balls['loc']) * 100)
        # Print the actions
        my_time = now() - net.metrics['start']
        h, m, s = get_time(my_time)
        print("\rRun Time = {}:{}:{:.2f} - Process: analyzing: {} %                 "
              .format(int(h), int(m), round(s, 2), percentage), end="")

        # Get the ball surfs
        ball_surfs = net.surfs.iloc[ball['surfs']].to_dict(orient='records')

        # Initial test for completeness
        if len(ball['surfs']) == 0:
            (b_vols, b_sas, b_cell, b_max_curvs, b_avg_surf_curvs, b_sphrctys, b_isopmqs, num_nbors, near_nbors,
             near_nbor_dists, nbor_dst_rmsds, num_olaps, nbor_dst_avgs, b_min_spikes, b_max_spikes, contact_areas,
             olap_vols, non_olap_vols, coms, mois, b_boxs) = (
                append_0(b_vols, b_sas, b_cell, b_max_curvs, b_avg_surf_curvs, b_sphrctys, b_isopmqs, num_nbors,
                         near_nbors, near_nbor_dists, nbor_dst_rmsds, num_olaps, nbor_dst_avgs, b_min_spikes,
                         b_max_spikes, contact_areas, olap_vols, non_olap_vols, coms, mois, b_boxs))
            continue

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
        neighbors, neighbor_dists = [], []
        for surf in ball_surfs:
            neighbor = net.balls.iloc[[_ for _ in surf['balls'] if _ != ball['num']][0]]
            neighbor_dists.append(calc_dist(ball['loc'], neighbor['loc']) - ball['rad'] - neighbor['rad'])
            neighbors.append(neighbor)

        # Add the number of neighbors and the nearest neighbor
        num_nbors.append(len(neighbors))
        near_nbor_dists.append(min(neighbor_dists))
        near_nbors.append(neighbors[neighbor_dists.index(near_nbor_dists[-1])])
        nbor_dist_avg = sum(neighbor_dists) / len(neighbor_dists)
        nbor_dst_avgs.append(nbor_dist_avg)
        nbor_dst_rmsds.append(sqrt(sum([(_ - nbor_dist_avg) ** 2 for _ in neighbor_dists]) / len(neighbor_dists)))

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
            contact_area, vdw_vol = calc_contacts(ball['loc'], ball['rad'], ball_surfs)
            num_olaps.append(len(contact_area))
            contact_areas.append(sum(contact_area))
            non_olap_vols.append(vdw_vol)
            olap_vols.append((4/3) * pi * ball['rad'] ** 3 - vdw_vol)

            time7 = time.perf_counter()
            timer['contacts'] += time7 - time6

            # get the center of mass
            coms.append(calc_cell_com(ball['loc'], ball_surfs, volume))
            # Get the Moment of inertia
            mois.append(calc_cell_moi(ball['loc'], ball_surfs, volume))

            time8 = time.perf_counter()
            timer['physical'] += time8 - time7

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

    net.balls = net.balls.assign(vol=b_vols, sa=b_sas, max_curv=b_max_curvs, complete=b_cell,
                                 avg_surf_curv=b_avg_surf_curvs, sphericity=b_sphrctys, isometric_quotient=b_isopmqs,
                                 neighbor_number=num_nbors, near_neighbor=near_nbors,
                                 near_neighbor_distance=near_nbor_dists, neighbor_distance_average=nbor_dst_avgs,
                                 neighbor_distance_rmsd=nbor_dst_rmsds, num_olaps=num_olaps, min_spike=b_min_spikes,
                                 max_spike=b_max_spikes, contact_area=contact_areas, non_olap_vol=non_olap_vols,
                                 vdw_vol=olap_vols, com=coms, moi=mois, b_box=b_boxs)

    for _ in timer:
        print(_, timer[_])

    net.metrics['anal'] = now() - net.metrics['start'] - net.metrics['surf'] - net.metrics['con'] - net.metrics['vert']

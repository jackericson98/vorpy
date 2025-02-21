import time
import pandas as pd
from System.sys_funcs.calcs.sorting import global_vars, box_search, get_balls
from System.sys_funcs.calcs.calcs import calc_dist
from System.Network.verts.find_verts import find_verts
from System.sys_funcs.output.net import write_verts


def find_net_verts(net):

    # Not sure what this does
    # vert_list_real = net.get_real_verts()
    # Create the group indices
    if net.group is None:
        net.group = [_['num'] for i, _ in net.balls.iterrows()]
    sphere_check_list = net.group.copy()
    # Get the indices of the balls in the network to keep track of the balls that haven't been visited
    my_guuy = find_verts(locs=net.balls['loc'].to_numpy(), rads=net.balls['rad'].to_numpy(),
                         max_vert=net.settings['max_vert'], net_type=net.settings['net_type'], check_ndxs=sphere_check_list,
                         my_group=net.group, start_time=net.metrics['start'], print_metrics=net.settings['print_metrics'],
                         vert_box=net.settings['foam_box'], box=net.box['verts'])
    if my_guuy is not None:
        vert_ndxs, vlocs, vrads, vloc2s, vrad2s, sphere_check_list, averts = my_guuy
    # Check to see if any of the balls are encapsulated
    if len(sphere_check_list) > 0:
        skip_nums = []
        for sphere in sphere_check_list:
            sphere_rad, sphere_loc = net.balls['rad'][sphere], net.balls['loc'][sphere]
            sphere_box = box_search(sphere_loc)
            close_spheres = get_balls(sphere_box, dist=max(net.balls['rad']) - sphere_rad)
            for sphere2 in close_spheres:
                if calc_dist(sphere_loc, net.balls['loc'][sphere2]) < abs(net.balls['rad'][sphere2] - sphere_rad):
                    print("\nUh oh! Ball # {} is fully encapsulated by ball # {}! Skipping {}"
                          .format(sphere, sphere2, sphere))
                    skip_nums.append(sphere)
                    break
        for _ in skip_nums:
            sphere_check_list.pop(sphere_check_list.index(_))
    # Check for disconnects in the network
    while len(sphere_check_list) > 0:
        a0 = sphere_check_list.pop()
        my_guuy = find_verts(b0=a0, locs=net.balls['loc'].to_numpy(), rads=net.balls['rad'].to_numpy(),
                             max_vert=net.settings['max_vert'], net_type=net.settings['net_type'], check_ndxs=sphere_check_list,
                             my_group=net.group, vert_ndxs=vert_ndxs, vlocs=vlocs, vrads=vrads,
                             vloc2s=vloc2s, vrad2s=vrad2s, start_time=net.metrics['start'],
                             vert_box=net.settings['foam_box'], b_verts=averts, box=net.box['verts'])
        if my_guuy is not None:
            vert_ndxs, vlocs, vrads, vloc2s, vrad2s, sphere_check_list, averts = my_guuy
        if net.settings['ball_type'] == 'foam' and len(sphere_check_list) <= 0.25 * len(net.balls['loc']):
            print(f'Missing Ball Indices:\n{sphere_check_list}\n')
            break
    # # Create the doublets list
    # if vert_list_real is not None and net.type == 'aw':
    #     missing_verts = [_ for _ in vert_list_real if _ not in vert_ndxs]
    #     print(missing_verts)
    #     extra_verts = [_ for _ in vert_ndxs if _ not in vert_list_real]
    #     print(extra_verts)
    doublets = [0 for _ in range(len(vert_ndxs))]
    # Incorporate the doublets into the v_locs, balls, v_rads lists and lose the v_loc2s and v_rad2s
    i = 0
    while i < len(vlocs):
        # Check for doubletness
        if vrad2s[i] is not None:
            # Insert the relevant information into their respective lists
            vert_ndxs.insert(i + 1, vert_ndxs[i])
            vlocs.insert(i + 1, vloc2s[i])
            vrads.insert(i + 1, vrad2s[i])
            doublets[i] = 2
            doublets.insert(i + 1, 1)
            # Preserve the relational aspects of vrad2s and vloc2s
            vrad2s.insert(i + 1, None)
            vloc2s.insert(i + 1, [None, None, None])
        i += 1

    # Make the dataframe
    net.verts = pd.DataFrame({"balls": vert_ndxs, 'loc': vlocs, 'rad': vrads, 'dub': doublets})
    # Clear the print statement
    print("\r                                                                  ", end="")
    net.metrics['vert'] = time.perf_counter() - net.metrics['start']
    write_verts(net)

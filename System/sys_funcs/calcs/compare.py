import os
import csv
import time
from os import path
from System.sys_objs.interface import Interface
from System.sys_funcs.calcs.calcs import calc_dist


def compare_networks(sys, group1, group2, data_file=None):
    """
    The goal is to take the comparison instructions and make two separate groups with their networks and compare
    their results
    """
    start = time.perf_counter()
    # Create the data storage
    data = {'vdn1': [], 'sdn1': [], 'vdn2': [], 'sdn2': [], 'rads': []}
    com_counter = 0
    # Compare the networks
    for i, ball1 in group1.net.balls.iterrows():
        # Get the equivalent ball from the second group
        ball2 = group2.net.balls.iloc[i]
        # Make sure both cells are complete
        if ball1['complete'] and ball2['complete']:
            com_counter += 1
            # Calculate the differences in volume and surface area for each network as the standard
            vdn1, sdn1, vdn2, sdn2, rads = ((ball2['vol'] - ball1['vol']) / ball1['vol'],
                                            (ball2['sa'] - ball1['sa']) / ball1['sa'],
                                            (ball1['vol'] - ball2['vol']) / ball2['vol'],
                                            (ball1['sa'] - ball2['sa']) / ball2['sa'], ball1['rad'])
            # Check for outliers
            if any([_ > 25 for _ in [vdn1, sdn1, vdn2, sdn2]]):
                print('Outlier in comparison detected: {} - Off by {} %'.format(ball1['name'], 100 * vdn1))
                continue

            overlaps = []
            for surf in ball1['surfs']:
                surfster = group1.net.surfs.iloc[surf]
                neighbor = group1.net.balls.iloc[[_ for _ in surfster['balls'] if _ != ball1['num']]].to_dict(orient='records')[0]
                # print(neighbor)
                overlap_distance = calc_dist(ball1['loc'], neighbor['loc']) - ball1['rad'] - neighbor['rad']
                if overlap_distance < 0:
                    percenty = abs(overlap_distance) / min(neighbor['rad'], ball1['rad'])
                else:
                    percenty = 0.0
                overlaps.append(percenty)
            cwd = os.getcwd()
            os.chdir(sys.files['dir'])
            os.chdir('..')

            with open(os.getcwd() + '/overlaps.csv', 'a') as poopster_mccalister:
                livvydunne = csv.writer(poopster_mccalister)
                livvydunne.writerow([sys.files['dir'], ball1['num']] + overlaps)
            os.chdir(cwd)

            # Record the overlaps per ball
            # Add the data
            data['vdn1'].append(vdn1)
            data['sdn1'].append(sdn2)
            data['vdn2'].append(vdn2)
            data['sdn2'].append(sdn2)
            data['rads'].append(ball1['rad'])

    # Create the data line to be added to the data file
    nbs, my_line = len(data['vdn1']), []
    if sys.foam_data is None:
        sys.foam_data = []
    if nbs > 0:
        my_line = ("\r{}".format(sys.files['dir']), *sys.foam_data,
                   round(sum([abs(_) for _ in data['vdn1']]) / nbs, 5),  # Mean absolute difference
                   round(sum([abs(_) for _ in data['sdn1']]) / nbs, 5),  # Mean absolute difference
                   round(sum([abs(_) for _ in data['vdn2']]) / nbs, 5),  # Mean absolute difference
                   round(sum([abs(_) for _ in data['sdn2']]) / nbs, 5),  # Mean absolute difference
                   round(sum(data['vdn1']) / nbs, 5),  # Percent Difference
                   round(sum(data['sdn1']) / nbs, 5),  # Percent Difference
                   round(sum(data['vdn2']) / nbs, 5),  # Percent Difference
                   round(sum(data['sdn2']) / nbs, 5),  # Percent Difference
                   # round(np.polyfit(data['rads'], data['vdn1'], 1)[0], 5),  # Slope of the val by radius
                   # round(np.polyfit(data['rads'], data['sdn1'], 1)[0], 5),  # Slope of the val by radius
                   # round(np.polyfit(data['rads'], data['vdn2'], 1)[0], 5),  # Slope of the val by radius
                   # round(np.polyfit(data['rads'], data['sdn2'], 1)[0], 5),  # Slope of the val by radius
                   nbs, round((time.perf_counter() - sys.start), 3), com_counter, group1.settings['max_vert'])
    print(*my_line, end="")
    print('\n')

    # Make the data file location
    if data_file is None or not path.exists(data_file):
        cwd = os.getcwd()
        os.chdir(sys.files['dir'])
        os.chdir('..')
        data_file = os.getcwd() + '/foam_data.csv'
        os.chdir(cwd)
        # data_file = self.files['root_dir'] + '/Data/user_data/foam_data.csv'

    try:
        with open(data_file, 'a') as foam_file:
            foam_writer = csv.writer(foam_file)
            foam_writer.writerow(my_line)
    except PermissionError:
        with open(data_file[:-4] + '1.csv', 'a') as foam_file:
            foam_writer = csv.writer(foam_file)
            foam_writer.writerow(my_line)


def make_interfaces(sys):
    """
    We want to go through the groups in the system and see if an interface exists and if it does gather the materials
    from the networks of the groups to make interface networks
    """
    # First make sure that there is at least two groups in the system
    if len(sys.groups) < 2:
        return
    # Instantiate the interfaces attribute
    if sys.ifaces is None:
        sys.ifaces = []
    # Group1s that have been made tracker for not doing reverse
    group1_trackers = []
    # Loop through the groups in the system
    for group1 in sys.groups:
        # Loop through the groups again
        for group2 in sys.groups:
            # Skip when the groups are the same or when the balls are the same
            if group1 == group2 or group1.ball_ndxs == group2.ball_ndxs or group2 in group1_trackers:
                continue

            # Check that there are no overlapping ball ndxs
            olap_ndxs = []
            for ball_ndx in group1.ball_ndxs:
                if ball_ndx in group2.ball_ndxs:
                    olap_ndxs.append(ball_ndx)

            # Create a set out of the group2. ball ndxs
            g2_bndxs = set(group2.ball_ndxs)
            # Get the overlapping surfaces
            possible_surfs = group1.net.surfs[group1.net.surfs['balls'].apply(lambda balls: any(ball in g2_bndxs for ball in balls))]
            # Check that there are any overlapping surfaces at all
            if len(possible_surfs) == 0:
                continue
            # Finally add the Interface to the system's list of interfaces
            sys.ifaces.append(Interface(group1, group2, surfs=possible_surfs))
        # Add the group to group1 trackers
        group1_trackers.append(group1)


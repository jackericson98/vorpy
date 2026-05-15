import os
import csv
import time
from os import path
from vorpy.src.calculations.calcs import calc_dist


def compare_networks(sys, group1, group2, data_file=None):
    """
    Compares two network groups by analyzing their structural and geometric properties.

    This function performs a detailed comparison between two network groups by:
    1. Calculating relative differences in volume and surface area between corresponding balls
    2. Identifying and filtering out statistical outliers
    3. Computing overlap metrics between neighboring balls
    4. Recording comparison data for analysis

    Parameters:
    -----------
    sys : object
        System object containing file paths and configuration
    group1 : object
        First network group to compare
    group2 : object
        Second network group to compare
    data_file : str, optional
        Path to file where comparison data will be stored

    Returns:
    --------
    dict
        Dictionary containing comparison metrics:
        - 'vdn1': Volume differences normalized by group1 volumes
        - 'sdn1': Surface area differences normalized by group1 surface areas
        - 'vdn2': Volume differences normalized by group2 volumes
        - 'sdn2': Surface area differences normalized by group2 surface areas
        - 'rads': Radii of the balls being compared
    """
    # Start the timer
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

            # Calculate overlap metrics
            overlaps = []
            # Loop through the surfaces of the first ball
            for surf in ball1['surfs']:
                # Get the surface data
                surfster = group1.net.surfs.iloc[surf]
                # Get the neighbor ball
                neighbor = group1.net.balls.iloc[[_ for _ in surfster['balls'] if _ != ball1['num']]].to_dict(orient='records')[0]
                # Calculate the overlap distance
                overlap_distance = calc_dist(ball1['loc'], neighbor['loc']) - ball1['rad'] - neighbor['rad']
                # Calculate the overlap percentage
                if overlap_distance < 0:
                    # Calculate the percentage of overlap
                    percenty = abs(overlap_distance) / min(neighbor['rad'], ball1['rad'])
                else:
                    # No overlap
                    percenty = 0.0
                overlaps.append(percenty)
            # Make the data file location
            cwd = os.getcwd()
            os.chdir(sys.files['dir'])
            os.chdir('..')
            # Write the overlaps to the data file
            with open(os.getcwd() + '/overlaps.csv', 'a') as overlaps_file:
                # Create the writer
                overlaps_writer = csv.writer(overlaps_file)
                # Write the data
                overlaps_writer.writerow([sys.files['dir'], ball1['num']] + overlaps)
            # Change back to the original directory
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
    # If the foam data is not None
    if sys.foam_data is None:
        # Set the foam data to an empty list
        sys.foam_data = []
    # If there are data to add
    if nbs > 0:
        # Create the data line
        my_line = ("\r{}".format(sys.files['dir']), *sys.foam_data,
                   round(sum([abs(_) for _ in data['vdn1']]) / nbs, 5),  # Mean absolute difference
                   round(sum([abs(_) for _ in data['sdn1']]) / nbs, 5),  # Mean absolute difference
                   round(sum([abs(_) for _ in data['vdn2']]) / nbs, 5),  # Mean absolute difference
                   round(sum([abs(_) for _ in data['sdn2']]) / nbs, 5),  # Mean absolute difference
                   round(sum(data['vdn1']) / nbs, 5),  # Percent Difference
                   round(sum(data['sdn1']) / nbs, 5),  # Percent Difference
                   round(sum(data['vdn2']) / nbs, 5),  # Percent Difference
                   round(sum(data['sdn2']) / nbs, 5),  # Percent Difference
                   nbs, round((time.perf_counter() - sys.start), 3), com_counter, group1.settings['max_vert'])
        # Print the data line
        print(*my_line, end="")
        # Print a new line
        print('\n')

    # Make the data file location
    if data_file is None or not path.exists(data_file):
        # Make the data file location
        cwd = os.getcwd()
        os.chdir(sys.files['dir'])
        os.chdir('..')
        data_file = os.getcwd() + '/foam_data.csv'
        # Change back to the original directory
        os.chdir(cwd)

    # Try to open the data file
    try:
        # Open the data file    
        with open(data_file, 'a') as foam_file:
            # Create the writer
            foam_writer = csv.writer(foam_file)
            # Write the data
            foam_writer.writerow(my_line)
    # If there is a permission error
    except PermissionError:
        # Open the data file
        with open(data_file[:-4] + '1.csv', 'a') as foam_file:
            # Create the writer
            foam_writer = csv.writer(foam_file)
            # Write the data
            foam_writer.writerow(my_line)

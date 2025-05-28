import numpy as np
import matplotlib.pyplot as plt


import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.analyze.Part_I_Statistical_Ensembles.Review1.NeighborsAverageRadii import get_data, select_file

"""
Splits the main ball containers into subsections and calculates the local density and polydispersity of each and compares them to sphericity and volume difference.

We are looking to compare a ball's local polydispersity, local density, radius, sphericity, and volume difference

"""

def sort_balls(box_size, balls, num_splits=10):
    """
    Splits the main ball containers into subsections and calculates the local density and polydispersity of each and compares them to sphericity and volume difference.
    """
    # Sub box length
    sub_box_length = box_size / num_splits
    # Create the sub boxes dictionary which are identified by their tuple of indices
    sub_boxes = {(i, j, k): [] for i in range(num_splits) for j in range(num_splits) for k in range(num_splits)}
    # Loop through the balls and sort them into the correct subsection for polydispersity
    for ball in balls:
        # Check to make sure the ball is in the original 1000
        if ball >= 1000:
            continue
        # Find the index of the sub box the ball is in 
        sub_box_index = tuple([int(np.floor(coord / sub_box_length)) for coord in balls[ball]['loc']])
        # Add the ball to the correct sub box
        sub_boxes[sub_box_index].append(ball)
    # Return the sub boxes
    return sub_boxes


def calculate_local_density(balls, sub_boxes, sub_box_length, num_probe_points=1000):
    """
    Calculates the local density of each sub box and assigns the value to each ball. 
    """
    # Create a dictionary to store the local densities
    local_densities = {}
    # Loop through the sub boxes and get the surrounding boxes' balls 
    for sub_box in sub_boxes:
        # Get the surrounding sub boxes
        surrounding_boxes = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    # Skip the current box
                    if i == 0 and j == 0 and k == 0:
                        continue
                    # Get the surrounding box index
                    surrounding_box = (sub_box[0] + i, sub_box[1] + j, sub_box[2] + k)
                    # Check if the surrounding box exists
                    if surrounding_box in sub_boxes:
                        surrounding_boxes.append(surrounding_box)
        
        # Gather all balls that could overlap with this subbox
        all_balls = []
        # Add balls from current subbox
        all_balls.extend(sub_boxes[sub_box])
        # Add balls from surrounding boxes that overlap
        for surrounding_box in surrounding_boxes:
            for ball in sub_boxes[surrounding_box]:
                # Calculate distance from ball center to subbox boundaries
                min_dist = 0
                for dim in range(3):
                    # Distance to nearest subbox boundary
                    dist_to_min = abs(balls[ball]['loc'][dim] - sub_box[dim] * sub_box_length)
                    dist_to_max = abs(balls[ball]['loc'][dim] - (sub_box[dim] + 1) * sub_box_length)
                    min_dist = max(min_dist, min(dist_to_min, dist_to_max))
                # If ball radius is greater than distance to boundary, it overlaps
                if balls[ball]['radius'] > min_dist:
                    all_balls.append(ball)

        # Place the probe points in the sub box
        probe_points = np.random.uniform(
            [sub_box[0] * sub_box_length, sub_box[1] * sub_box_length, sub_box[2] * sub_box_length],
            [(sub_box[0] + 1) * sub_box_length, (sub_box[1] + 1) * sub_box_length, (sub_box[2] + 1) * sub_box_length],
            (num_probe_points, 3)
        )
        # Create a variable for the count of inside balls
        inside_balls = 0
        # Loop through the probe points and measure against the balls to see if they are inside or outside of a sub_box
        for probe_point in probe_points:
            # Check if the probe point is inside or outside of a ball
            for ball in all_balls:
                # print([sub_box[0] * sub_box_length, sub_box[1] * sub_box_length, sub_box[2] * sub_box_length],
                #       [(sub_box[0] + 1) * sub_box_length, (sub_box[1] + 1) * sub_box_length, (sub_box[2] + 1) * sub_box_length], balls[ball]['loc'])
                # Check if the probe point is inside or outside of a ball
                if np.linalg.norm(probe_point - balls[ball]['loc']) < balls[ball]['radius']:
                    inside_balls += 1
        # Calculate the local density
        local_density = inside_balls / num_probe_points
        # Add the local density to the sub box
        local_densities[sub_box] = local_density
        # Assign the local density to each ball
        for ball in sub_boxes[sub_box]:
            balls[ball]['local_density'] = local_density
        


def calculate_local_polydispersity(balls, what_neighbors='pw'):
    """
    Calculates the local polydispersity of each ball in the balls dictionary.
    """
    # Loop through the balls and calculate the local polydispersity
    for ball in balls.values():
        # Create a list for recording the radii of the neighbors   
        neighbor_radii = []
        # Check if the ball has an 'aw_neighbors' or 'pw_neighbors' key
        if 'aw_neighbors' not in ball or 'pw_neighbors' not in ball:
            continue
        if what_neighbors == 'aw':
            # Loop through the neighbors of the ball
            for neighbor in ball['aw_neighbors']:
                # Record their radii
                neighbor_radii.append(balls[neighbor % 1000]['radius'])
        elif what_neighbors == 'pw':
            # Loop through the neighbors of the ball
            for neighbor in ball['pw_neighbors']:
                # Record their radii
                neighbor_radii.append(balls[neighbor % 1000]['radius'])
        # Calculate the polydispersity
        ball['local_polydispersity'] = np.std(neighbor_radii) / np.mean(neighbor_radii)


def calculate_std_neighbor_distance(balls, what_neighbors='pw'):
    """
    Calculates the standard deviation of the neighbor distance of each ball in the balls dictionary.
    """
    # Loop through the balls and calculate the standard deviation of the neighbor distance
    for ball in balls.values():
        # Check if the ball has an 'aw_neighbors' or 'pw_neighbors' key
        if 'aw_neighbors' not in ball or 'pw_neighbors' not in ball:
            continue
        # Get the neighbors
        if what_neighbors == 'aw':
            neighbors = ball['aw_neighbors']
        elif what_neighbors == 'pw':
            neighbors = ball['pw_neighbors']
        # Calculate the standard deviation of the neighbor distance
        ball['std_neighbor_distance'] = np.std([np.linalg.norm(balls[neighbor]['loc'] - ball['loc']) for neighbor in neighbors])


if __name__ == "__main__":

    # Get the data
    data = get_data(select_file(title='Select the PDB file'), select_file(title='Select the AW log file'), select_file(title='Select the PW log file'))
    # Sort the balls into sub boxes
    sub_boxes = sort_balls(data['box_size'], data['balls'], 5)
    # Calculate the local density
    calculate_local_density(data['balls'], sub_boxes, data['box_size'] / 5)
    # Calculate the local polydispersity
    calculate_local_polydispersity(data['balls'], what_neighbors='pw')
    # Calculate the standard deviation of the neighbor distance
    calculate_std_neighbor_distance(data['balls'], what_neighbors='pw')
    
    # Extract only the balls with required keys
    valid_balls = [ball for ball in data['balls'].values() if all(key in ball for key in ['local_density', 'local_polydispersity', 'sphericity_diff', 'std_neighbor_distance'])]
    
    # Make a surface plot showing the local density vs polydispersity with the height being the volume difference
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Extract the data for plotting
    x = [ball['std_neighbor_distance'] for ball in valid_balls]
    y = [ball['local_polydispersity'] for ball in valid_balls]
    z = [ball['pw_sphericity'] for ball in valid_balls]

    for i in range(len(x)):
        print(x[i], y[i], z[i])
    colors = [ball['radius'] for ball in valid_balls]  # Color by radius
    # for i in range(len(x)):
    #     print(x[i], y[i], z[i], colors[i])
    # Create scatter plot with colorbar
    scatter = ax.scatter(x, y, z, c=colors, cmap='viridis')
    plt.colorbar(scatter, label='Ball Radius')
    
    # Add labels
    ax.set_xlabel('Std Neighbor Distance')
    ax.set_ylabel('Local Polydispersity')
    ax.set_zlabel('Sphericity')
    
    plt.show()

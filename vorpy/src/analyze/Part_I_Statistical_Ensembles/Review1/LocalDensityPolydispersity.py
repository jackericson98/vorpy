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

"""

def sort_balls(box_size, balls, num_splits=10):
    """
    Splits the main ball containers into subsections and calculates the local density and polydispersity of each and compares them to sphericity and volume difference.
    """
    # Sub box length
    sub_box_length = box_size / num_splits
    # Create the sub boxes dictionary which are identified by their tuple of indices
    sub_boxes = {(i, j, k): [] for i in range(num_splits) for j in range(num_splits) for k in range(num_splits)}
    print(balls)
    # Loop through the balls and sort them into the correct subsection for polydispersity
    for ball in balls.values():
        # Find the index of the sub box the ball is in 
        sub_box_index = tuple([np.floor(coord / sub_box_length) for coord in ball['loc']])
        # Add the ball to the correct sub box
        sub_boxes[sub_box_index].append(ball)
    # Return the sub boxes
    return sub_boxes


def calculate_local_density(sub_boxes, sub_box_length, num_probe_points=100):
    """
    Calculates the local density of each sub box.
    """
    # Create a dictionary to store the local densities
    local_densities = {}
    # Loop through the sub boxes and calculate the local density
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
                    dist_to_min = abs(ball['loc'][dim] - sub_box[dim] * sub_box_length)
                    dist_to_max = abs(ball['loc'][dim] - (sub_box[dim] + 1) * sub_box_length)
                    min_dist = max(min_dist, min(dist_to_min, dist_to_max))
                # If ball radius is greater than distance to boundary, it overlaps
                if ball['radius'] > min_dist:
                    all_balls.append(ball)

        # Place the probe points in the sub box
        probe_points = np.random.uniform(0, sub_box_length, (num_probe_points, 3))
        # Create a variable for the count of inside balls
        inside_balls = 0
        # Loop through the probe points and measure against the balls to see if they are inside or outside of a sub_box
        for probe_point in probe_points:
            # Check if the probe point is inside or outside of a ball
            for ball in all_balls:
                # Check if the probe point is inside or outside of a ball
                if np.linalg.norm(probe_point - ball['loc']) < ball['rad']:
                    inside_balls += 1
        # Calculate the local density
        local_density = inside_balls / num_probe_points
        # Add the local density to the sub box
        local_densities[sub_box] = local_density
    # Return the sub boxes
    return local_densities




if __name__ == "__main__":

    # Get the data
    data = get_data(select_file(), select_file(), select_file())
    # Sort the balls into sub boxes
    sub_boxes = sort_balls(data['box_size'], data)
    # Calculate the local density
    local_densities = calculate_local_density(sub_boxes, data['box_size'])

    # Plot the local densities
    plt.scatter(list(local_densities.keys()), list(local_densities.values()))
    plt.show()

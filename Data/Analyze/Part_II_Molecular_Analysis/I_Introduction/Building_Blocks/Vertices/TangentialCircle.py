import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


def plot_circles(centers, radii, colors=None, fig_size=(8, 6), linewidth=4):
    """
    Plot multiple circles on a single figure.

    :param centers: List of tuples, each representing the (x, y) center of a circle.
    :param radii: List of radii corresponding to each center.
    :param colors: List of colors for each circle (optional).
    :param fig_size: Size of the figure (default is (8, 6)).
    """
    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_aspect('equal')  # Ensure the circles aren't distorted

    for center, radius, color in zip(centers, radii, colors):
        circle = Circle(center, radius, color=color, fill=False, linewidth=linewidth)
        ax.add_patch(circle)

    # Set limits based on the circles
    all_x = [x for x, y in centers]
    all_y = [y for x, y in centers]
    buffer = max(radii) * 1.1  # Add a little buffer to ensure circles aren't cut off

    ax.set_xlim(min(all_x) - buffer, max(all_x) + buffer)
    ax.set_ylim(min(all_y) - buffer, max(all_y) + buffer)
    ax.axis('off')  # Turn off the axis
    plt.show()

# Example usage with specific locations and radii
centers = [(2, 3), (2, 7), (6, 7), (3.65, 4.65)]
radii = [1, 1.53, 2, 1.3]
# Assign colors: first three blue, last one red
colors = ['blue', 'blue', 'blue', 'red']

plot_circles(centers, radii, colors)

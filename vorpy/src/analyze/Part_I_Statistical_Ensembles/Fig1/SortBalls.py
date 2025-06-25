import matplotlib.pyplot as plt
import numpy as np
import random

import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

from vorpy.src.calculations.calcs import calc_dist

# Set grid and circle parameters
grid_size = 5  # Grid dimensions
circle_count = 20  # Number of circles to place
overlap_allowance = 0.2  # Allowed overlap factor

# Initialize plot
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xticks(np.arange(0, grid_size + 1))  # Major ticks at grid lines
ax.set_yticks(np.arange(0, grid_size + 1))  # Major ticks at grid lines
ax.set_xticks(np.arange(0.5, grid_size), minor=True)  # Minor ticks between grid lines
ax.set_yticks(np.arange(0.5, grid_size), minor=True)  # Minor ticks between grid lines
ax.set_xticklabels([])  # Turn off major tick labels
ax.set_yticklabels([])  # Turn off major tick labels
ax.set_xticklabels(np.arange(grid_size), minor=True)  # Label minor ticks only
ax.set_yticklabels(np.arange(grid_size), minor=True)  # Label minor ticks only
ax.grid(True, which='major')  # Show grid lines only for major ticks


# Seed for reproducibility
np.random.seed(0)
circles = []
colors = ['r', 'g', 'y', 'b']
# Randomly place circles in the grid ensuring full containment within grid
while len(circles) < circle_count:
    if len(circles) < 4:
        color = colors[len(circles)]
    else:
        color = 'grey'
    # Random coordinates and radius within adjusted boundaries
    radius = np.random.uniform(0.1, 1)
    x, y = random.uniform(radius, grid_size - radius), random.uniform(radius, grid_size - radius)

    # Check overlap with existing circles
    if any(
        calc_dist([x, y], circ[0]) < max(circ[1], radius) * (1 - overlap_allowance) + min(circ[1], radius)
        for circ in circles
    ):
        continue  # Skip if overlap is too high

    # Add circle to plot and list
    ax.add_patch(plt.Circle((x, y), radius, color=color, alpha=0.5))
    ax.scatter(x, y, s=0.1, c='k')  # Mark circle center
    circles.append(([x, y], radius))

# Configure plot display
ax.set_xlim(0, grid_size)
ax.set_ylim(0, grid_size)
ax.set_aspect('equal')

plt.show()

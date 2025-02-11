import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def is_point_inside_circles(x, y, circles):
    """
    Determine if a point is inside any of the given circles.
    """
    for cx, cy, cr in circles:
        if (x - cx)**2 + (y - cy)**2 < cr**2:
            return True
    return False


def can_place_circle(x, y, r, circles, max_overlap=1):
    """
    Check if a new circle can be placed without violating the overlap condition.
    """
    for cx, cy, cr in circles:
        distance = np.sqrt((x - cx)**2 + (y - cy)**2)
        min_distance = r + cr - min(cr, r) * max_overlap
        if distance < min_distance:
            return False
    return True


# Number of circles and points, and square size
num_circles = 60
num_points = 20000
square_size = 10

# Lists to store circle properties
circles = []

# Generate circles ensuring overlap conditions
while len(circles) < num_circles:
    r = np.random.uniform(0.1, 1.0)  # Random radius
    x = np.random.uniform(r, square_size - r)  # X position
    y = np.random.uniform(r, square_size - r)  # Y position
    if can_place_circle(x, y, r, circles):
        circles.append((x, y, r))

# Generate random points
points_x = np.random.uniform(0, square_size, num_points)
points_y = np.random.uniform(0, square_size, num_points)
colors = ['green' if is_point_inside_circles(x, y, circles) else 'red' for x, y in zip(points_x, points_y)]
num_in = len([_ for _ in colors if _ == 'green'])

# Create a plot
fig, axes = plt.subplots(1, 2, figsize=(8, 4))  # Adjust the figure size
ax, ax1 = axes
# Scatter points
ax.scatter(points_x, points_y, c=colors, s=0.3, alpha=0.5)

# Add circles to the plot
for x, y, r in circles:
    circle = plt.Circle((x, y), r, color='black', fill=False, linewidth=1)
    ax.add_artist(circle)

# Set the circle plot stuff
ax.set_xlim([0, square_size])
ax.set_ylim([0, square_size])
ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect('equal')
# Make all spines (borders) thicker
for spine in ax.spines.values():
    spine.set_linewidth(1)  # Set the thickness of the border here

# Calculate the circle area
circle_area = sum(np.pi * (r ** 2) for x, y, r in circles)

# Text with box for better readability
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# Calculations
square_area_val = square_size ** 2
initial_density = circle_area / square_area_val
real_density = num_in / num_points

# Formatting string
textstr = (
    "{:<20s}{:10s}{:<20s}\n"
    "{:<20s}{:10s}{:<20s}\n"
    "{:<20s}{:10s}{:<20s}\n"
    "{:<20s}{:10s}{:<20s}\n"
    "{:<20s}{:10s}{:<20s}\n"
    "{:<20s}{:10s}{:<20s}"
).format(
    'Square Area', '=', str(square_area_val),
    'Total Circle Area', '=', str(round(circle_area, 3)),
    'Initial Density', '=', str(round(initial_density, 3)),
    'Points Inside', '=', str(num_in),
    'Points Outside', '=', str(num_points - num_in),
    'Output Density', '=', str(round(real_density, 3))
)

labels = ("{:<20s}\n{:<20s}\n{:<20s}\n{:<20s}\n{:<20s}\n{:<20s}").format(
    'Square Area', 'Total Circle Area', 'Initial Density', 'Points Inside', 'Points Outside', 'Output Density')
equals = "=\n=\n=\n=\n=\n="
values = "{:<20s}\n{:<20s}\n{:<20s}\n{:<20s}\n{:<20s}\n{:<20s}".format(
    str(square_area_val), str(round(circle_area, 3)), str(round(initial_density, 3)), str(num_in),
    str(num_points - num_in), str(round(real_density, 3)))

# mpl.rcParams['font.size'] = 12
# mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['font.weight'] = 'bold'

# Position the text outside the plotting area
ax1.text(-0.05, 0.5, labels, fontsize=15, va='center')
ax1.text(0.65, 0.5, equals, fontsize=15, va='center')
ax1.text(0.8, 0.5, values, fontsize=15, va='center')
ax1.axis('off')


# Show the plot
plt.show()

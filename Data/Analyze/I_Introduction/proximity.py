import numpy as np
import matplotlib.pyplot as plt

# Define circle parameters
radius1 = 1
radius2 = 2
distance_between_circles = 4

num_points = 1000

# Sort the radii
if radius1 > radius2:
    radius1, radius2 = radius2, radius1
# Get the areas
s_ar = np.pi * radius1 ** 2  # Small Area
l_ar = np.pi * radius2 ** 2  # Large Area
f_ar = np.pi * (2 * radius2) ** 2  # Double the radius of the large area or the total area for both circles reach

# Get the relative ratios

half_points = num_points // 2
small_in_num, small_out_num = int(half_points * s_ar / f_ar), int(half_points * (f_ar - s_ar) / f_ar)
large_in_num, large_out_num = int(half_points * l_ar / f_ar), int(half_points * (f_ar - l_ar) / f_ar)


# Generate points around the circles
theta = np.linspace(0, 2*np.pi, num_points)
x_circle1 = radius1 * np.cos(theta)
y_circle1 = radius1 * np.sin(theta)
x_circle2 = radius2 * np.cos(theta) + distance_between_circles
y_circle2 = radius2 * np.sin(theta)

# Generate random points both inside and outside the circles
np.random.seed()  # for reproducibility
angle_rand_inside = np.random.uniform(0, 2 * np.pi, small_in_num)
radius_rand_inside = np.random.uniform(0, radius1, small_in_num)
x1_rand_inside = np.concatenate([radius_rand_inside * np.cos(angle_rand_inside),
                                distance_between_circles + radius_rand_inside * np.cos(angle_rand_inside)])
y1_rand_inside = np.concatenate([radius_rand_inside * np.sin(angle_rand_inside),
                                radius_rand_inside * np.sin(angle_rand_inside)])

angle_rand_outside = np.random.uniform(0, 2 * np.pi, small_out_num)
radius_rand_outside = np.random.uniform(radius1, 2 * radius2, small_out_num)
x1_rand_outside = np.concatenate([radius_rand_outside * np.cos(angle_rand_outside),
                                 distance_between_circles + radius_rand_outside * np.cos(angle_rand_outside)])
y1_rand_outside = np.concatenate([radius_rand_outside * np.sin(angle_rand_outside),
                                 radius_rand_outside * np.sin(angle_rand_outside)])

# Generate random points both inside and outside the circles
np.random.seed()  # for reproducibility
angle_rand_inside = np.random.uniform(0, 2 * np.pi, large_in_num)
radius_rand_inside = np.random.uniform(0, radius1, large_in_num)
x2_rand_inside = np.concatenate([radius_rand_inside * np.cos(angle_rand_inside),
                                distance_between_circles + radius_rand_inside * np.cos(angle_rand_inside)])
y2_rand_inside = np.concatenate([radius_rand_inside * np.sin(angle_rand_inside),
                                radius_rand_inside * np.sin(angle_rand_inside)])

angle_rand_outside = np.random.uniform(0, 2 * np.pi, large_out_num)
radius_rand_outside = np.random.uniform(radius1, 2 * radius2, large_out_num)
x2_rand_outside = np.concatenate([radius_rand_outside * np.cos(angle_rand_outside),
                                 distance_between_circles + radius_rand_outside * np.cos(angle_rand_outside)])
y2_rand_outside = np.concatenate([radius_rand_outside * np.sin(angle_rand_outside),
                                 radius_rand_outside * np.sin(angle_rand_outside)])

# Create figure and axis
fig, ax = plt.subplots()

# Plot circles
ax.plot(x_circle1, y_circle1, 'r')
ax.plot(x_circle2, y_circle2, 'b')

# Plot random points inside and outside the circles
# ax.scatter(x_rand_inside, y_rand_inside, color='blue', label='Points Closer to Blue Circle')
# ax.scatter(x_rand_outside, y_rand_outside, color='red', label='Points Closer to Red Circle')

# Color the points based on the closest circle
for x, y in zip(x1_rand_inside, y1_rand_inside):
    distance_to_circle1 = np.sqrt((x - x_circle1)**2 + (y - y_circle1)**2)
    distance_to_circle2 = np.sqrt((x - x_circle2)**2 + (y - y_circle2)**2)
    if np.min(distance_to_circle1) < np.min(distance_to_circle2):
        ax.scatter(x, y, color='r', marker='.')
    else:
        ax.scatter(x, y, color='b', marker='.')

for x, y in zip(x2_rand_outside, y2_rand_outside):
    distance_to_circle1 = np.sqrt((x - x_circle1)**2 + (y - y_circle1)**2)
    distance_to_circle2 = np.sqrt((x - x_circle2)**2 + (y - y_circle2)**2)
    if np.min(distance_to_circle1) < np.min(distance_to_circle2):
        ax.scatter(x, y, color='r', marker='.')
    else:
        ax.scatter(x, y, color='b', marker='.')

# Set aspect ratio to equal and add legend
ax.set_aspect('equal', 'box')
# ax.legend()
plt.xticks([])
plt.yticks([])
plt.axis('off')

# Show plot
plt.show()


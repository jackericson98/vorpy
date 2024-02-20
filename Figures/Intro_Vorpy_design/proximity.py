import numpy as np
import matplotlib.pyplot as plt

# Define circle parameters
radius = 1
distance_between_circles = 3

# Generate points around the circles
num_points = 100
theta = np.linspace(0, 2*np.pi, num_points)
x_circle1 = radius * np.cos(theta)
y_circle1 = radius * np.sin(theta)
x_circle2 = radius * np.cos(theta) + distance_between_circles
y_circle2 = radius * np.sin(theta)

# Generate random points both inside and outside the circles
np.random.seed(0)  # for reproducibility
angle_rand_inside = np.random.uniform(0, 2 * np.pi, num_points // 2)
radius_rand_inside = np.random.uniform(0, radius, num_points // 2)
x_rand_inside = np.concatenate([radius_rand_inside * np.cos(angle_rand_inside),
                                distance_between_circles + radius_rand_inside * np.cos(angle_rand_inside)])
y_rand_inside = np.concatenate([radius_rand_inside * np.sin(angle_rand_inside),
                                radius_rand_inside * np.sin(angle_rand_inside)])

angle_rand_outside = np.random.uniform(0, 2 * np.pi, num_points // 2)
radius_rand_outside = np.random.uniform(radius, 2 * radius, num_points // 2)
x_rand_outside = np.concatenate([radius_rand_outside * np.cos(angle_rand_outside),
                                 distance_between_circles + radius_rand_outside * np.cos(angle_rand_outside)])
y_rand_outside = np.concatenate([radius_rand_outside * np.sin(angle_rand_outside),
                                 radius_rand_outside * np.sin(angle_rand_outside)])

# Create figure and axis
fig, ax = plt.subplots()

# Plot circles
ax.plot(x_circle1, y_circle1, 'r')
ax.plot(x_circle2, y_circle2, 'b')

# Plot random points inside and outside the circles
ax.scatter(x_rand_inside, y_rand_inside, color='blue', label='Points Closer to Blue Circle')
ax.scatter(x_rand_outside, y_rand_outside, color='red', label='Points Closer to Red Circle')

# Color the points based on the closest circle
for x, y in zip(x_rand_inside, y_rand_inside):
    distance_to_circle1 = np.sqrt((x - x_circle1)**2 + (y - y_circle1)**2)
    distance_to_circle2 = np.sqrt((x - x_circle2)**2 + (y - y_circle2)**2)
    if np.min(distance_to_circle1) < np.min(distance_to_circle2):
        ax.scatter(x, y, color='r')
    else:
        ax.scatter(x, y, color='b')

for x, y in zip(x_rand_outside, y_rand_outside):
    distance_to_circle1 = np.sqrt((x - x_circle1)**2 + (y - y_circle1)**2)
    distance_to_circle2 = np.sqrt((x - x_circle2)**2 + (y - y_circle2)**2)
    if np.min(distance_to_circle1) < np.min(distance_to_circle2):
        ax.scatter(x, y, color='r')
    else:
        ax.scatter(x, y, color='b')

# Set aspect ratio to equal and add legend
ax.set_aspect('equal', 'box')
# ax.legend()
plt.xticks([])
plt.yticks([])

# Show plot
plt.show()

